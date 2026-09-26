// $Id$
//==============================================================================
//!
//! \file MultigridTransfer.C
//!
//! \date Sep 23 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Grid transfer operators for geometric multigrid.
//!
//==============================================================================

#include "MultigridTransfer.h"

#include "ASMbase.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "DomainDecomposition.h"
#include "GaussQuadrature.h"
#include "IFEM.h"
#include "LogStream.h"
#include "MatVec.h"
#include "ProcessAdm.h"
#include "SAM.h"
#include "SIMbase.h"
#include "SparseMatrix.h"
#include "Utilities.h"

#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

#ifdef HAS_LRSPLINE
#include "LR/ASMu2D.h"
#include "LR/ASMu3D.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/Element.h"
#include "LRSpline/LRSpline.h"
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/MeshRectangle.h"
#include "LRSpline/Meshline.h"
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>


namespace // anonymous namespace for local helpers
{
  /*!
    \brief The two spline bases a transfer operator is built between.

    \details The coarse and the fine basis are always of the same kind, both
    tensor product or both locally refined, and what has to be done with them
    differs more between the two kinds than it agrees. A tensor basis is a
    product of univariate ones and both operators factor along with it, while
    a locally refined basis has no such structure and has to be taken function
    by function. What they have in common is the two questions asked of them
    here, and holding the pair rather than one basis at a time is what lets
    the kind be settled once, where the pair is made.
  */

  class BasisPair
  {
  public:
    //! \brief Empty destructor.
    virtual ~BasisPair() = default;

    //! \brief Returns the number of parameter directions.
    virtual int nVariate() const = 0;
    //! \brief Returns whether the two bases are of the same polynomial order.
    virtual bool sameOrder() const = 0;
    //! \brief Returns the number of functions in the coarse basis.
    virtual size_t nCoarse() const = 0;
    //! \brief Returns the number of functions in the fine basis.
    virtual size_t nFine() const = 0;

    //! \brief Expresses each coarse basis function in the fine basis.
    //! \param[out] rows For each fine function, its coefficient against each
    //! coarse function it takes part in
    //! \return \e false if the coarse space is not contained in the fine one
    virtual bool insertKnots(std::vector<std::map<int,Real>>& rows) const = 0;

    /*!
      \brief Integrates the two bases against each other and the fine one
      against itself.
      \param[in] fPch The patch the fine basis belongs to
      \param[in] mine Elements this process integrates, all of them if empty
      \param mass Mass matrix of the fine basis of the patch
      \param B The two bases of the patch against each other
      \return \e false if the two meshes do not cover each other

      \details Both come out as sparse as the bases they belong to, which the
      projection they define is not: it is the one solved out of the other,
      and solving it is left to whoever applies it, so that the dense matrix
      the two multiply to is never formed.

      The integrals are over the elements of the fine mesh, each of which lies
      within one element of the coarse one, in the inner product of the
      parameter domain.
    */
    virtual bool integrate(const ASMbase& fPch, const std::set<int>& mine,
                           SparseMatrix& mass, SparseMatrix& B) const = 0;
  };


  //! \brief Values of a univariate basis, and the function the first is of.
  struct Values
  {
    int       first = 0; //!< Index of the function the first value belongs to
    RealArray val;       //!< One value per function from \a first onwards
  };


  /*!
    \brief A tensor product spline basis, as the transfer operators see it.

    \details What they ask of it is the univariate bases it is a product of
    and, where the patch is rational, the weight each of its functions
    carries. The functions are numbered with the first parameter direction
    running fastest, which is the order the patch numbers its nodes in.
  */

  struct TensorBasis
  {
    std::vector<const Go::BsplineBasis*> dirs; //!< The univariate bases
    RealArray weight; //!< Weight of each function, empty if not rational

    //! \brief Returns the number of functions in a direction.
    int size(int dir) const { return dirs[dir]->numCoefs(); }
    //! \brief Returns the polynomial order in a direction.
    int order(int dir) const { return dirs[dir]->order(); }
    //! \brief Returns the knot vector of a direction.
    RealArray knots(int dir) const
    { return RealArray(dirs[dir]->begin(),dirs[dir]->end()); }
    //! \brief Returns one knot of a direction.
    Real knot(int dir, int i) const { return *(dirs[dir]->begin()+i); }

    //! \brief Returns the total number of functions.
    size_t nFunctions() const
    {
      size_t n = 1;
      for (const Go::BsplineBasis* b : dirs) n *= b->numCoefs();
      return n;
    }

    //! \brief Returns the weight of a function, one where it carries none.
    Real w(int i) const { return weight.empty() ? Real(1) : weight[i]; }
  };


  //! \brief Picks the univariate bases and the weights out of a spline patch.
  //! \param[in] geo The spline surface or volume
  //! \param[in] nsd Number of parameter directions it has
  template<class Spline>
  TensorBasis tensorBasis (const Spline* geo, int nsd)
  {
    TensorBasis basis;
    for (int d = 0; d < nsd; d++)
      basis.dirs.push_back(&geo->basis(d));

    // A rational patch stores its coefficients scaled by their weights, with
    // the weight last of each. The weights are what the transfer needs: the
    // functions of a rational basis are the weighted B-splines over the one
    // weight function, which refinement leaves as it is, so it cancels
    // between the two levels and what remains is the two sets of weights.
    if (geo->rational())
    {
      const int nsdg = geo->dimension();
      std::vector<double>::const_iterator c = geo->rcoefs_begin();
      basis.weight.reserve(basis.nFunctions());
      for (size_t i = 0; i < basis.nFunctions(); i++)
        basis.weight.push_back(c[(nsdg+1)*i + nsdg]);
    }

    return basis;
  }


  /*!
    \brief Multiplies univariate values out into the tensor basis.
    \param[in] basis The basis the functions belong to
    \param[in] dir Values of each univariate basis
    \param[out] terms The functions reached and the value of each

    \details The functions are numbered with the first direction running
    fastest, and each is weighted as the basis weights it, so that what comes
    out are the values of the functions the patch carries a coefficient for.
  */

  void tensorValues (const TensorBasis& basis,
                     const std::vector<const Values*>& dir,
                     std::vector<std::pair<int,Real>>& terms)
  {
    terms.assign(1,{0,Real(1)});

    std::vector<std::pair<int,Real>> next;
    int stride = 1;
    for (size_t d = 0; d < dir.size(); d++)
    {
      next.clear();
      next.reserve(terms.size()*dir[d]->val.size());
      for (const std::pair<int,Real>& t : terms)
        for (size_t i = 0; i < dir[d]->val.size(); i++)
          next.emplace_back(t.first + (dir[d]->first+i)*stride,
                            t.second*dir[d]->val[i]);
      terms.swap(next);
      stride *= basis.size(d);
    }

    for (std::pair<int,Real>& t : terms)
      t.second *= basis.w(t.first);
  }


  /*!
    \brief Expresses a univariate B-spline basis in a refinement of itself.
    \param[in] s Knot vector of the coarse basis
    \param[in] t Knot vector of the fine basis
    \param[in] k Order of both bases
    \param[out] rows For each fine function, the coarse functions it takes
    part in and with which coefficient
    \return \e false if the coarse knots are not among the fine ones

    \details These are the discrete B-splines of the Oslo algorithm. The
    recurrence computing them mirrors the one defining the B-splines
    themselves, and its coefficients are the ones knot insertion gives: over
    a refinement each lies between zero and one, so nothing is subtracted
    anywhere and no precision is lost however far apart the two meshes are.

    The coefficients of a fine function against the coarse ones sum to one
    whatever the two knot vectors are, so the partition of unity the operator
    is checked against downstream says nothing about nestedness. That is why
    it is checked here instead, where the knots are.
  */

  bool oslo (const RealArray& s, const RealArray& t, int k,
             std::vector<Values>& rows)
  {
    const Real eps = Real(1.0e-12);

    // Every knot of the coarse vector has to be among the fine ones, with at
    // least the multiplicity the coarse vector gives it, or the coarse space
    // is not contained in the fine one and there is nothing to compute.
    for (size_t i = 0, j = 0; i < s.size(); i++, j++)
    {
      while (j < t.size() && t[j] < s[i]-eps)
        ++j;
      if (j >= t.size() || t[j] > s[i]+eps)
        return false;
    }

    const int nc = static_cast<int>(s.size()) - k;
    const int nf = static_cast<int>(t.size()) - k;
    if (nc < k || nf < k)
      return false;

    rows.assign(nf,Values());
    for (int j = 0; j < nf; j++)
    {
      // The last coarse knot interval starting at or before this function
      int mu = k-1;
      while (mu+1 < nc && s[mu+1] <= t[j]+eps)
        ++mu;

      RealArray a(1,Real(1));
      for (int r = 1; r < k; r++)
      {
        RealArray b(r+1,Real(0));
        for (int l = 0; l <= r; l++)
        {
          const int i = mu-r+l;
          if (l > 0)
            if (Real d = s[i+r] - s[i]; d > eps)
              b[l] += (t[j+r] - s[i])/d * a[l-1];
          if (l < r)
            if (Real d = s[i+r+1] - s[i+1]; d > eps)
              b[l] += (s[i+r+1] - t[j+r])/d * a[l];
        }
        a.swap(b);
      }

      rows[j].first = mu-k+1;
      rows[j].val.swap(a);
    }

    return true;
  }


  /*!
    \brief The transfer between two tensor product spline bases.

    \details Both operators factor along the parameter directions, which is
    what makes this cheaper than the same thing on a locally refined mesh
    rather than merely a special case of it. The change of basis between two
    nested spaces is the tensor product of the univariate ones, each of which
    the Oslo algorithm gives directly from the two knot vectors.
  */

  class TensorPair : public BasisPair
  {
  public:
    //! \brief The constructor takes over the two bases.
    TensorPair(TensorBasis&& c, TensorBasis&& f)
      : coarse(std::move(c)), fine(std::move(f)) {}

    //! \brief Returns the number of parameter directions.
    int nVariate() const override { return fine.dirs.size(); }
    //! \brief Returns the number of functions in the coarse basis.
    size_t nCoarse() const override { return coarse.nFunctions(); }
    //! \brief Returns the number of functions in the fine basis.
    size_t nFine() const override { return fine.nFunctions(); }

    //! \brief Returns whether the two bases are of the same polynomial order.
    bool sameOrder() const override
    {
      for (int d = 0; d < this->nVariate(); d++)
        if (coarse.order(d) != fine.order(d))
          return false;

      return true;
    }

    //! \brief Expresses each coarse basis function in the fine basis.
    bool insertKnots(std::vector<std::map<int,Real>>& rows) const override;
    //! \brief Integrates the two bases against each other and itself.
    bool integrate(const ASMbase& fPch, const std::set<int>& mine,
                   SparseMatrix& mass, SparseMatrix& B) const override;

  private:
    //! \brief Evaluates the univariate bases of one of the two in a point.
    //! \param[in] basis The basis to evaluate
    //! \param[in] X The parameter point
    //! \param[out] val Values of each univariate basis
    static void evaluate(const TensorBasis& basis, const RealArray& X,
                         std::vector<Values>& val)
    {
      for (size_t d = 0; d < basis.dirs.size(); d++)
      {
        val[d].val.resize(basis.order(d));
        val[d].first = basis.dirs[d]->knotInterval(X[d]) - basis.order(d) + 1;
        basis.dirs[d]->computeBasisValues(X[d],val[d].val.data());
      }
    }

    TensorBasis coarse; //!< The coarse basis
    TensorBasis fine;   //!< The fine basis
  };


  bool TensorPair::insertKnots (std::vector<std::map<int,Real>>& rows) const
  {
    const int nsd = this->nVariate();

    // The change of basis of each direction on its own
    std::vector<std::vector<Values>> along(nsd);
    for (int d = 0; d < nsd; d++)
      if (!oslo(coarse.knots(d),fine.knots(d),fine.order(d),along[d]))
        return false;

    // Coefficients below this are roundoff rather than structure
    const Real dropTol = Real(1.0e-12);

    std::vector<const Values*> dir(nsd);
    std::vector<std::pair<int,Real>> terms;
    IntVec jd(nsd,0);
    for (size_t j = 0; j < fine.nFunctions(); j++)
    {
      for (int d = 0; d < nsd; d++)
        dir[d] = &along[d][jd[d]];

      // The coefficients come out scaled by the weights of the coarse basis,
      // and the fine functions they belong to carry their own.
      tensorValues(coarse,dir,terms);
      for (const std::pair<int,Real>& t : terms)
        if (fabs(t.second) > dropTol)
          rows[j][t.first] += t.second/fine.w(j);

      for (int d = 0; d < nsd; d++)
        if (++jd[d] < fine.size(d) || d == nsd-1)
          break;
        else
          jd[d] = 0;
    }

    return true;
  }


  bool TensorPair::integrate (const ASMbase& fPch, const std::set<int>& mine,
                              SparseMatrix& mass, SparseMatrix& B) const
  {
    const int nsd = this->nVariate();

    // A rule which integrates the product of the two bases exactly
    IntVec nG(nsd);
    int nGP = 1;
    for (int d = 0; d < nsd; d++)
    {
      nG[d] = (fine.order(d) + coarse.order(d))/2 + 1;
      if (!GaussQuadrature::getCoord(nG[d]))
      {
        std::cerr <<" *** MG::prolongation: No Gauss rule with "<< nG[d]
                  <<" points, needed to integrate an order "<< fine.order(d)
                  <<" basis against an order "<< coarse.order(d) <<" one."
                  << std::endl;
        return false;
      }
      nGP *= nG[d];
    }

    // The elements are the knot spans, counted the way the patch counts them
    // so that the ones this process was given can be told apart. Those of
    // zero measure are counted as well and integrate to nothing.
    IntVec nElm(nsd);
    int nel = 1;
    for (int d = 0; d < nsd; d++)
      nel *= nElm[d] = fine.size(d) - fine.order(d) + 1;

    std::vector<Values> fVal(nsd), cVal(nsd);
    std::vector<const Values*> fDir(nsd), cDir(nsd);
    for (int d = 0; d < nsd; d++)
    {
      fDir[d] = &fVal[d];
      cDir[d] = &cVal[d];
    }

    std::vector<std::pair<int,Real>> fTerm, cTerm;
    RealArray X(nsd), x0(nsd), x1(nsd);
    IntVec ed(nsd,0);
    for (int iel = 0; iel < nel; iel++)
    {
      Real vol = Real(1);
      for (int d = 0; d < nsd; d++)
      {
        const int mu = ed[d] + fine.order(d) - 1;
        x0[d] = fine.knot(d,mu);
        x1[d] = fine.knot(d,mu+1);
        vol *= Real(0.5)*(x1[d] - x0[d]);
      }

      // A partitioned mesh has each process integrate the elements it was
      // given, and what they leave is added to what the others do.
      if (vol > Real(0) &&
          (mine.empty() || mine.find(fPch.getElmID(1+iel)) != mine.end()))
      {
        IntVec ig(nsd,0);
        for (int ip = 0; ip < nGP; ip++)
        {
          Real w = vol;
          for (int d = 0; d < nsd; d++)
          {
            const double* xg = GaussQuadrature::getCoord(nG[d]);
            const double* wg = GaussQuadrature::getWeight(nG[d]);
            X[d] = Real(0.5)*((x1[d]-x0[d])*xg[ig[d]] + x1[d] + x0[d]);
            w *= wg[ig[d]];
          }

          this->evaluate(fine,X,fVal);
          this->evaluate(coarse,X,cVal);
          tensorValues(fine,fDir,fTerm);
          tensorValues(coarse,cDir,cTerm);

          for (const std::pair<int,Real>& fi : fTerm)
          {
            const Real Ni = w*fi.second;
            for (const std::pair<int,Real>& fj : fTerm)
              mass(1+fi.first,1+fj.first) += Ni*fj.second;
            for (const std::pair<int,Real>& cj : cTerm)
              B(1+fi.first,1+cj.first) += Ni*cj.second;
          }

          for (int d = 0; d < nsd; d++)
            if (++ig[d] < nG[d] || d == nsd-1)
              break;
            else
              ig[d] = 0;
        }
      }

      for (int d = 0; d < nsd; d++)
        if (++ed[d] < nElm[d] || d == nsd-1)
          break;
        else
          ed[d] = 0;
    }

    return true;
  }


#ifdef HAS_LRSPLINE
  //! \brief A knot line of a mesh, over the part of the mesh it covers.
  struct KnotLine
  {
    int       dir  = 0;   //!< Direction the knot is inserted in
    Real      par  = 0.0; //!< Value of the knot
    int       mult = 1;   //!< Multiplicity the mesh gives it
    RealArray start;      //!< Lower corner of what it covers, other directions
    RealArray stop;       //!< Upper corner of what it covers, other directions
  };


  //! \brief Collects the knot lines of a mesh.
  //! \param[in] basis The mesh to collect the lines of
  //! \param[out] lines The knot lines of that mesh
  //! \return \e false if the mesh is of a kind not covered here
  bool knotLines (const LR::LRSpline* basis, std::vector<KnotLine>& lines)
  {
    if (const LR::LRSplineSurface* srf =
        dynamic_cast<const LR::LRSplineSurface*>(basis); srf)
    {
      for (const LR::Meshline* m : srf->getAllMeshlines())
      {
        KnotLine& line = lines.emplace_back();
        // A line spanning u is a line of constant v, so it is a knot in v
        line.dir   = m->is_spanning_u() ? 1 : 0;
        line.par   = m->const_par_;
        line.mult  = m->multiplicity();
        line.start = { m->start_ };
        line.stop  = { m->stop_ };
      }
      return true;
    }

    if (const LR::LRSplineVolume* vol =
        dynamic_cast<const LR::LRSplineVolume*>(basis); vol)
    {
      for (const LR::MeshRectangle* m : vol->getAllMeshRectangles())
      {
        KnotLine& line = lines.emplace_back();
        line.dir  = m->constDirection();
        line.par  = m->constParameter();
        line.mult = m->multiplicity();
        for (int d = 0; d < 3; d++)
          if (d != line.dir)
          {
            line.start.push_back(m->start_[d]);
            line.stop.push_back(m->stop_[d]);
          }
      }
      return true;
    }

    return false;
  }


  //! \brief Checks whether a knot line splits the support of a B-spline.
  //! \param[in] line The knot line
  //! \param[in] knots Local knot vectors of the B-spline
  bool splits (const KnotLine& line, const std::vector<RealArray>& knots)
  {
    const Real eps = Real(1.0e-12);
    const RealArray& x = knots[line.dir];
    if (line.par <= x.front()+eps || line.par >= x.back()-eps)
      return false; // the knot lies outside the support

    int have = 0;
    for (Real k : x)
      if (fabs(k-line.par) < eps)
        ++have;
    if (have >= line.mult)
      return false; // the knot is there as many times as the mesh has it

    // The line has to cover the support in the directions it spans, or it
    // stops short of the function and leaves it whole.
    for (size_t d = 0, k = 0; d < knots.size(); d++)
      if (static_cast<int>(d) != line.dir)
      {
        if (knots[d].front() < line.start[k]-eps ||
            knots[d].back()  > line.stop[k]+eps)
          return false;
        ++k;
      }

    return true;
  }


  //! \brief Splits a B-spline in two by inserting a knot in one direction.
  //! \param[in] knots Local knot vectors of the B-spline
  //! \param[in] line The knot line to insert
  //! \param[out] lo The B-spline covering the lower part, and its coefficient
  //! \param[out] hi The B-spline covering the upper part, and its coefficient
  //!
  //! \details Both coefficients lie between zero and one and they sum to one,
  //! which is what makes this stable where solving for them is not: nothing
  //! cancels, so no precision is lost however many times it is done.
  void splitBspline (const std::vector<RealArray>& knots, const KnotLine& line,
                     std::pair<std::vector<RealArray>,Real>& lo,
                     std::pair<std::vector<RealArray>,Real>& hi)
  {
    const RealArray& x = knots[line.dir];
    const size_t n = x.size()-1; // index of the last knot, the order

    RealArray a(x.begin(),x.end()-1);
    a.insert(std::upper_bound(a.begin(),a.end(),line.par),line.par);
    RealArray b(x.begin()+1,x.end());
    b.insert(std::upper_bound(b.begin(),b.end(),line.par),line.par);

    lo.first = knots; lo.first[line.dir] = a;
    hi.first = knots; hi.first[line.dir] = b;

    lo.second = line.par >= x[n-1] ? Real(1)
                                   : (line.par - x[0])/(x[n-1] - x[0]);
    hi.second = line.par <= x[1]   ? Real(1)
                                   : (x[n] - line.par)/(x[n] - x[1]);
  }


  //! \brief Evaluates a basis function in a parameter point.
  //! \param[in] f The basis function
  //! \param[in] X The parameter point
  double evalBasis (const LR::Basisfunction* f, const RealArray& X)
  {
    return X.size() == 2 ? f->evaluate(X[0],X[1]) : f->evaluate(X[0],X[1],X[2]);
  }


  //! \brief Keys a B-spline by its local knot vectors.
  RealArray knotKey (const std::vector<RealArray>& knots)
  {
    RealArray key;
    for (const RealArray& k : knots)
      key.insert(key.end(),k.begin(),k.end());
    return key;
  }


  /*!
    \brief The transfer between two locally refined spline bases.

    \details Neither operator factors here, so both are taken function by
    function and element by element over meshes which have no structure to
    lean on.
  */

  class LRPair : public BasisPair
  {
  public:
    //! \brief The constructor holds on to the two bases.
    LRPair(const LR::LRSpline* c, const LR::LRSpline* f) : cB(c), fB(f) {}

    //! \brief Returns the number of parameter directions.
    int nVariate() const override { return fB->nVariate(); }
    //! \brief Returns the number of functions in the coarse basis.
    size_t nCoarse() const override { return cB->nBasisFunctions(); }
    //! \brief Returns the number of functions in the fine basis.
    size_t nFine() const override { return fB->nBasisFunctions(); }

    //! \brief Returns whether the two bases are of the same polynomial order.
    bool sameOrder() const override
    {
      for (int d = 0; d < this->nVariate(); d++)
        if (cB->order(d) != fB->order(d))
          return false;

      return true;
    }

    //! \brief Expresses each coarse basis function in the fine basis.
    bool insertKnots(std::vector<std::map<int,Real>>& rows) const override;
    //! \brief Integrates the two bases against each other and itself.
    bool integrate(const ASMbase& fPch, const std::set<int>& mine,
                   SparseMatrix& mass, SparseMatrix& B) const override;

  private:
    const LR::LRSpline* cB; //!< The coarse basis
    const LR::LRSpline* fB; //!< The fine basis
  };


  /*!
    A coarse function is a B-spline on its own local knot vectors, and the
    fine mesh cuts its support with knot lines the coarse mesh does not have.
    Inserting one of those splits it in two B-splines whose knot vectors have
    the knot, weighted so that the two together are what was split. Repeating
    until every piece is a function of the fine basis expresses the coarse
    function in that basis.

    The weights are the ones knot insertion gives, each between zero and one
    and summing to one over a split. Nothing is subtracted anywhere, so no
    precision is lost however deep the refinement goes; this is what the
    alternative of solving for the coefficients on each element cannot do,
    the systems there conditioning like the Bernstein basis of the order and
    passing that on from element to element.
  */

  bool LRPair::insertKnots (std::vector<std::map<int,Real>>& rows) const
  {
    std::vector<KnotLine> lines;
    if (!knotLines(fB,lines))
      return false;

    // The functions of the fine basis, looked up by their knot vectors
    std::map<RealArray,const LR::Basisfunction*> fine;
    for (const LR::Basisfunction* f : fB->getAllBasisfunctions())
    {
      std::vector<RealArray> knots(f->nVariate());
      for (int d = 0; d < f->nVariate(); d++)
        knots[d] = (*f)[d];
      fine[knotKey(knots)] = f;
    }

    // Coefficients below this are roundoff rather than structure
    const Real dropTol = Real(1.0e-12);

    typedef std::pair<std::vector<RealArray>,Real> Term;
    std::vector<Term> todo, split(2);
    for (const LR::Basisfunction* c : cB->getAllBasisfunctions())
    {
      todo.clear();
      Term& first = todo.emplace_back();
      first.first.resize(c->nVariate());
      for (int d = 0; d < c->nVariate(); d++)
        first.first[d] = (*c)[d];
      first.second = c->w();

      while (!todo.empty())
      {
        const Term term = todo.back();
        todo.pop_back();
        if (fabs(term.second) < dropTol)
          continue;

        std::map<RealArray,const LR::Basisfunction*>::const_iterator it =
          fine.find(knotKey(term.first));
        if (it != fine.end())
        {
          // The weights scale the B-splines into the functions of the basis
          rows[it->second->getId()][c->getId()] += term.second/it->second->w();
          continue;
        }

        const std::vector<KnotLine>::const_iterator line =
          std::find_if(lines.begin(),lines.end(),
                       [&term](const KnotLine& l)
                       { return splits(l,term.first); });
        if (line == lines.end())
          return false; // the two meshes are not nested

        splitBspline(term.first,*line,split[0],split[1]);
        for (Term& half : split)
        {
          half.second *= term.second;
          todo.push_back(half);
        }
      }
    }

    return true;
  }


  bool LRPair::integrate (const ASMbase& fPch, const std::set<int>& mine,
                          SparseMatrix& mass, SparseMatrix& B) const
  {
    const int nsd = fB->nVariate();

    // A rule which integrates the product of the two bases exactly
    IntVec nG(nsd);
    int nGP = 1;
    for (int d = 0; d < nsd; d++)
    {
      nG[d] = (fB->order(d) + cB->order(d))/2 + 1;
      if (!GaussQuadrature::getCoord(nG[d]))
      {
        std::cerr <<" *** MG::prolongation: No Gauss rule with "<< nG[d]
                  <<" points, needed to integrate an order "<< fB->order(d)
                  <<" basis against an order "<< cB->order(d) <<" one."
                  << std::endl;
        return false;
      }
      nGP *= nG[d];
    }

    RealArray X(nsd);
    for (int iel = 0; iel < fB->nElements(); iel++)
    {
      // A partitioned mesh has each process integrate the elements it was
      // given, and what they leave is added to what the others do.
      if (!mine.empty() && mine.find(fPch.getElmID(1+iel)) == mine.end())
        continue;

      const LR::Element* fEl = fB->getElement(iel);
      for (int d = 0; d < nsd; d++)
        X[d] = 0.5*(fEl->getParmin(d) + fEl->getParmax(d));

      const int cel = cB->getElementContaining(X);
      if (cel < 0)
      {
        std::cerr <<" *** MG::prolongation: No coarse element contains the"
                  <<" midpoint of fine element "<< 1+iel <<"."<< std::endl;
        return false;
      }
      const LR::Element* cEl = cB->getElement(cel);

      Real vol = Real(1);
      for (int d = 0; d < nsd; d++)
        vol *= Real(0.5)*(fEl->getParmax(d) - fEl->getParmin(d));

      IntVec ig(nsd,0);
      for (int ip = 0; ip < nGP; ip++)
      {
        Real w = vol;
        for (int d = 0; d < nsd; d++)
        {
          const double* xg = GaussQuadrature::getCoord(nG[d]);
          const double* wg = GaussQuadrature::getWeight(nG[d]);
          const Real x0 = fEl->getParmin(d), x1 = fEl->getParmax(d);
          X[d] = Real(0.5)*((x1-x0)*xg[ig[d]] + x1 + x0);
          w *= wg[ig[d]];
        }

        for (const LR::Basisfunction* fi : fEl->support())
        {
          const Real Ni = w*evalBasis(fi,X);
          const size_t ri = 1 + fi->getId();
          for (const LR::Basisfunction* fj : fEl->support())
            mass(ri,1+fj->getId()) += Ni*evalBasis(fj,X);
          for (const LR::Basisfunction* cj : cEl->support())
            B(ri,1+cj->getId()) += Ni*evalBasis(cj,X);
        }

        for (int d = 0; d < nsd; d++)
          if (++ig[d] < nG[d] || d == nsd-1)
            break;
          else
            ig[d] = 0;
      }
    }

    return true;
  }


  //! \brief Returns the LR-spline basis of a patch, or null if it has none.
  //! \param[in] pch The patch to obtain the basis from
  //! \param[in] basis One-based index of the basis
  const LR::LRSpline* getLRBasis (const ASMbase* pch, int basis)
  {
    if (const ASMu2D* p2 = dynamic_cast<const ASMu2D*>(pch); p2)
      return p2->getBasis(basis);
    if (const ASMu3D* p3 = dynamic_cast<const ASMu3D*>(pch); p3)
      return p3->getBasis(basis);

    return nullptr;
  }
#endif


  //! \brief Pairs up the bases of a coarse and a fine patch.
  //! \param[in] cPch The coarse patch
  //! \param[in] fPch The fine patch
  //! \param[in] basis One-based index of the basis
  //! \return The pair, or null if the patches carry no basis of that kind
  std::unique_ptr<BasisPair> basisPair (const ASMbase* cPch,
                                        const ASMbase* fPch, int basis)
  {
#ifdef HAS_LRSPLINE
    if (const LR::LRSpline* cL = getLRBasis(cPch,basis); cL)
      if (const LR::LRSpline* fL = getLRBasis(fPch,basis); fL)
        return std::make_unique<LRPair>(cL,fL);
#endif

    if (const ASMs2D* c2 = dynamic_cast<const ASMs2D*>(cPch); c2)
      if (const ASMs2D* f2 = dynamic_cast<const ASMs2D*>(fPch); f2)
        if (const Go::SplineSurface* cs = c2->getBasis(basis); cs)
          if (const Go::SplineSurface* fs = f2->getBasis(basis); fs)
            return std::make_unique<TensorPair>(tensorBasis(cs,2),
                                                tensorBasis(fs,2));

    if (const ASMs3D* c3 = dynamic_cast<const ASMs3D*>(cPch); c3)
      if (const ASMs3D* f3 = dynamic_cast<const ASMs3D*>(fPch); f3)
        if (const Go::SplineVolume* cv = c3->getBasis(basis); cv)
          if (const Go::SplineVolume* fv = f3->getBasis(basis); fv)
            return std::make_unique<TensorPair>(tensorBasis(cv,3),
                                                tensorBasis(fv,3));

    std::cerr <<" *** MG::prolongation: The two patches carry no spline basis "
              << basis <<" the transfer operators can be built between."
              << std::endl;
    return nullptr;
  }


  //! \brief Returns the local DOF indices selected by a component mask.
  //! \param[in] comps Component mask on the form 1, 2, 12, 123, ..., 0 for all
  //! \param[in] nDofs Number of DOFs at the node
  //!
  //! \details The mask follows the convention of LinSolParams::BlockParams,
  //! where the digits of the number name the components to include.
  IntVec selectDofs (size_t comps, int nDofs)
  {
    IntVec dofs;
    if (comps == 0)
    {
      dofs.resize(nDofs);
      std::iota(dofs.begin(),dofs.end(),1);
      return dofs;
    }

    for (size_t c = comps; c > 0; c /= 10)
      if (int d = c%10; d >= 1 && d <= nDofs)
        dofs.push_back(d);

    std::sort(dofs.begin(),dofs.end());
    return dofs;
  }


  //! \brief Returns the local node number offset of a basis within a patch.
  //! \param[in] pch The patch to consider
  //! \param[in] basis One-based index of the basis
  size_t nodeOffset (const ASMbase& pch, int basis)
  {
    size_t ofs = 0;
    for (int b = 1; b < basis; b++)
      ofs += pch.getNoNodes(b);

    return ofs;
  }


  /*!
    \brief Numbering of the free DOFs an operator is defined on.
    \details The DOFs of a multigrid hierarchy are in general only a subset of
    the equations of the simulator, for instance the pressure equations of a
    Stokes problem. They are numbered consecutively in order of increasing
    global equation number, which is the same order PETSc uses for the index
    set of a matrix block, so that the transfer operators built here line up
    with the blocks of the system matrix.
  */

  //! \brief Maps a local equation to its global number within a block.
  //! \param[in] dd The decomposition holding the numbering
  //! \param[in] iBlk Index of the block in the decomposition, zero for the
  //! whole system
  //! \param[in] eq Local equation number of the whole system
  int globalEq (const DomainDecomposition& dd, size_t iBlk, int eq)
  {
    if (iBlk == 0)
      return dd.getGlobalEq(eq);

    // A block numbers its own equations, so the equation of the system has
    // to be looked up among them before it can be made global.
    const std::map<int,int>& g2l = dd.getG2LEQ(iBlk);
    std::map<int,int>::const_iterator it = g2l.find(eq);
    return it == g2l.end() ? 0 : dd.getGlobalEq(it->second,iBlk);
  }


  class DofNumbering
  {
  public:
    //! \brief The constructor enumerates the free DOFs of an operator.
    //! \param[in] sim The simulator to enumerate the DOFs of
    //! \param[in] op The operator defining the DOF subset
    DofNumbering(const SIMbase& sim, const MG::Operator& op)
    {
      const SAM* sam = sim.getSAM();
      if (!sam) return;

      const DomainDecomposition& dd = sim.getProcessAdm().dd;
      const size_t iBlk = dd.getNoBlocks() > 0 ? op.block+1 : 0;

      const int* madof = sam->getMADOF();
      std::map<int,int> byGlobal;
      for (const ASMbase* pch : sim.getFEModel())
      {
        if (!pch || pch->empty()) continue;

        size_t ofs = nodeOffset(*pch,op.basis);
        size_t nnod = pch->getNoNodes(op.basis);
        for (size_t i = 1; i <= nnod; i++)
        {
          int inod = pch->getNodeID(ofs+i);
          if (inod < 1) continue;

          for (int d : selectDofs(op.comps,madof[inod]-madof[inod-1]))
            if (int eq = sam->getEquation(inod,d); eq > 0)
              if (int geq = globalEq(dd,iBlk,eq); geq > 0)
                byGlobal[geq] = eq;
        }
      }

      // The DOFs are numbered in order of increasing global equation number,
      // which is the order PETSc uses for the index set of a matrix block.
      // The equations a process owns are a contiguous stretch of the global
      // ones, so the DOFs it owns are a contiguous stretch of this numbering,
      // which is what lets the operator be laid out over the processes.
      int idx = 0;
      for (const std::pair<const int,int>& dof : byGlobal)
      {
        index[dof.second] = ++idx;
        if (dof.first >= dd.getMinEq(iBlk) && dof.first <= dd.getMaxEq(iBlk))
          ++nOwned;
      }
    }

    //! \brief Returns the number of free DOFs.
    size_t size() const { return index.size(); }

    //! \brief Returns the number of free DOFs this process owns.
    int owned() const { return nOwned; }

    //! \brief Returns the one-based index of an equation, or zero if not in.
    int operator[](int eq) const
    {
      std::map<int,int>::const_iterator it = index.find(eq);
      return it == index.end() ? 0 : it->second;
    }

  private:
    std::map<int,int> index; //!< Maps local equation number to DOF index
    int nOwned = 0;          //!< DOFs of \a index this process owns
  };
}


//! \brief Maps the basis functions of a patch onto the DOFs of its level.
//! \param[in] pch The patch
//! \param[in] sam Assembly handler of the level
//! \param[in] op The operator the hierarchy is built for
//! \param[in] num DOF numbering of the level
//! \param[out] dofs For each basis function, the DOFs it carries
//!
//! \details A basis function carries one DOF for each component the operator
//! takes, and the operator acts on each of them the same way, so what is
//! built once over the functions is laid out over the DOFs this gives.

static void mapBasisToDofs (const ASMbase& pch, const SAM& sam,
                            const MG::Operator& op, const DofNumbering& num,
                            std::vector<IntVec>& dofs)
{
  const int* madof = sam.getMADOF();
  const size_t ofs = nodeOffset(pch,op.basis);
  dofs.resize(pch.getNoNodes(op.basis));
  for (size_t i = 0; i < dofs.size(); i++)
  {
    dofs[i].clear();
    const int inod = pch.getNodeID(ofs+1+i);
    if (inod < 1) continue;

    for (int d : selectDofs(op.comps,madof[inod]-madof[inod-1]))
      dofs[i].push_back(num[sam.getEquation(inod,d)]);
  }
}


//! \brief Adds the factors of a patch into those of its level.
//! \param[in] Mp Mass matrix of the patch
//! \param[in] Bp The two bases of the patch against each other
//! \param[in] fDof DOFs of the fine basis functions
//! \param[in] cDof DOFs of the coarse basis functions
//! \param mass Mass matrix of the level, added to
//! \param B The two bases of the level against each other, added to

static void scatterFactors (const SparseMatrix& Mp, const SparseMatrix& Bp,
                            const std::vector<IntVec>& fDof,
                            const std::vector<IntVec>& cDof,
                            SparseMatrix& mass, SparseMatrix& B)
{
  for (const auto& [ij,v] : Mp.getValues())
  {
    const IntVec& ri = fDof[ij.first-1];
    const IntVec& ci = fDof[ij.second-1];
    for (size_t d = 0; d < ri.size() && d < ci.size(); d++)
      if (ri[d] > 0 && ci[d] > 0)
        mass(ri[d],ci[d]) += v;
  }

  for (const auto& [ij,v] : Bp.getValues())
  {
    const IntVec& ri = fDof[ij.first-1];
    const IntVec& ci = cDof[ij.second-1];
    for (size_t d = 0; d < ri.size() && d < ci.size(); d++)
      if (ri[d] > 0 && ci[d] > 0)
        B(ri[d],ci[d]) += v;
  }
}


//! \brief Maps coefficients in basis function numbering onto the equations.
//! \param[in] cPch The coarse patch
//! \param[in] fPch The fine patch
//! \param[in] cSam Assembly handler of the coarse level
//! \param[in] fSam Assembly handler of the fine level
//! \param[in] op The operator the hierarchy is built for
//! \param[in] cNum DOF numbering of the coarse level
//! \param[in] fNum DOF numbering of the fine level
//! \param[in] rows Coefficient of each coarse function in each fine function
//! \param P The operator to fill in

static bool mapOntoEquations (const ASMbase& cPch, const ASMbase& fPch,
                              const SAM& cSam, const SAM& fSam,
                              const MG::Operator& op,
                              const DofNumbering& cNum,
                              const DofNumbering& fNum,
                              const std::vector<std::map<int,Real>>& rows,
                              SparseMatrix& P)
{
  const int* fMad = fSam.getMADOF();
  const size_t cOfs = nodeOffset(cPch,op.basis);
  const size_t fOfs = nodeOffset(fPch,op.basis);
  for (size_t i = 0; i < rows.size(); i++)
  {
    int fnod = fPch.getNodeID(fOfs+1+i);
    if (fnod < 1) continue;

    for (int d : selectDofs(op.comps,fMad[fnod]-fMad[fnod-1]))
    {
      int row = fNum[fSam.getEquation(fnod,d)];
      if (row < 1) continue;

      for (const auto& [j,v] : rows[i])
      {
        int cnod = cPch.getNodeID(cOfs+1+j);
        if (cnod < 1) continue;

        if (int col = cNum[cSam.getEquation(cnod,d)]; col > 0)
          P(row,col) = v;
      }
    }
  }

  return true;
}


/*!
  Refinement only inserts knot lines, so the coarse spline space is contained
  in the fine one and every coarse basis function has a unique representation
  in the fine basis. The operator holding those coefficients is what a
  multigrid cycle prolongates with.

  Knot insertion gives those coefficients directly, and how it is carried out
  is up to the pair of bases: a tensor product basis has the change of basis
  of each parameter direction multiply out into the one of the patch, while a
  locally refined one is taken function by function.

  The weights of an insertion lie between zero and one and sum to one, so
  nothing is ever subtracted and no precision is lost however deep the
  refinement goes. That the rows of the operator sum to one is what the two
  bases summing to one leaves behind, and it is checked below.
*/

static bool addPatchTerms (const ASMbase& cPch, const ASMbase& fPch,
                           const SAM& cSam, const SAM& fSam,
                           const MG::Operator& op, const BasisPair& bases,
                           const DofNumbering& cNum, const DofNumbering& fNum,
                           SparseMatrix& P)
{
  // Coefficients of the coarse functions, indexed by fine function. They are
  // kept in basis function numbering here, and mapped onto equations below.
  std::vector<std::map<int,Real>> rows(bases.nFine());

  if (!bases.insertKnots(rows))
  {
    std::cerr <<" *** MG::prolongation: A piece of a coarse basis function is"
              <<" neither a function\n     of the fine basis nor split by any"
              <<" line of the fine mesh. The two\n     meshes are not nested,"
              <<" so there is no transfer operator between them.\n     Levels"
              <<" have to be built by inserting knots into a common geometry,"
              <<"\n     not by removing them from the finest."<< std::endl;
    return false;
  }

  // Every coarse function is a sum of fine ones with weights summing to one,
  // and the coarse functions sum to one, so the fine ones inherit it. Anything
  // else means the pieces were not put back together as they were taken apart.
  for (const std::map<int,Real>& row : rows)
  {
    if (row.empty()) continue; // a function no coarse one reaches

    Real sum = Real(0);
    for (const std::pair<const int,Real>& c : row)
      sum += c.second;

    if (fabs(sum-Real(1)) > Real(1.0e-10))
    {
      std::cerr <<" *** MG::prolongation: The coefficients of a fine basis"
                <<" function sum to "<< sum <<",\n     not to one. The"
                <<" partition of unity the two bases share is not carried"
                <<"\n     over by the operator between them."<< std::endl;
      return false;
    }
  }

  return mapOntoEquations(cPch,fPch,cSam,fSam,op,cNum,fNum,rows,P);
}


std::unique_ptr<MG::Prolongation> MG::prolongation (const SIMbase& coarse,
                                                    const SIMbase& fine,
                                                    const MG::Operator& op,
                                                    MG::Transfer method)
{
  const SAM* cSam = coarse.getSAM();
  const SAM* fSam = fine.getSAM();
  if (!cSam || !fSam)
  {
    std::cerr <<" *** MG::prolongation: The simulators are not preprocessed."
              << std::endl;
    return nullptr;
  }

  const ASM::PatchVec& cModel = coarse.getFEModel();
  const ASM::PatchVec& fModel = fine.getFEModel();
  if (cModel.size() != fModel.size())
  {
    std::cerr <<" *** MG::prolongation: The two levels have a different"
              <<" number of patches, "<< cModel.size() <<" and "
              << fModel.size() <<"."<< std::endl;
    return nullptr;
  }

  // The bases of each patch, paired up coarse with fine
  std::vector<std::unique_ptr<BasisPair>> bases(fModel.size());
  for (size_t i = 0; i < fModel.size(); i++)
    if (cModel[i] && fModel[i] && !cModel[i]->empty() && !fModel[i]->empty())
      if (!(bases[i] = basisPair(cModel[i],fModel[i],op.basis)))
        return nullptr;

  DofNumbering cNum(coarse,op), fNum(fine,op);

  std::unique_ptr<Prolongation> res = std::make_unique<Prolongation>();
  res->rowsOwned = fNum.owned();
  res->colsOwned = cNum.owned();

  // Spaces of a different polynomial order are not nested, whatever their
  // meshes, so there is nothing to insert knots into. The coarse basis is
  // projected onto the fine one instead, and the two factors of that
  // projection are what is kept, the projection itself being dense.
  bool project = method == MG::Transfer::L2_PROJECTION;
  for (size_t i = 0; i < bases.size() && !project; i++)
    if (bases[i] && !bases[i]->sameOrder())
      project = true;

  if (!project)
  {
    res->P = std::make_unique<SparseMatrix>(fNum.size(),cNum.size());
    for (size_t i = 0; i < bases.size(); i++)
      if (bases[i])
        if (!addPatchTerms(*cModel[i],*fModel[i],*cSam,*fSam,op,*bases[i],
                           cNum,fNum,*res->P))
          return nullptr;

    IFEM::cout <<"\tProlongation for \""<< op.name <<"\": "<< res->P->rows()
               <<" x "<< res->P->cols() <<", "<< res->P->size()
               <<" non-zeroes"<< std::endl;
    return res;
  }

  res->B = std::make_unique<SparseMatrix>(fNum.size(),cNum.size());
  res->mass = std::make_unique<SparseMatrix>(fNum.size(),fNum.size());
  res->distributed = true;

  // Each process integrates the elements of the fine level it was given, so
  // that neither the work nor what it produces is repeated on all of them.
  const IntVec& myElms = fine.getProcessAdm().dd.getElms();
  const std::set<int> mine(myElms.begin(),myElms.end());
  for (size_t i = 0; i < bases.size(); i++)
  {
    if (!bases[i])
      continue;

    SparseMatrix Mp(bases[i]->nFine(),bases[i]->nFine());
    SparseMatrix Bp(bases[i]->nFine(),bases[i]->nCoarse());
    if (!bases[i]->integrate(*fModel[i],mine,Mp,Bp))
      return nullptr;

    std::vector<IntVec> fDof, cDof;
    mapBasisToDofs(*fModel[i],*fSam,op,fNum,fDof);
    mapBasisToDofs(*cModel[i],*cSam,op,cNum,cDof);
    scatterFactors(Mp,Bp,fDof,cDof,*res->mass,*res->B);
  }

  IFEM::cout <<"\tProlongation for \""<< op.name <<"\": "<< res->B->rows()
             <<" x "<< res->B->cols() <<", projected through "
             << res->B->size() <<" and "<< res->mass->size()
             <<" non-zeroes"<< std::endl;
  return res;
}
