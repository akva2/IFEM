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
#include <numeric>
#include <set>
#include <vector>


namespace // anonymous namespace for local helpers
{
#ifdef HAS_LRSPLINE
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


  /*!
    \brief Projects each coarse basis function onto the fine basis.
    \param[in] cB The coarse basis
    \param[in] fB The fine basis
    \param[out] rows For each fine function, its coefficient against each
    coarse function it takes part in
    \return \e false if the projection could not be formed

    \details Two spline spaces of different polynomial order on the same
    knots are not nested, whatever their meshes: a function of the lower
    order is continuous where one of the higher order is differentiable, so
    it is not among them. There is then no change of basis to be had, and the
    coarse function is instead represented by the function of the fine space
    closest to it, which is its projection.

    That is \f${\bf M}_f{\bf P}={\bf B}\f$, with \f${\bf M}_f\f$ the mass
    matrix of the fine basis and \f$B_{ij}=\int\phi^f_i\phi^c_j\f$, both
    integrated over the elements of the fine mesh, each of which lies within
    one element of the coarse one. The inner product is the one of the
    parameter domain rather than of the physical one, which is a projection
    all the same, and a cheaper one to form.
  */

  bool projectOntoFine (const LR::LRSpline* cB, const LR::LRSpline* fB,
                        std::vector<std::map<int,Real>>& rows)
  {
    const int nsd = fB->nVariate();
    const size_t nF = fB->nBasisFunctions();
    const size_t nC = cB->nBasisFunctions();

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

    Matrix M(nF,nF), B(nF,nC);
    RealArray X(nsd);
    for (int iel = 0; iel < fB->nElements(); iel++)
    {
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
          for (const LR::Basisfunction* fj : fEl->support())
            M(1+fi->getId(),1+fj->getId()) += Ni*evalBasis(fj,X);
          for (const LR::Basisfunction* cj : cEl->support())
            B(1+fi->getId(),1+cj->getId()) += Ni*evalBasis(cj,X);
        }

        for (int d = 0; d < nsd; d++)
          if (++ig[d] < nG[d] || d == nsd-1)
            break;
          else
            ig[d] = 0;
      }
    }

    // The mass matrix is factorized once and every coarse function solved
    // against it, they being right-hand sides of the one system.
    RealArray rhs(nF*nC);
    for (size_t j = 0; j < nC; j++)
      for (size_t i = 0; i < nF; i++)
        rhs[j*nF+i] = B(1+i,1+j);

    std::vector<int> pivot(nF,0);
    if (!utl::solve(M,rhs,&pivot))
      return false;

    const Real dropTol = Real(1.0e-12);
    for (size_t j = 0; j < nC; j++)
      for (size_t i = 0; i < nF; i++)
        if (Real v = rhs[j*nF+i]; fabs(v) > dropTol)
          rows[i][j] = v;

    return true;
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
    \brief Expresses each coarse basis function in the fine basis.
    \param[in] cB The coarse basis
    \param[in] fB The fine basis
    \param[out] rows For each fine function, its coefficient against each
    coarse function it takes part in
    \return \e false if a coarse function was not reached

    \details A coarse function is a B-spline on its own local knot vectors,
    and the fine mesh cuts its support with knot lines the coarse mesh does
    not have. Inserting one of those splits it in two B-splines whose knot
    vectors have the knot, weighted so that the two together are what was
    split. Repeating until every piece is a function of the fine basis
    expresses the coarse function in that basis.

    The weights are the ones knot insertion gives, each between zero and one
    and summing to one over a split. Nothing is subtracted anywhere, so no
    precision is lost however deep the refinement goes; this is what the
    alternative of solving for the coefficients on each element cannot do,
    the systems there conditioning like the Bernstein basis of the order and
    passing that on from element to element.
  */

  bool insertKnots (const LR::LRSpline* cB, const LR::LRSpline* fB,
                    std::vector<std::map<int,Real>>& rows)
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
          return false; // no line reaches it, so leave this to the fallback

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
#endif


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


#ifdef HAS_LRSPLINE
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
  Adaptive refinement only inserts knot lines, so the coarse spline space is
  contained in the fine one and every coarse basis function has a unique
  representation in the fine basis. The operator holding those coefficients is
  what a multigrid cycle prolongates with.

  Knot insertion gives those coefficients directly. A coarse function is a
  B-spline on its own local knot vectors, and each knot line of the fine mesh
  cutting its support splits it into two B-splines which carry that knot,
  weighted so that the two together are what was split. Repeating until every
  piece is a function of the fine basis expresses the coarse function in it.

  Both weights of a split lie between zero and one and sum to one, so nothing
  is ever subtracted and no precision is lost however deep the refinement
  goes. That the rows of the operator sum to one is what the two bases summing
  to one leaves behind, and it is checked below.
*/

static bool addPatchTerms (const ASMbase& cPch, const ASMbase& fPch,
                           const SAM& cSam, const SAM& fSam,
                           const MG::Operator& op,
                           const DofNumbering& cNum, const DofNumbering& fNum,
                           SparseMatrix& P)
{
  const LR::LRSpline* cB = getLRBasis(&cPch,op.basis);
  const LR::LRSpline* fB = getLRBasis(&fPch,op.basis);
  if (!cB || !fB)
  {
    std::cerr <<" *** MG::prolongation: Patch has no LR-spline basis "
              << op.basis <<". Geometric multigrid needs an adaptive"
              <<" discretization."<< std::endl;
    return false;
  }

  // Coefficients of the coarse functions, indexed by fine function. They are
  // kept in basis function numbering here, and mapped onto equations below.
  std::vector<std::map<int,Real>> rows(fB->nBasisFunctions());

  // Spaces of a different polynomial order are not nested, whatever their
  // meshes, so there is nothing to insert knots into and the coarse basis is
  // projected onto the fine one instead.
  bool sameOrder = true;
  for (int d = 0; d < fB->nVariate() && sameOrder; d++)
    sameOrder = cB->order(d) == fB->order(d);

  if (!sameOrder)
  {
    if (!projectOntoFine(cB,fB,rows))
      return false;

    return mapOntoEquations(cPch,fPch,cSam,fSam,op,cNum,fNum,rows,P);
  }

  if (!insertKnots(cB,fB,rows))
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
#endif


std::unique_ptr<SparseMatrix> MG::prolongation (const SIMbase& coarse,
                                                const SIMbase& fine,
                                                const MG::Operator& op,
                                                MG::Transfer method,
                                                int* rowsOwned, int* colsOwned)
{
  if (method == MG::Transfer::L2_PROJECTION)
  {
    std::cerr <<" *** MG::prolongation: The L2-projection transfer operator"
              <<" is not implemented yet.\n     It is only needed for spaces"
              <<" which are not nested, whereas adaptive refinement always"
              <<"\n     gives nested spaces, where the change of basis is"
              <<" both exact and cheaper."<< std::endl;
    return nullptr;
  }

#ifdef HAS_LRSPLINE
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

  DofNumbering cNum(coarse,op), fNum(fine,op);
  if (rowsOwned) *rowsOwned = fNum.owned();
  if (colsOwned) *colsOwned = cNum.owned();

  std::unique_ptr<SparseMatrix> P = std::make_unique<SparseMatrix>(fNum.size(),
                                                                   cNum.size());
  for (size_t i = 0; i < fModel.size(); i++)
    if (cModel[i] && fModel[i] && !cModel[i]->empty() && !fModel[i]->empty())
      if (!addPatchTerms(*cModel[i],*fModel[i],*cSam,*fSam,op,cNum,fNum,*P))
        return nullptr;

  IFEM::cout <<"\tProlongation for \""<< op.name <<"\": "<< P->rows()
             <<" x "<< P->cols() <<", "<< P->size() <<" non-zeroes"
             << std::endl;

  return P;
#else
  std::cerr <<" *** MG::prolongation: Built without LR-spline support."
            << std::endl;
  return nullptr;
#endif
}
