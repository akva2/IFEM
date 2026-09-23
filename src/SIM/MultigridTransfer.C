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
#include "GaussQuadrature.h"
#include "IFEM.h"
#include "LogStream.h"
#include "MatVec.h"
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
#endif

#include <algorithm>
#include <cmath>
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


  //! \brief Evaluates an LR-spline basis function in a point.
  //! \param[in] f The basis function to evaluate
  //! \param[in] X The parameter point
  double evalBasis (const LR::Basisfunction* f, const RealArray& X)
  {
    return X.size() == 2 ? f->evaluate(X[0],X[1]) : f->evaluate(X[0],X[1],X[2]);
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

      const int* madof = sam->getMADOF();
      std::set<int> eqs;
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
              eqs.insert(eq);
        }
      }

      int idx = 0;
      for (int eq : eqs)
        index[eq] = ++idx;
    }

    //! \brief Returns the number of free DOFs.
    size_t size() const { return index.size(); }

    //! \brief Returns the one-based index of an equation, or zero if not in.
    int operator[](int eq) const
    {
      std::map<int,int>::const_iterator it = index.find(eq);
      return it == index.end() ? 0 : it->second;
    }

  private:
    std::map<int,int> index; //!< Maps global equation number to DOF index
  };
}


#ifdef HAS_LRSPLINE
/*!
  Adaptive refinement only inserts knot lines, so the coarse spline space is
  contained in the fine one and every coarse basis function has a unique
  representation in the fine basis. The operator holding those coefficients is
  what a multigrid cycle prolongates with.

  A fine element lies entirely within one coarse element, and both bases
  restrict to polynomials of the same degree there. Evaluating them in a set of
  points which is unisolvent for that polynomial space therefore gives
  \f${\bf A}_f{\bf P} = {\bf A}_c\f$ for the rows of \b P belonging to the
  functions on the element. When the element carries exactly \a p+1 functions
  per parameter direction, \f${\bf A}_f\f$ is square and invertible, and the
  rows come straight out.

  LR-spline meshes may however be \e overloaded: an element can carry more
  functions than that, in which case their restrictions to it are linearly
  dependent and \f${\bf A}_f\f$ is singular. Structured mesh refinement does
  not prevent this. It only guarantees that the mesh as a whole is linearly
  independent, which is the weaker property tested by
  LRSplineSurface::isLinearIndepByOverloading. Overloaded elements do occur in
  practice, roughly a quarter of them after a few rounds of local refinement.

  The way out uses the locality of knot insertion: a fine function contributing
  to a coarse one has its support contained in the support of that coarse
  function. A single element of \f$\mbox{supp}(N^f_i)\f$ which is not
  overloaded therefore determines the whole of row \a i, and since every row
  found this way is a row of the unique global operator, rows found on
  different elements never disagree.

  So the rows are peeled off: the elements which are not overloaded are
  resolved first, and the overloaded ones are then revisited with the rows
  already known moved to the right-hand side, which removes the functions
  causing the dependency. This is the same argument the overloading test itself
  makes, so it terminates on exactly those meshes the test accepts.
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

  const int nsd = fB->nVariate();

  // Coefficients below this are roundoff, not structure. Dropping them keeps
  // the sparsity pattern of the operator, and of the Galerkin products formed
  // from it, down to the entries which are really there. The coefficients of
  // a prolongation between spline spaces are O(1) and its rows sum to one, so
  // an absolute tolerance is well scaled here.
  const Real dropTol = Real(1.0e-12);

  // One Gauss point per polynomial order in each direction. This is unisolvent
  // for the polynomials living on the element, and keeps the evaluation points
  // well inside it.
  int nGP = 1;
  IntVec nG(nsd);
  for (int d = 0; d < nsd; d++)
  {
    if (!GaussQuadrature::getCoord(nG[d] = fB->order(d)))
    {
      std::cerr <<" *** MG::prolongation: No Gauss rule with "<< nG[d]
                <<" points, needed for a basis of order "<< nG[d] <<"."
                << std::endl;
      return false;
    }
    nGP *= nG[d];
  }

  // Locate the coarse element containing each fine element. Since the fine
  // mesh refines the coarse one, the element midpoint identifies it.
  RealArray X(nsd);
  IntVec parent(fB->nElements(),-1);
  for (int iel = 0; iel < fB->nElements(); iel++)
  {
    const LR::Element* fEl = fB->getElement(iel);
    for (int d = 0; d < nsd; d++)
      X[d] = 0.5*(fEl->getParmin(d) + fEl->getParmax(d));

    if ((parent[iel] = cB->getElementContaining(X)) < 0)
    {
      std::cerr <<" *** MG::prolongation: No coarse element contains the"
                <<" midpoint of fine element "<< 1+iel <<"."<< std::endl;
      return false;
    }
  }

  // Coefficients of the coarse functions, indexed by fine function. They are
  // kept in basis function numbering here, and mapped onto equations below.
  std::vector<std::map<int,Real>> rows(fB->nBasisFunctions());
  std::vector<bool> known(fB->nBasisFunctions(),false);

  size_t nLeft = rows.size();
  Matrix Af, Ac, Au, AtA, AtB, B, Psub;
  IntVec cIdx, fIdx, unknown;
  while (nLeft > 0)
  {
    const size_t nLeftBefore = nLeft;
    for (int iel = 0; iel < fB->nElements() && nLeft > 0; iel++)
    {
      const LR::Element* fEl = fB->getElement(iel);
      const LR::Element* cEl = cB->getElement(parent[iel]);

      fIdx.clear();
      for (const LR::Basisfunction* f : fEl->support())
        fIdx.push_back(f->getId());
      cIdx.clear();
      for (const LR::Basisfunction* f : cEl->support())
        cIdx.push_back(f->getId());

      unknown.clear();
      for (size_t i = 0; i < fIdx.size(); i++)
        if (!known[fIdx[i]])
          unknown.push_back(i);

      // Nothing left to gain from this element, or more unknown rows than the
      // values sampled on it can determine. Such an element is overloaded and
      // is revisited once some of its functions are known from elsewhere.
      if (unknown.empty() || unknown.size() > static_cast<size_t>(nGP))
        continue;

      const size_t nLoc = fIdx.size();
      const size_t nCo = cIdx.size();
      Af.resize(nGP,nLoc);
      Ac.resize(nGP,nCo);

      IntVec ig(nsd,0);
      for (int ip = 1; ip <= nGP; ip++)
      {
        for (int d = 0; d < nsd; d++)
        {
          const double* xg = GaussQuadrature::getCoord(nG[d]);
          double x0 = fEl->getParmin(d), x1 = fEl->getParmax(d);
          X[d] = 0.5*((x1-x0)*xg[ig[d]] + x1 + x0);
        }

        int i = 0;
        for (const LR::Basisfunction* f : fEl->support())
          Af(ip,++i) = evalBasis(f,X);
        i = 0;
        for (const LR::Basisfunction* f : cEl->support())
          Ac(ip,++i) = evalBasis(f,X);

        // Advance the tensorial Gauss point counter
        for (int d = 0; d < nsd; d++)
          if (++ig[d] < nG[d] || d == nsd-1)
            break;
          else
            ig[d] = 0;
      }

      // Move the contributions of the rows already known over to the
      // right-hand side, leaving a system in the unknown rows only
      B = Ac;
      for (size_t i = 0; i < nLoc; i++)
        if (known[fIdx[i]])
          for (const auto& [j,v] : rows[fIdx[i]])
          {
            int jc = utl::findIndex(cIdx,j);
            if (jc < 0) continue;

            for (int ip = 1; ip <= nGP; ip++)
              B(ip,1+jc) -= Af(ip,1+i)*v;
          }

      Au.resize(nGP,unknown.size());
      for (size_t i = 0; i < unknown.size(); i++)
        for (int ip = 1; ip <= nGP; ip++)
          Au(ip,1+i) = Af(ip,1+unknown[i]);

      if (unknown.size() == static_cast<size_t>(nGP))
      {
        // Square system, solved directly rather than through the normal
        // equations, which would square the conditioning
        if (!utl::invert(Au))
          continue; // singular, so the element is overloaded after all
        Psub.multiply(Au,B);
      }
      else
      {
        // Overdetermined, but consistent since the coarse space is contained
        // in the fine one, so the normal equations give the exact solution
        AtA.multiply(Au,Au,true,false);
        if (!utl::invert(AtA))
          continue;
        AtB.multiply(Au,B,true,false);
        Psub.multiply(AtA,AtB);
      }

      for (size_t i = 0; i < unknown.size(); i++)
      {
        std::map<int,Real>& row = rows[fIdx[unknown[i]]];
        for (size_t j = 0; j < nCo; j++)
          if (Real v = Psub(1+i,1+j); fabs(v) > dropTol)
            row[cIdx[j]] = v;

        known[fIdx[unknown[i]]] = true;
        --nLeft;
      }
    }

    if (nLeft == nLeftBefore)
    {
      std::cerr <<" *** MG::prolongation: Could not resolve "<< nLeft <<" of "
                << rows.size() <<" basis functions.\n     The mesh is"
                <<" overloaded beyond what peeling can resolve, which means"
                <<" it is not\n     linearly independent either."<< std::endl;
      return false;
    }
  }

  // Map the coefficients onto the equations of the two levels
  const int* cMad = cSam.getMADOF();
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
#endif


std::unique_ptr<SparseMatrix> MG::prolongation (const SIMbase& coarse,
                                                const SIMbase& fine,
                                                const MG::Operator& op,
                                                MG::Transfer method)
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
