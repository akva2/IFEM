//==============================================================================
//!
//! \file SIMPC.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Use a simulator as a preconditioner for iterative solvers.
//!
//==============================================================================
#ifndef SIM_PC_H_
#define SIM_PC_H_

#ifdef HAS_PETSC

#include "AlgEqSystem.h"
#include "ASMs2D.h"
#include "SAM.h"
#include "PETScSolParams.h"
#include <GoTools/geometry/SplineSurface.h>


/*! \brief Template for K-cycle preconditioning.
*/

template<class Sim>
class SIMKCyclePC : public PETScPC
{
public:
  //! \brief Default constructor.
  //! \param fromSim The model for the solution
  //! \param pcSim The model for the preconditioner
  SIMKCyclePC(Sim& fromSim, Sim& pcSim) : 
    S1(pcSim), fSim(fromSim),
    B(SparseMatrix::SUPERLU), BT(SparseMatrix::SUPERLU) {}

  bool setupPC()
  {
    std::array<RealArray,2> gPrm;
    const ASMs2D* fpch = static_cast<const ASMs2D*>(fSim.getPatch(1));
    const ASMs2D* pch = static_cast<const ASMs2D*>(S1.getPatch(1));
    fpch->getGrevilleParameters(gPrm[0], 0);
    fpch->getGrevilleParameters(gPrm[1], 1);
    N.redim(gPrm[0].size()*gPrm[1].size(), S1.getNoDOFs());
    NT.redim(S1.getNoDOFs(), gPrm[0].size()*gPrm[1].size());
    std::vector<Go::BasisDerivsSf> spline;
    pch->getBasis(1)->computeBasisGrid(gPrm[0], gPrm[1], spline);
    int p1 = pch->getBasis(1)->order_u();
    int p2 = pch->getBasis(1)->order_v();
    int n1 = pch->getBasis(1)->numCoefs_u();
    int n2 = pch->getBasis(1)->numCoefs_v();
    for (size_t ip = 0; ip < gPrm[0].size()*gPrm[1].size(); ++ip) {
      IntVec idx;
      ASMs2D::scatterInd(n1,n2,p1,p2,spline[ip].left_idx,idx);
      for (size_t k = 0; k < spline[ip].basisValues.size(); ++k) {
        N(ip+1, idx[k]+1) += spline[ip].basisValues[k];
        NT(idx[k]+1, ip+1) += spline[ip].basisValues[k];
      }
    }
    fpch->getBasis(1)->computeBasisGrid(gPrm[0], gPrm[1], spline);
    p1 = fpch->getBasis(1)->order_u();
    p2 = fpch->getBasis(1)->order_v();
    n1 = fpch->getBasis(1)->numCoefs_u();
    n2 = fpch->getBasis(1)->numCoefs_v();
    B.redim(fSim.getNoDOFs(), fSim.getNoDOFs());
    BT.redim(fSim.getNoDOFs(), fSim.getNoDOFs());

    for (size_t ip = 0; ip < gPrm[0].size()*gPrm[1].size(); ++ip) {
      IntVec idx;
      ASMs2D::scatterInd(n1,n2,p1,p2,spline[ip].left_idx,idx);
      for (size_t k = 0; k < spline[ip].basisValues.size(); ++k) {
        B(ip+1, idx[k]+1) += spline[ip].basisValues[k];
        BT(idx[k]+1, ip+1) += spline[ip].basisValues[k];
      }
    }

    return true;
  }

  //! \brief Evaluate preconditioner.
  bool eval(Vec& x, Vec& y) override
  {
    int len;
    VecGetSize(x, &len);
    StdVector tmp(len);
    std::vector<PetscInt> idx(len);
    std::iota(idx.begin(), idx.end(), 0);
    VecGetValues(x, len, idx.data(), tmp.getPtr());
    Vector sol(fSim.getNoDOFs());
    fSim.getSAM()->expandSolution(tmp, sol);

    StdVector tmp2(sol.size());
    for (size_t i = 1; i <= sol.size(); ++i)
      tmp2(i) = sol(i);

    BT.solve(tmp2);
    StdVector tmp3(NT.rows());
    NT.multiply(tmp2, tmp3);
    Vector solCoarse;
    solCoarse.resize(tmp3.dim());

    // copy to solution vector
    PETScVector* eqVec = dynamic_cast<PETScVector*>(S1.getAlgEqSystem()->getVector(0));
    StdVector* eqsVec = dynamic_cast<StdVector*>(S1.getAlgEqSystem()->getVector(0));
    const int* meqn = S1.getSAM()->getMEQN();
    for (size_t i = 1; i <= S1.getPatch(1)->getNoNodes(); ++i) {
      int eq = meqn[i-1];
      if (eq > 0)
        if (eqVec)
          VecSetValue(eqVec->getVector(), eq-1, tmp3(i), INSERT_VALUES);
        else
          (*eqsVec)(eq) = tmp3(i);
    }
    double cond;
    if (!this->S1.solveSystem(solCoarse, 0, &cond, "displacement", false))
      return false;

    for (size_t i = 1; i <= solCoarse.size(); ++i)
      tmp3(i) = solCoarse(i);

    N.multiply(tmp3, tmp2);
    B.solve(tmp2);

    meqn = fSim.getSAM()->getMEQN();
    for (size_t i = 1; i <= fSim.getPatch(1)->getNoNodes(); ++i) {
      int eq = meqn[i-1];
      if (eq > 0)
        VecSetValue(y, eq-1, tmp2(i), INSERT_VALUES);
    }

    return true;
  }

  Sim& S1;   //!< Preconditioner simulator
  Sim& fSim; //!< Solution model
  SparseMatrix B, BT; //!< Greville mass matrix for fSim
  SparseMatrix N; //!< Basis functions of S1 in greville of fSim
  SparseMatrix NT; //!< Stupid no-transpose multiply
};

#endif

#endif
