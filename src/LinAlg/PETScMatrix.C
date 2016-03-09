// $Id$
//==============================================================================
//!
//! \file PETScMatrix.C
//!
//! \date Jan 15 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the system matrix in PETSc format.
//!
//==============================================================================

#include "PETScMatrix.h"
#include "LinSolParams.h"
#include "LinAlgInit.h"
#include "SAMpatchPETSc.h"
#include "ProcessAdm.h"
#include "SIMenums.h"
#include "ASMstruct.h"
#include "DomainDecomposition.h"
#include "Profiler.h"


PETScVector::PETScVector(const ProcessAdm& padm) : adm(padm)
{
  VecCreate(*padm.getCommunicator(),&x);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, size_t n) :
  StdVector(n), adm(padm)
{
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,adm.dd.getMaxEq()-adm.dd.getMinEq()+1,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, const Real* values, size_t n) :
  StdVector(values, n), adm(padm)
{
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,adm.dd.getMaxEq()-adm.dd.getMinEq() + 1,PETSC_DECIDE);
  VecSetFromOptions(x);
  this->restore(values);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const PETScVector& vec) :
  StdVector(vec), adm(vec.adm)
{
  VecDuplicate(vec.x,&x);
  VecCopy(vec.x,x);
  LinAlgInit::increfs();
}


PETScVector::~PETScVector()
{
  VecDestroy(&x);
  LinAlgInit::decrefs();
}


void PETScVector::init(Real value)
{
  StdVector::init(value);
  VecSet(x,value);
}


void PETScVector::redim(size_t n)
{
  VecDestroy(&x);
  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,adm.dd.getMaxEq()-adm.dd.getMinEq() + 1,PETSC_DECIDE);
  VecSetFromOptions(x);
  StdVector::redim(n);
}


bool PETScVector::beginAssembly()
{
  for (size_t i = 0; i < size(); ++i)
    VecSetValue(x , adm.dd.getGlobalEq(i+1)-1, (*this)[i], ADD_VALUES);

  VecAssemblyBegin(x);
  return true;
}


bool PETScVector::endAssembly()
{
  VecAssemblyEnd(x);
  return true;
}


Real PETScVector::L1norm() const
{
  PetscReal val;

  VecNorm(x,NORM_1,&val);
  return val;
}


Real PETScVector::L2norm() const
{
  PetscReal val;

  VecNorm(x,NORM_2,&val);
  return val;
}


Real PETScVector::Linfnorm() const
{
  PetscReal val;

  VecNorm(x,NORM_INFINITY,&val);
  return val;
}


PETScMatrix::PETScMatrix (const ProcessAdm& padm, const LinSolParams& spar,
                          LinAlg::LinearSystemType ltype) :
 SparseMatrix(SUPERLU, 1), nsp(nullptr), adm(padm), solParams(spar, adm),
 linsysType(ltype)
{
  // Create matrix object, by default the matrix type is AIJ
  MatCreate(*adm.getCommunicator(),&A);

  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  LinAlgInit::increfs();

  if (spar.getNoBlocks() > 1) {
    isvec.resize(spar.getNoBlocks());
    matvec.resize(spar.getNoBlocks()*spar.getNoBlocks());
    for (auto& it : matvec) {
      MatCreate(*adm.getCommunicator(), &it);
      MatSetFromOptions(it);
    }
  }

  setParams = true;
  ISsize = 0;
  nLinSolves = 0;
}


PETScMatrix::~PETScMatrix ()
{
  // Deallocation of linear solver object.
  KSPDestroy(&ksp);

  // Deallocation of matrix object.
  MatDestroy(&A);
  LinAlgInit::decrefs();
  for (auto& it : matvec)
    MatDestroy(&it);
  matvec.clear();
  for (auto& it : isvec)
    ISDestroy(&it);
}


void PETScMatrix::initAssembly (const SAM& sam, bool b)
{
  SparseMatrix::initAssembly(sam, b);
  SparseMatrix::preAssemble(sam, b);

  const SAMpatchPETSc* samp = dynamic_cast<const SAMpatchPETSc*>(&sam);
  if (!samp)
    return;

  // Get number of equations in linear system
  const PetscInt neq  = adm.dd.getMaxEq()- adm.dd.getMinEq() + 1;

  size_t nx = solParams.getBlock(0).subdomains[0];
  size_t ny = solParams.getBlock(0).subdomains[1];
  size_t nz = solParams.getBlock(0).subdomains[2];
  int overlap = solParams.getBlock(0).overlap;

  if (nx+ny+nz > 0) {
    locSubdDofs.resize(nx*ny*samp->getPatches().size());
    size_t d = 0;
    for (const auto& it : samp->getPatches()) {
      const ASMstruct* pch = dynamic_cast<const ASMstruct*>(it);
      if (!pch)
        break;
      int n1, n2, n3;
      pch->getNoStructElms(n1,n2,n3);
      const_cast<DomainDecomposition&>(adm.dd).calcAppropriateGroups(n1, n2, n3, nx, ny, nz, overlap);
      for (size_t g = 0; g < adm.dd.getNoSubdomains(); ++g, ++d) {
        std::set<int> eqnSet;
        for (const auto& iEl : adm.dd[g]) {
          IntVec eqns;
          sam.getElmEqns(eqns, iEl+1);
          for (auto& it : eqns)
            it--;
          eqnSet.insert(eqns.begin(), eqns.end());
        }
        std::copy_if(eqnSet.begin(), eqnSet.end(),
                     std::back_inserter(locSubdDofs[d]), [](const int& a) { return a > -1;});
      }
    }
    subdDofs = locSubdDofs;
  }

  // Set correct number of rows and columns for matrix.
  MatSetSizes(A,neq,neq,PETSC_DECIDE,PETSC_DECIDE);
  MPI_Barrier(*adm.getCommunicator());

  MatSetFromOptions(A);

  // Allocate sparsity pattern
  std::vector<std::set<int>> dofc;
  sam.getDofCouplings(dofc);

  if (matvec.empty()) {
    // Set correct number of rows and columns for matrix.
    MatSetSizes(A,neq,neq,PETSC_DECIDE,PETSC_DECIDE);
    MPI_Barrier(*adm.getCommunicator());

    MatSetFromOptions(A);

    // Allocate sparsity pattern
    if (adm.isParallel()) {
      int ifirst = adm.dd.getMinEq();
      int ilast  = adm.dd.getMaxEq();
      int maxfill = std::min(100, ilast-ifirst+1);
      // !!! SERIOUS OVER/UNDERALLOCATION !!!
      PetscIntVec d_nnz(ilast-ifirst+1, maxfill), o_nnz(ilast-ifirst+1, maxfill);

  //    for (size_t i = 0; dofc.size(); ++i) {
  //      for (const auto& it : dofc[i]) {
  //        if (
  //        if (adm.dd.getGlobalEq(it) >= ifirst && adm.dd.getGlobalEq(it) < ilast)
  //          ++d_nnz[i];
  //        else
  //          ++o_nnz[i];
  //      }
  //    }

  //    Vec nnz;
  //    VecCreate(*adm.getCommunicator(),&nnz);
  //    VecSetSizes(nnz, ilast-ifirst);

      MatMPIAIJSetPreallocation(A,PETSC_DEFAULT,d_nnz.data(),
                                  PETSC_DEFAULT,o_nnz.data());
    } else {
      PetscIntVec Nnz;
      for (const auto& it : dofc)
        Nnz.push_back(it.size());

      MatSeqAIJSetPreallocation(A,PETSC_DEFAULT,Nnz.data());
      MatSeqAIJSetPreallocation(A,PETSC_DEFAULT,Nnz.data());

      PetscIntVec col;
      for (const auto& it2 : dofc)
        for (const auto& it : it2)
          col.push_back(it-1);

      MatSeqAIJSetColumnIndices(A,&col[0]);
      MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
      MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }

    MatSetUp(A);

    if (linsysType == LinAlg::SPD)
      MatSetOption(A, MAT_SPD, PETSC_TRUE);
    if (linsysType == LinAlg::SYMMETRIC)
      MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);

#ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
#endif
  } else {
    // get equations for block
    std::vector<std::vector<int>> blockEqs(solParams.getNoBlocks());
    for (size_t i = 0; i < solParams.getNoBlocks(); ++i) {
      char dofType = solParams.getBlock(i).basis == 1 ? 'D' :
                            'P'+solParams.getBlock(i).basis-2;
      // grab DOFs of given type(s)
      if (solParams.getBlock(i).comps != 0) {
        char temp[32];
        sprintf(temp,"%lu",solParams.getBlock(i).comps);
        for (char* t = temp; *t != 0; ++t) {
          std::vector<int> tmp = adm.dd.getSAM()->getEquations(dofType, (int)(*t - '0'));
          blockEqs[i].insert(blockEqs[i].end(), tmp.begin(), tmp.end());
        }
      } else
        blockEqs[i] = adm.dd.getSAM()->getEquations(dofType);
    }

    // nnz per block
    std::vector<std::vector<PetscInt>> nnz(solParams.getNoBlocks());
    for (size_t i = 0; i < solParams.getNoBlocks(); ++i)
      for (const auto& it : blockEqs[i])
        nnz[i].push_back(dofc[it].size());

    auto it = matvec.begin();
    for (size_t i = 0; i < solParams.getNoBlocks(); ++i)
      for (size_t j = 0; j < solParams.getNoBlocks(); ++j, ++it) {
        MatSetSizes(*it, blockEqs[i].size(), blockEqs[j].size(),
                    PETSC_DETERMINE, PETSC_DETERMINE);
        MatSetFromOptions(*it);
        MatSeqAIJSetPreallocation(*it, PETSC_DEFAULT, nnz[i].data());
      }

    // index sets
    for (size_t i = 0; i < isvec.size(); ++i)
      ISCreateGeneral(*adm.getCommunicator(),blockEqs[i].size(),
                      blockEqs[i].data(),PETSC_COPY_VALUES,&isvec[i]);

    MatCreateNest(*adm.getCommunicator(),solParams.getNoBlocks(),isvec.data(),
                  solParams.getNoBlocks(),isvec.data(),matvec.data(),&A);

 #ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    for (auto& it : matvec)
      MatSetOption(it,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
 #endif
  }
}


bool PETScMatrix::beginAssembly()
{
  this->optimiseSLU();
  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      MatSetValue(A, adm.dd.getGlobalEq(JA[i]+1)-1,
                  adm.dd.getGlobalEq(j+1)-1, SparseMatrix::A[i], ADD_VALUES);

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  return true;
}


bool PETScMatrix::endAssembly()
{
  // Finalizes parallel assembly process
  PROFILE("extra assembly");
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  return true;
}


void PETScMatrix::init ()
{
  SparseMatrix::init();

  // Set all matrix elements to zero
  MatZeroEntries(A);
}


bool PETScMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&B);
        PETScVector* Cptr = dynamic_cast<PETScVector*>(&C);

  if ((!Bptr) || (!Cptr))
    return false;

  MatMult(A,Bptr->getVector(),Cptr->getVector());
  return true;
}


bool PETScMatrix::solve (SystemVector& B, bool newLHS, Real*)
{
  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  StdVector* Bsptr = dynamic_cast<StdVector*>(&B);
  if (!Bsptr)
    return false;

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);

  bool result = this->solve(x,Bptr->getVector(),newLHS,true);
  if (result) {
    if (!glob2LocEq) {
      std::vector<int> mlgeq(adm.dd.getMLGEQ());
      for (auto& it : mlgeq)
        --it;

      ISCreateGeneral(*adm.getCommunicator(),adm.dd.getMLGEQ().size(),
                      mlgeq.data(), PETSC_COPY_VALUES, &glob2LocEq);
    }

    if (adm.isParallel()) {
#ifdef PARALLEL_PETSC
      Vec solution;
      VecCreateSeqWithArray(PETSC_COMM_SELF, 1, adm.dd.getMLGEQ().size(), Bsptr->getPtr(), &solution);
      VecScatter ctx;
      VecScatterCreate(Bptr->getVector(), glob2LocEq, solution, nullptr, &ctx);
      VecScatterBegin(ctx, Bptr->getVector(), solution, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(ctx, Bptr->getVector(),solution,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterDestroy(&ctx);
      VecDestroy(&solution);
#endif
    } else {
      PetscScalar* data;
      VecGetArray(Bptr->getVector(), &data);
      std::copy(data, data + Bptr->dim(), Bptr->getPtr());
      VecRestoreArray(Bptr->getVector(), &data);
    }
  }
  VecDestroy(&x);

  return result;
}


bool PETScMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&b);
  if (!Bptr)
    return false;

  PETScVector* Xptr = dynamic_cast<PETScVector*>(&x);
  if (!Xptr)
    return false;

  bool result = this->solve(Bptr->getVector(),Xptr->getVector(),newLHS,false);
  if (result) {
    PetscScalar* data;
    VecGetArray(Xptr->getVector(), &data);
    std::copy(data, data + Xptr->dim(), Xptr->getPtr());
    VecRestoreArray(Xptr->getVector(), &data);
  }

  return result;
}


bool PETScMatrix::solve (const Vec& b, Vec& x, bool newLHS, bool knoll)
{
  // Reset linear solver
  if (nLinSolves && solParams.getResetSolver())
    if (nLinSolves%solParams.getResetSolver() == 0) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
    }

  if (setParams) {
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,A,A, newLHS ? SAME_NONZERO_PATTERN : SAME_PRECONDITIONER);
#else
    KSPSetOperators(ksp,A,A);
    KSPSetReusePreconditioner(ksp, newLHS ? PETSC_FALSE : PETSC_TRUE);
#endif
    if (!setParameters())
      return false;
    setParams = false;
  }
  if (knoll)
    KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  else
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSolve(ksp,b,x);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason < 0) {
    PetscPrintf(PETSC_COMM_WORLD, "\n Linear solve failed with reason %s",KSPConvergedReasons[reason]);
    return false;
  }

  if (solParams.getMessageLevel() > 1) {
    PetscInt its;
    KSPGetIterationNumber(ksp,&its);
    PetscPrintf(PETSC_COMM_WORLD,"\n Iterations for %s = %D\n",solParams.getMethod().c_str(),its);
  }
  nLinSolves++;

  return true;
}


bool PETScMatrix::solveEig (PETScMatrix& B, RealArray& val,
			    Matrix& vec, int nv, Real shift, int iop)
{
#ifdef HAS_SLEPC
  ST          st;
  PetscInt    m, n, nconv;
  PetscScalar kr, ki;
  PetscScalar *xrarr;
  Vec         xr, xi;

  EPS eps;
  EPSCreate(*adm.getCommunicator(),&eps);

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B.A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B.A,MAT_FINAL_ASSEMBLY);

  EPSSetOperators(eps,A,B.A);
  EPSSetProblemType(eps,EPS_GHEP);
  EPSSetType(eps,EPSKRYLOVSCHUR);
  EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);
  EPSGetST(eps,&st);
  STSetShift(st,shift);
  EPSSetDimensions(eps,nv,4*nv,PETSC_NULL);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetConverged(eps,&nconv);

  MatGetSize(A,&m,&n);
  if (m != n) return false;

  VecCreate(*adm.getCommunicator(),&xr);
  VecSetSizes(xr,n,PETSC_DETERMINE);
  VecSetFromOptions(xr);
  VecDuplicate(xr,&xi);

  val.resize(nv);
  vec.resize(n,nv);
  for (int i = 0; i < nv; i++) {
    EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
    VecGetArray(xr,&xrarr);
    val[i] = kr;
    vec.fillColumn(i+1,xrarr);
    VecRestoreArray(xr,&xrarr);
  }

  VecDestroy(&xi);
  VecDestroy(&xr);

  EPSDestroy(&eps);

  return true;
#endif
  return false;
}


Real PETScMatrix::Linfnorm () const
{
  PetscReal norm;
  MatNorm(A,NORM_INFINITY,&norm);
  return norm;
}


bool PETScMatrix::setParameters(PETScMatrix* P, PETScVector* Pb)
{
  solParams.setParams(ksp,locSubdDofs,subdDofs,coords,dirIndexSet);
  return true;
}


PETScVector operator*(const SystemMatrix& A, const PETScVector& b)
{
  PETScVector results(b.getAdm());
  A.multiply(b, results);
  return results;
}


PETScVector operator/(SystemMatrix& A, const PETScVector& b)
{
  PETScVector results(b.getAdm());
  A.solve(b, results);
  return results;
}
