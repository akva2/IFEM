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
#include "PETScSchurPC.h"
#include "ProcessAdm.h"
#include "LinAlgInit.h"
#include "SAM.h"

#include <algorithm>
#include <numeric>

namespace {

/*!
  \brief This is a C++ version of the F77 subroutine ADDEM2 (SAM library).
  \details It performs exactly the same tasks, except that \a NRHS always is 1
*/
void assemSparse (const Matrix& eM, PETScMatrix& SM, StdVector* SV,
                  const DomainDecomposition& dd,
                  const std::vector<std::array<int,2>>& glb2Blk,
                  const IntVec& meen, const int* meqn,
                  const int* mpmceq, const int* mmceq, const Real* ttcc)
{
  // Get block for equation
  auto getBlk = [&glb2Blk, nBlock = dd.getNoBlocks()](const int ieq, const int jeq)
  {
    if (nBlock > 1)
      return glb2Blk[ieq-1][0] * nBlock + glb2Blk[jeq-1][0];
    else
      return size_t{0};
  };

  // Get equation number in block
  auto getEq = [&glb2Blk, &dd](const int ieq)
  {
    if (dd.getNoBlocks() < 2)
      return dd.getGlobalEq(ieq) - 1;
    else
      return glb2Blk[ieq-1][1] - 1;
  };

  int nedof = meen.size();
  auto A = SM.getBlockMatrices();
  if (A.empty())
    A.push_back(SM.getMatrix());
  for (int j = 1; j <= nedof; ++j) {
    int jeq = meen[j-1];
    if (jeq < 1)
      continue;

    MatSetValue(A[getBlk(jeq, jeq)], getEq(jeq), getEq(jeq), eM(j,j), ADD_VALUES);

    for (int i = 1; i < j; ++i) {
      int ieq = meen[i-1];
      if (ieq < 1)
        continue;

      MatSetValue(A[getBlk(ieq, jeq)], getEq(ieq), getEq(jeq), eM(i,j), ADD_VALUES);
      MatSetValue(A[getBlk(jeq, ieq)], getEq(jeq), getEq(ieq), eM(j,i), ADD_VALUES);
    }
  }

  // Add (appropriately weighted) elements corresponding to constrained
  // (dependent and prescribed) dofs in eM into SM and/or SV
  for (int j = 1; j <= nedof; ++j) {
    int jceq = -meen[j-1];
    if (jceq < 1)
      continue;

    int jp = mpmceq[jceq-1];
    Real c0 = ttcc[jp-1];

    // Add contributions to SV (right-hand-side)
    if (SV)
      for (int i = 1; i <= nedof; ++i) {
        int ieq = meen[i-1];
        int iceq = -ieq;
        if (ieq > 0)
          (*SV)(ieq) -= c0*eM(i,j);
        else if (iceq > 0)
          for (int ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ++ip)
            if (mmceq[ip] > 0) {
              ieq = meqn[mmceq[ip]-1];
              (*SV)(ieq) -= c0*ttcc[ip]*eM(i,j);
            }
      }

    // Add contributions to SM
    for (jp = mpmceq[jceq-1]; jp < mpmceq[jceq]-1; ++jp)
      if (mmceq[jp] > 0) {
        int jeq = meqn[mmceq[jp]-1];
        for (int i = 1; i <= nedof; ++i) {
          int ieq = meen[i-1];
          int iceq = -ieq;
          if (ieq > 0) {
            MatSetValue(A[getBlk(ieq, jeq)], getEq(ieq), getEq(jeq),
                        ttcc[jp]*eM(i,j), ADD_VALUES);
            MatSetValue(A[getBlk(jeq, ieq)], getEq(jeq), getEq(ieq),
                        ttcc[jp]*eM(j,i), ADD_VALUES);
          }
          else if (iceq > 0)
            for (int ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ++ip)
              if (mmceq[ip] > 0) {
                ieq = meqn[mmceq[ip]-1];
                MatSetValue(A[getBlk(ieq, jeq)], getEq(ieq), getEq(jeq),
                            ttcc[ip]*ttcc[jp]*eM(i,j), ADD_VALUES);
              }
        }
      }
  }
}

}


PETScVector::PETScVector(const ProcessAdm& padm) : adm(padm)
{
  VecCreate(*padm.getCommunicator(),&x);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, size_t n)
  : StdVector(n), adm(padm)
{
  if (adm.isParallel())
    n = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;

  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
  LinAlgInit::increfs();
}


PETScVector::PETScVector(const ProcessAdm& padm, const Real* values, size_t n)
  : StdVector(values,n), adm(padm)
{
  if (adm.isParallel())
    n = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;

  VecCreate(*adm.getCommunicator(),&x);
  VecSetSizes(x,n,PETSC_DECIDE);
  VecSetFromOptions(x);
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


bool PETScVector::endAssembly()
{
  // Poor man's assembleDirect
  if (!adm.isParallel() && adm.dd.getMaxDOF() == 0)
    for (size_t i = 0; i < this->size(); ++i)
      VecSetValue(x, i, (*this)[i], ADD_VALUES);
  else
    for (size_t i = 0; i < this->size(); ++i)
      VecSetValue(x, adm.dd.getGlobalEq(i+1)-1, (*this)[i], ADD_VALUES);

  VecAssemblyBegin(x);
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


PETScMatrix::PETScMatrix (const ProcessAdm& padm, const LinSolParams& spar)
  : SparseMatrix(SUPERLU, 1), nsp(nullptr), adm(padm), solParams(spar, adm)
{
  // Create matrix object, by default the matrix type is AIJ
  MatCreate(*adm.getCommunicator(),&pA);

  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  LinAlgInit::increfs();

  if (spar.getNoBlocks() > 1) {
    matvec.resize(spar.getNoBlocks()*spar.getNoBlocks());
    for (Mat& m : matvec) {
      MatCreate(*adm.getCommunicator(), &m);
      MatSetFromOptions(m);
    }
  }

  setParams = true;
  ISsize = 0;
  nLinSolves = 0;
  assembled = false;
}


PETScMatrix::PETScMatrix(const ProcessAdm& padm, const PETScSolParams& spar,
                         const SparseMatrix& A)
    : SparseMatrix(A), nsp(nullptr), adm(padm), solParams(spar)
{
  // Create linear solver object
  KSPCreate(*adm.getCommunicator(),&ksp);

  LinAlgInit::increfs();

  setParams = true;
  ISsize = 0;
  nLinSolves = 0;
  assembled = false;
}


PETScMatrix::~PETScMatrix ()
{
  // Deallocation of linear solver object.
  KSPDestroy(&ksp);

  // Deallocation of matrix object.
  MatDestroy(&pA);
  LinAlgInit::decrefs();
  for (Mat& m : matvec)
    MatDestroy(&m);

  for (IS& v : isvec)
    ISDestroy(&v);

  matvec.clear();
}


SystemMatrix* PETScMatrix::copy() const
{
  PETScMatrix* result = new PETScMatrix(this->adm, this->solParams,
                                        static_cast<const SparseMatrix&>(*this));
  if (this->assembled) {
    MatDuplicate(this->pA, MAT_COPY_VALUES, &result->pA);
    result->assembled = true;
  }
  return result;
}


void PETScMatrix::initAssembly (const SAM& sam, char)
{
  this->resize(sam.neq,sam.neq);
  if (!adm.dd.isPartitioned() && solParams.useSparseMatrix())
    this->preAssemble(sam,false);

  // Get number of local equations in linear system
  PetscInt neq = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;
  // Set correct number of rows and columns for matrix.
  MatSetSizes(pA,neq,neq,PETSC_DETERMINE,PETSC_DETERMINE);

  // Allocate sparsity pattern
  if (matvec.empty()) {
    MatSetFromOptions(pA);

    // Allocate sparsity pattern
    if (adm.dd.isPartitioned())
      this->setupSparsity(adm.dd.getElms(), sam);
    else {
      std::vector<int> elms(sam.nel);
      std::iota(elms.begin(), elms.end(), 0);
      this->setupSparsity(elms, sam);
    }

    MatSetUp(pA);

    switch (solParams.getLinSysType()) {
    case LinAlg::SPD:
      MatSetOption(pA, MAT_SPD, PETSC_TRUE);
      break;
    case LinAlg::SYMMETRIC:
      MatSetOption(pA, MAT_SYMMETRIC, PETSC_TRUE);
      break;
    default:
      break;
    }

#ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    MatSetOption(pA,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
#endif
  } else {
    isvec.resize(adm.dd.getNoBlocks());
    // index sets
    for (size_t i = 0; i < isvec.size(); ++i) {
      IntVec blockEq;
      blockEq.reserve(adm.dd.getMaxEq(i+1) - adm.dd.getMinEq(i+1) + 1);
      for (int leq : adm.dd.getBlockEqs(i)) {
        int eq = adm.dd.getGlobalEq(leq);
        if (eq >= adm.dd.getMinEq() && eq <= adm.dd.getMaxEq())
          blockEq.push_back(eq-1);
      }
      if (adm.dd.isPartitioned())
        std::sort(blockEq.begin(), blockEq.end());

      ISCreateGeneral(*adm.getCommunicator(),blockEq.size(),
                      blockEq.data(),PETSC_COPY_VALUES,&isvec[i]);
    }

    if (adm.dd.isPartitioned())
      this->setupBlockSparsity(adm.dd.getElms(), sam);
    else {
      std::vector<int> elms(sam.nel);
      std::iota(elms.begin(), elms.end(), 0);
      this->setupBlockSparsity(elms, sam);
    }

    MatCreateNest(*adm.getCommunicator(),solParams.getNoBlocks(),isvec.data(),
                  solParams.getNoBlocks(),isvec.data(),matvec.data(),&pA);

 #ifndef SP_DEBUG
    // Do not abort program for allocation error in release mode
    for (Mat& m : matvec)
      MatSetOption(m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
 #endif
  }

  assembled = false;
}


Mat PETScMatrix::preAllocator (const int nrows, const int ncols) const
{
  Mat prealloc;
  MatCreate(*adm.getCommunicator(), &prealloc);
  MatSetType(prealloc, MATPREALLOCATOR);
  MatSetSizes(prealloc, nrows, ncols > 0 ? ncols : nrows,
              PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetUp(prealloc);

  return prealloc;
}


std::vector<Mat>
PETScMatrix::preAllocators () const
{
  const size_t blocks = solParams.getNoBlocks();
  std::vector<Mat> prealloc;
  prealloc.resize(blocks*blocks);
  auto itPre = prealloc.begin();
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++itPre) {
      const int nrows = adm.dd.getMaxEq(i+1) - adm.dd.getMinEq(i+1) + 1;
      const int ncols = adm.dd.getMaxEq(j+1) - adm.dd.getMinEq(j+1) + 1;
      *itPre = preAllocator(nrows, ncols);
    }

  return prealloc;
}


void PETScMatrix::setupSparsity (const std::vector<int>& elms,
                                 const SAM& sam)
{
  PetscInt neq = adm.dd.getMaxEq() - adm.dd.getMinEq() + 1;
  Mat prealloc = preAllocator(neq);

  std::swap(pA, prealloc);
  std::for_each(elms.begin(), elms.end(),
                [this, &sam](const int elm)
                {
                  IntVec meen;
                  sam.getElmEqns(meen, elm+1);
                  this->assemble(Matrix(meen.size(), meen.size()), sam, elm+1);
                });
  this->endAssembly();
  std::swap(pA, prealloc);

  MatPreallocatorPreallocate(prealloc, PETSC_TRUE, pA);

  MatDestroy(&prealloc);
  MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  MatSetOption(pA, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
}


void PETScMatrix::setupBlockSparsity (const std::vector<int>& elms,
                                      const SAM& sam)
{
  auto prealloc = preAllocators();
  this->setupGlb2Blk(sam);

  std::swap(matvec, prealloc);
  std::for_each(elms.begin(), elms.end(),
                [this, &sam](const int elm)
                {
                  IntVec meen;
                  sam.getElmEqns(meen, elm+1);
                  this->assemble(Matrix(meen.size(), meen.size()), sam, elm+1);
                });
  std::swap(matvec, prealloc);

  for (Mat& pmat : prealloc) {
    MatAssemblyBegin(pmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pmat, MAT_FINAL_ASSEMBLY);
  }

  const size_t blocks = solParams.getNoBlocks();

  auto it = matvec.begin();
  auto itPre = prealloc.begin();
  for (size_t i = 0; i < blocks; ++i)
    for (size_t j = 0; j < blocks; ++j, ++it, ++itPre) {
        const int nrows = adm.dd.getMaxEq(i+1) - adm.dd.getMinEq(i+1) + 1;
        const int ncols = adm.dd.getMaxEq(j+1) - adm.dd.getMinEq(j+1) + 1;
        MatSetSizes(*it, nrows, ncols,
                    PETSC_DETERMINE, PETSC_DETERMINE);
        MatPreallocatorPreallocate(*itPre, PETSC_TRUE, *it);
        MatSetUp(*it);

        MatDestroy(&*itPre);
        MatSetOption(*it, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
        MatSetOption(*it, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }
}


void PETScMatrix::setupGlb2Blk (const SAM& sam)
{
  // map from SAM indices to block matrix indices
  size_t blocks = solParams.getNoBlocks();
  glb2Blk.resize(sam.neq, {});
  const DomainDecomposition& dd = adm.dd;

  for (int ieq = 1; ieq <= sam.neq; ++ieq)
    for (size_t b = 0; b < blocks; ++b) {
      if (const auto it = dd.getG2LEQ(b+1).find(ieq);
          it != dd.getG2LEQ(b+1).end())
      {
        glb2Blk[ieq-1][0] = b;
        if (adm.isParallel())
          glb2Blk[ieq-1][1] = adm.dd.isPartitioned()
                                  ? adm.dd.getGlobalEq(ieq, b+1)
                                  : adm.dd.getGlobalEq(it->second, b+1);
        else
          glb2Blk[ieq-1][1] = it->second;
        break;
      }
    }
}


bool PETScMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  if (solParams.useSparseMatrix())
    return this->SparseMatrix::assemble(eM, sam, e);

  IntVec meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;

#pragma omp critical
  assemSparse(eM,*this,nullptr,adm.dd,glb2Blk,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return this->flagNonZeroEqs(meen);
}


bool PETScMatrix::assemble (const Matrix& eM, const SAM& sam,
                            SystemVector& B, int e)
{
  if (solParams.useSparseMatrix())
    return this->SparseMatrix::assemble(eM, sam, B, e);

  StdVector* Bptr = dynamic_cast<StdVector*>(&B);
  if (!Bptr) return false;

  IntVec meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;

#pragma omp critical
  assemSparse(eM,*this,Bptr,adm.dd,glb2Blk,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);

  return this->flagNonZeroEqs(meen);
}


bool PETScMatrix::endAssembly ()
{
  if (solParams.useSparseMatrix()) {
    if (!this->SparseMatrix::endAssembly())
      return false;

    if (IA.empty() && !assembled)
      return this->assembleDirect();
  } else if (!this->getValues().empty())
    return this->assembleDirect();

  for (size_t j = 0; j < cols() && solParams.useSparseMatrix(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      if (matvec.empty())
        MatSetValue(pA,
                    adm.dd.getGlobalEq(JA[i]+1)-1,
                    adm.dd.getGlobalEq(j+1)-1,
                    A[i], ADD_VALUES);
      else {
        const int b = glb2Blk[JA[i]][0] * adm.dd.getNoBlocks() + glb2Blk[j][0];
        MatSetValue(matvec[b],
                    glb2Blk[JA[i]][1] - 1,
                    glb2Blk[j][1] - 1,
                    A[i], ADD_VALUES);
      }

  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);

  assembled = true;

  return true;
}


void PETScMatrix::init ()
{
  this->SparseMatrix::init();

  // Set all matrix elements to zero
  if (matvec.empty())
    MatZeroEntries(pA);
  else for (Mat& m : matvec)
    MatZeroEntries(m);

  assembled = false;
}


bool PETScMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&B);
        PETScVector* Cptr = dynamic_cast<PETScVector*>(&C);

  if ((!Bptr) || (!Cptr))
    return false;

  MatMult(pA,Bptr->getVector(),Cptr->getVector());
  return true;
}


bool PETScMatrix::solve (SystemVector& B, Real*)
{
  PETScVector* Bptr = dynamic_cast<PETScVector*>(&B);
  if (!Bptr)
    return false;

  if (!A.empty() && !assembled)
    return this->solveDirect(*Bptr);

  Vec x;
  VecDuplicate(Bptr->getVector(),&x);
  VecCopy(Bptr->getVector(),x);

  bool result = this->solve(x, Bptr->getVector(),
                            solParams.getStringValue("type") != "preonly");
  VecDestroy(&x);

  return result;
}


bool PETScMatrix::solve (const SystemVector& b, SystemVector& x)
{
  const PETScVector* Bptr = dynamic_cast<const PETScVector*>(&b);
  if (!Bptr)
    return false;

  PETScVector* Xptr = dynamic_cast<PETScVector*>(&x);
  if (!Xptr)
    return false;

  return this->solve(Bptr->getVector(),Xptr->getVector(),false);
}


bool PETScMatrix::solve (const Vec& b, Vec& x, bool knoll)
{
  // Reset linear solver
  if (nLinSolves && solParams.hasValue("reset_pc")) {
    const std::string string_val = solParams.getStringValue("reset_pc");
    int val = solParams.getIntValue("reset_pc");
    if (string_val == "all" ||
        (string_val == "first" && nLinSolves == 1) ||
        (val > 0 && nLinSolves % val == 0)) {
      KSPDestroy(&ksp);
      KSPCreate(*adm.getCommunicator(),&ksp);
      setParams = true;
      factored = false;
      adm.cout << "Resetting preconditioner" << std::endl;
    }
  }

  if (setParams) {
#if PETSC_VERSION_MINOR < 5
    KSPSetOperators(ksp,pA,pA, factored ? SAME_PRECONDITIONER : SAME_NONZERO_PATTERN);
#else
    KSPSetOperators(ksp,pA,pA);
    KSPSetReusePreconditioner(ksp, factored ? PETSC_TRUE : PETSC_FALSE);
#endif
    if (!setParameters(true))
      return false;
    setParams = false;
  }
  if (knoll)
    KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
  else
    KSPSetInitialGuessNonzero(ksp,solParams.getStringValue("type") == "preonly" ?
                                   PETSC_FALSE : PETSC_TRUE);
  KSPSolve(ksp,b,x);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason < 0) {
    adm.cout << "\n Linear solve failed with reason " << KSPConvergedReasons[reason] << std::endl;
    return false;
  }

  if (solParams.getIntValue("verbosity") > 1) {
    PetscInt its;
    KSPGetIterationNumber(ksp,&its);
    adm.cout << "\n Iterations for " << solParams.getStringValue("type")
             << " = " << its << std::endl;
  }
  nLinSolves++;
  factored = true;

  return true;
}


bool PETScMatrix::assembleDirect()
{
  MatSetSizes(pA, PETSC_DETERMINE, PETSC_DETERMINE, this->dim(1), this->dim(2));

  if (this->adm.isParallel())
    MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
  else {
    IntVec iA, jA;
    this->calcCSR(iA,jA);
    MatMPIAIJSetPreallocationCSR(pA, iA.data(), jA.data(), nullptr);
    MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  }

  MatSetUp(pA);
  PetscInt low, high;

  MatGetOwnershipRange(pA, &low, &high);

  for (const auto& e : this->getValues())
    if (static_cast<PetscInt>(e.first.first-1) >= low &&
        static_cast<PetscInt>(e.first.first-1) < high)
      MatSetValue(pA, e.first.first-1, e.first.second-1, e.second, INSERT_VALUES);

  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);
  this->assembled = true;

  return true;
}


bool PETScMatrix::solveDirect(PETScVector& B)
{
  // the sparsity pattern has been grown in-place, we need to init PETsc state.
  // this is currently only used for patch-global L2 systems.
  if (A.empty() && !this->optimiseCols())
    return false;

  // Set correct number of rows and columns for matrix.
  size_t nrow = IA.size()-1;
  if (nrow == 0 || IA.empty())
    return false;

  MatSetSizes(pA, nrow, nrow, PETSC_DECIDE, PETSC_DECIDE);
  MatSetFromOptions(pA);
  PetscInt max = 0;
  for (size_t i = 0; i < nrow; ++i) // symmetric so row/column sizes should be the same
    if (IA[i+1]-IA[i] > max)
      max = IA[i+1]-IA[i];
  MatSeqAIJSetPreallocation(pA, max, nullptr);
  MatSetOption(pA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
  MatSetUp(pA);

  for (size_t j = 0; j < nrow; ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      MatSetValue(pA, JA[i], j, A[i], INSERT_VALUES);

  MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);

  Vec B1, x;
  VecCreate(PETSC_COMM_SELF, &B1);
  VecCreate(PETSC_COMM_SELF, &x);
  VecSetSizes(B1, nrow, PETSC_DECIDE);
  VecSetSizes(x, nrow, PETSC_DECIDE);
  VecSetFromOptions(B1);
  VecSetFromOptions(x);

  size_t nrhs = B.dim() / nrow;
  PetscScalar* bv;
  VecGetArray(B.getVector(), &bv);
  for (size_t i = 0; i < nrhs; ++i) {
    for (size_t j = 0; j < nrow; ++j)
      VecSetValue(B1, j, bv[j*nrow+i], INSERT_VALUES);

    VecAssemblyBegin(B1);
    VecAssemblyEnd(B1);

    if (!this->solve(B1, x, false))
      return false;
    PetscScalar* aa;
    VecGetArray(x, &aa);
    std::copy(aa, aa+nrow, B.getPtr()+i*nrow);
    std::copy(aa, aa+nrow, bv+i*nrow);
    VecRestoreArray(x, &aa);
  }
  VecRestoreArray(B.getVector(), &bv);

  VecDestroy(&x);
  VecDestroy(&B1);

  return true;
}


bool PETScMatrix::solveEig (PETScMatrix& B, RealArray& val,
                            Matrix& vec, int nv, Real shift, int iop)
{
#ifdef HAS_SLEPC
  EPS eps;
  EPSCreate(*adm.getCommunicator(),&eps);

  const auto slepc_mode = std::array{EPS_HEP, EPS_NHEP, EPS_GHEP, EPS_GNHEP};
  EPSSetOperators(eps, pA, iop > 2 ? B.pA : nullptr);
  EPSSetProblemType(eps, slepc_mode[iop-1]);

  EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
  EPSSetDimensions(eps, nv, PETSC_DETERMINE, PETSC_DETERMINE);
  EPSSetFromOptions(eps);

  ST st;
  EPSGetST(eps, &st);
  STSetShift(st, shift);

  KSP oldKsp = ksp;
  STGetKSP(st, &ksp);
  this->setParameters(false);
  ksp = oldKsp;

  if (solParams.getIntValue("verbosity") > 0)
    EPSView(eps, PETSC_VIEWER_STDOUT_WORLD);

  EPSSolve(eps);

  PetscInt nconv;
  EPSGetConverged(eps, &nconv);

  PetscInt m, n;
  MatGetSize(pA, &m, &n);
  if (m != n)
    return false;

  Vec xr, xi;
  MatCreateVecs(pA, nullptr, &xr);
  VecDuplicate(xr, &xi);

  val.resize(nv);
  vec.resize(n, nv);

  Vec gr;
  VecScatter ctx;
  if (adm.isParallel())
    VecScatterCreateToAll(xr, &ctx, &gr);

  for (int i = 0; i < std::min(nv, nconv); ++i) {
    PetscScalar kr, ki;
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    val[i] = kr;
    if (adm.isParallel()) {
      VecScatterBegin(ctx, xr, gr, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(ctx, xr, gr, INSERT_VALUES, SCATTER_FORWARD);
      PetscScalar* grarr;
      VecGetArray(gr, &grarr);
      if (adm.dd.isPartitioned())
        for (const auto& it : adm.dd.getG2LEQ(0))
          vec(it.second,i+1) = grarr[it.first-1];
      VecRestoreArray(gr, &grarr);
    } else {
      PetscScalar* xrarr;
      VecGetArray(xr, &xrarr);
      vec.fillColumn(i+1,xrarr);
      VecRestoreArray(xr, &xrarr);
    }
  }

  VecDestroy(&xi);
  VecDestroy(&xr);

  if (adm.isParallel()) {
    VecDestroy(&gr);
    VecScatterDestroy(&ctx);
  }

  EPSDestroy(&eps);

  return true;
#else
  return false;
#endif
}


Real PETScMatrix::Linfnorm () const
{
  PetscReal norm;
  MatNorm(pA,NORM_INFINITY,&norm);
  return norm;
}


bool PETScMatrix::setParameters (bool setup)
{
  // Set linear solver method
  KSPSetType(ksp,
             !forcedKSPType.empty() ? forcedKSPType.c_str()
                                    : solParams.getStringValue("type").c_str());
  KSPSetTolerances(ksp,solParams.getDoubleValue("rtol"),
                   solParams.getDoubleValue("atol"),
                   solParams.getDoubleValue("dtol"),
                   solParams.getIntValue("maxits"));
  PC pc;
  KSPGetPC(ksp,&pc);

  if (matvec.empty())
    solParams.setupPC(pc, 0, "", IntSet(), setup);
  else if (matvec.size() > 4) {
    std::cerr << "** PETSCMatrix ** Only two blocks supported for now." << std::endl;
    return false;
  }
  else {
    PCSetType(pc,PCFIELDSPLIT);
    PetscInt nsplit;
    KSP  *subksp;
    PC   subpc[2];

    PCFieldSplitSetIS(pc,"u",isvec[0]);
    PCFieldSplitSetIS(pc,"p",isvec[1]);
    PCFieldSplitSetType(pc,PC_COMPOSITE_SCHUR);
    if (solParams.getStringValue("schur") == "lower")
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    else if (solParams.getStringValue("schur") == "full")
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_FULL);
    else if (solParams.getStringValue("schur") == "diag")
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_DIAG);
    else
      PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_UPPER);

    PCFieldSplitSetSchurPre(pc,PC_FIELDSPLIT_SCHUR_PRE_SELFP,nullptr);

    PCSetFromOptions(pc);
    if (setup)
      PCSetUp(pc);
    PCFieldSplitGetSubKSP(pc,&nsplit,&subksp);

    // Preconditioner for blocks
    char pchar='1';
    for (PetscInt m = 0; m < nsplit; m++, pchar++) {
      std::string prefix;
      if (nsplit == 2) {
        if (m == 0)
          prefix = "fieldsplit_u";
        else
          prefix = "fieldsplit_p";
      } else
        prefix = std::string("fieldsplit_b")+pchar;

      KSPSetType(subksp[m],"preonly");
      KSPGetPC(subksp[m],&subpc[m]);
      if (solParams.getBlock(m).getStringValue("pc") == "schur")
        new PETScSchurPC(subpc[m], matvec, solParams.getBlock(m), adm);
      else
        solParams.setupPC(subpc[m], m, prefix, adm.dd.getBlockEqs(m), setup);
    }
  }

  KSPSetFromOptions(ksp);
  if (setup)
    KSPSetUp(ksp);

  if (setup && solParams.getIntValue("verbosity") >= 1)
    KSPView(ksp, PETSC_VIEWER_STDOUT_(*adm.getCommunicator()));

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
