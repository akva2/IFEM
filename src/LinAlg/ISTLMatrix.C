// $Id$
//==============================================================================
//!
//! \file ISTLMatrix.C
//!
//! \date Mar 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Representation of the system matrix in ISTL format.
//!
//==============================================================================

#include "ISTLMatrix.h"
#include "SAM.h"
#include "LinAlgInit.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


ISTLVector::ISTLVector(const ProcessAdm& padm) : adm(padm)
{
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, size_t n)
  : StdVector(n)
  , adm(padm)
{
  x.resize(n);
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, const Real* values, size_t n) : adm(padm)
{
  x.resize(n);
  this->restore(values);
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ISTLVector& vec) : StdVector(vec), adm(vec.adm)
{
  x = vec.x;
  LinAlgInit::increfs();
}


ISTLVector::~ISTLVector()
{
  LinAlgInit::decrefs();
}


void ISTLVector::init(Real value)
{
  x = value;
  this->StdVector::init(value);
}


void ISTLVector::redim(size_t n)
{
  x.resize(n);
  this->StdVector::redim(n);
}


bool ISTLVector::endAssembly()
{
  for (size_t i = 0; i < x.size(); ++i)
    x[i] = (*this)[i];

  return true;
}


Real ISTLVector::L1norm() const
{
  return x.one_norm();
}


Real ISTLVector::L2norm() const
{
  return x.two_norm();
}


Real ISTLVector::Linfnorm() const
{
  return x.infinity_norm();
}


ISTLMatrix::ISTLMatrix (const ProcessAdm& padm, const LinSolParams& spar)
  : SparseMatrix(SUPERLU,1), adm(padm), solParams(spar,adm)
{
  LinAlgInit::increfs();
}


ISTLMatrix::ISTLMatrix (const ISTLMatrix& B)
  : SparseMatrix(B), adm(B.adm), solParams(B.solParams.get(),B.adm)
{
  iA = B.iA;

  LinAlgInit::increfs();
}


ISTLMatrix::~ISTLMatrix ()
{
  LinAlgInit::decrefs();
}


void ISTLMatrix::preAssemble (const std::vector<IntVec>& MMNPC, size_t nel)
{
  this->SparseMatrix::preAssemble(MMNPC, nel);

  // Compute the nodal sparsity pattern
  std::vector<std::set<int>> dofc(rows());
  int inod, jnod;
  for (size_t iel = 0; iel < nel; iel++)
    for (size_t j = 0; j < MMNPC[iel].size(); j++)
      if ((jnod = MMNPC[iel][j]) > -1) {
        dofc[jnod].insert(jnod+1);
        for (size_t i = 0; i < j; i++)
          if ((inod = MMNPC[iel][i]) > -1) {
            dofc[inod].insert(jnod+1);
            dofc[jnod].insert(inod+1);
          }
      }

  this->setupSparsity(dofc);
}


void ISTLMatrix::initAssembly (const SAM& sam, bool delayLocking)
{
  SparseMatrix::initAssembly(sam, delayLocking);
  SparseMatrix::preAssemble(sam, delayLocking);

  std::vector<std::set<int>> dofc;
  sam.getDofCouplings(dofc);

  this->setupSparsity(dofc);
}

void ISTLMatrix::setupSparsity (const std::vector<std::set<int>>& dofc)
{
  // Set correct number of rows and columns for matrix.
  size_t sum = 0;
  for (const auto& it : dofc)
    sum += it.size();

  iA.setSize(rows(), cols(), sum);
  iA.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < dofc.size(); ++i)
    iA.setrowsize(i,dofc[i].size());
  iA.endrowsizes();

  for (size_t i = 0; i < dofc.size(); ++i)
    for (const auto& it : dofc[i])
      iA.addindex(i, it-1);

  iA.endindices();

  iA = 0;
}


bool ISTLMatrix::assembleDirect ()
{
  IntVec row_idx, col_idx;
  this->calcCSR(row_idx,col_idx);

  iA.setSize(rows(), cols(), A.size());
  iA.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < rows(); ++i)
    iA.setrowsize(i,row_idx[i+1] - row_idx[i]);
  iA.endrowsizes();

  for (size_t i = 0; i < rows(); ++i)
    for (int j = row_idx[i]; j < row_idx[i+1]; ++j)
      iA.addindex(i, col_idx[j]);

  iA.endindices();

  for (const ValueMap::value_type& val : this->elems())
    iA[val.first.first-1][val.first.second-1] = val.second;

  return true;
}


bool ISTLMatrix::endAssembly()
{
  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      iA[JA[i]][j] = A[i];

  assembled = true;
  return true;
}


void ISTLMatrix::init ()
{
  SparseMatrix::init();

  // Set all matrix elements to zero
  iA = 0;
  assembled = false;
}


bool ISTLMatrix::solve (SystemVector& B, Real*)
{
  ISTLVector* Bptr = dynamic_cast<ISTLVector*>(&B);
  if (!Bptr)
    return false;

  if (!assembled)
    return this->solveDirect(*Bptr);

  if (!pre)
    std::tie(solver, pre, op) = solParams.setupPC(iA);

  if (!solver || !pre)
    return false;

  try {
    Dune::InverseOperatorResult r;
    ISTL::Vec b(Bptr->getVector());
    Bptr->getVector() = 0;
    solver->apply(Bptr->getVector(), b, r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Bptr)(i+1) = Bptr->getVector()[i];

  return true;
}


bool ISTLMatrix::solveDirect (ISTLVector& B)
{
  if (A.empty()) {
    for (const auto& elem : this->elems())
      iA[elem.first.first-1][elem.first.second-1] = elem.second;
  } else
    this->endAssembly();

  if (!pre)
    std::tie(solver, pre, op) = solParams.setupPC(iA);

  if (!solver || !pre)
    return false;

  size_t nrhs = B.dim() / rows();

  try {
    for (size_t c = 0; c < nrhs; ++c) {
      Dune::InverseOperatorResult r;
      ISTL::Vec b(rows()), x(rows());
      for (size_t i = 0; i < rows(); ++i)
        b[i] = B(c*rows()+i+1);
      x = 0;
      solver->apply(x, b, r);
      for (size_t i = 0; i < rows(); ++i)
        B(c * rows() + i + 1) = x[i];
    }
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  return true;
}


bool ISTLMatrix::solve (const SystemVector& b, SystemVector& x)
{
  if (A.empty() && !this->assembleDirect())
    return false;

  if (!pre)
    std::tie(solver, pre, op) = solParams.setupPC(iA);

  const ISTLVector* Bptr = dynamic_cast<const ISTLVector*>(&b);
  if (!Bptr || ! solver || !pre)
    return false;

  ISTLVector* Xptr = dynamic_cast<ISTLVector*>(&x);
  if (!Xptr)
    return false;

  try {
    Dune::InverseOperatorResult r;
    solver->apply(Xptr->getVector(),
                  const_cast<ISTL::Vec&>(Bptr->getVector()), r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Xptr)(i+1) = Xptr->getVector()[i];

  return true;
}


Real ISTLMatrix::Linfnorm () const
{
  return iA.infinity_norm();
}
