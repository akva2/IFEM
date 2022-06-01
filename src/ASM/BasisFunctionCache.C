// $Id$
//==============================================================================
//!
//! \file BasisFunctionCache.C
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function cache.
//!
//==============================================================================

#include "BasisFunctionCache.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


void BasisFunctionCache::init (int nd)
{
  if (!values.empty() && nderiv == nd)
    return;

  values.clear();
  nderiv = nd;
  this->internalInit();

  if (policy == NO_CACHE) {
#ifdef USE_OPENMP
    size_t size = omp_get_max_threads();
#else
    size_t size = 1;
#endif
    values.resize(size);
    return;
  }
  this->calculateAll();
}


void BasisFunctionCache::finalizeAssembly ()
{
  if (policy == PRE_CACHE)
    values.clear();
}


const BasisFunctionVals& BasisFunctionCache::getVals (size_t el, size_t gp)
{
  if (policy == NO_CACHE) {
#ifdef USE_OPENMP
    size_t idx = omp_get_thread_num();
#else
    size_t idx = 0;
#endif
    values[idx] = this->calculatePt(el, gp);
    return values[idx];
  }

  size_t idx = this->index(el) + gp;
  if (policy == ON_THE_FLY && values[idx].N.empty())
    values[idx] = this->calculatePt(el, gp);

  return values[idx];
}
