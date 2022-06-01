// $Id$
//==============================================================================
//!
//! \file BasisFunctionCache.h
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function cache.
//!
//==============================================================================

#ifndef _BASIS_FUNCTION_CACHE_H
#define _BASIS_FUNCTION_CACHE_H

#include "matrixnd.h"
#include "MatVec.h"


struct BasisFunctionVals {
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;
  Matrix4D d3Ndu3;
};


class BasisFunctionCache
{
public:
  enum Policy {
    NO_CACHE,   //!< Cache is disabled - calculate on the fly
    PRE_CACHE,  //!< Cache basis function values up front, clear on assembly end
    ON_THE_FLY, //!< Cache basis functions on the fly
    FULL_CACHE  //!< Cache basis function values up front
  };

  void finalizeAssembly();

  void init(int nd);

  const BasisFunctionVals& getVals(size_t el, size_t gp);

  virtual double getParam(int dir, size_t el, size_t gp) = 0;

  virtual size_t index(size_t el) = 0;

protected:
  virtual bool internalInit() = 0;
  virtual BasisFunctionVals calculatePt(size_t el, size_t gp) = 0;
  virtual void calculateAll() = 0;

  Policy policy = PRE_CACHE;
  std::vector<BasisFunctionVals> values;
  int nderiv = 0;
};

#endif
