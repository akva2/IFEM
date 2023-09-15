// $Id$
//==============================================================================
//!
//! \file BasisFunctionVals.h
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function values container.
//!
//==============================================================================

#ifndef _BASIS_FUNCTION_VALS_H
#define _BASIS_FUNCTION_VALS_H

#include "MatVec.h"


/*!
  \brief Struct holding basis function values and derivatives.
*/

struct BasisFunctionVals
{
  Vector     N;    //!< Basis function values
  Matrix    dNdu;  //!< Basis function derivatives
  Matrix3D d2Ndu2; //!< Second order derivatives of basis functions
  Matrix4D d3Ndu3; //!< Third order derivatives of basis functions
};


//! \brief Convenience type alias
using BasisValuesPtrs = std::vector<const BasisFunctionVals*>;


//! \brief Utility class holding a vector of basis function values.
class BasisValues : public std::vector<BasisFunctionVals>
{
public:
  //! \brief Constructor resizing to a given size.
  BasisValues(size_t size) : std::vector<BasisFunctionVals>(size)
  {
    pointers.reserve(this->size());
    for (const BasisFunctionVals& val : *this)
      pointers.push_back(&val);
  }

  //! \brief Cast to a vector with pointers to our values.
  operator const BasisValuesPtrs& () const { return pointers; }

private:
  BasisValuesPtrs pointers; //!< Vector of pointers to our values
};

#endif
