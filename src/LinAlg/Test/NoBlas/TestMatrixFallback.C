//==============================================================================
//!
//! \file TestMatrixFallback.C
//!
//! \date Apr 11 2016
//!
//! \author Eivind Fonn / SINTEF
//!
//! \brief Unit tests for matrix and matrix3d.
//!
//==============================================================================

#undef USE_MKL
#undef USE_ACCELERATE
#undef USE_CBLAS
#undef HAS_BLAS

#include "matrix.h"

#include "../MatrixTests.h"


TEST(TestMatrixFallback, Multiply)
{
  multiplyTest<utl::vector, utl::matrix>();
}


TEST(TestMatrixFallback, Norm)
{
  normTest<utl::matrix>();
}
