//==============================================================================
//!
//! \file TestSparseatrix.C
//!
//! \date Apr 9 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for sparse matrices
//!
//==============================================================================

#include "SparseMatrix.h"

#include "gtest/gtest.h"


TEST(TestSparseMatrix, CalcCSR)
{
  SparseMatrix Mat1(2,3);
  Mat1(1,1) = 1.0;
  Mat1(2,1) = 2.0;
  Mat1(2,3) = 3.0;
  IntVec IA1, JA1;
  SparseMatrix::calcCSR(IA1, JA1, 0, 2, Mat1.getValues());

  SparseMatrix Mat1b(3,3);
  Mat1b(1,1) = 1.0;
  Mat1b(2,1) = 2.0;
  Mat1b(2,3) = 3.0;
  IntVec IA2, JA2;
  SparseMatrix::calcCSR(IA2, JA2, 0, 2, Mat1b.getValues());

  EXPECT_EQ(IA1.size(), IA2.size());
  for (size_t i = 0; i < IA1.size(); ++i)
    EXPECT_EQ(IA1[i], IA2[i]);

  EXPECT_EQ(JA1.size(), JA2.size());
  for (size_t j = 0; j < JA1.size(); ++j)
    EXPECT_EQ(JA1[j], JA2[j]);

  SparseMatrix Mat2(1,3);
  Mat2(1,1) = 4.0;
  Mat2(1,2) = 5.0;
  Mat2(1,3) = 6.0;
  SparseMatrix::calcCSR(IA1, JA1, 0, 1, Mat2.getValues());

  SparseMatrix Mat2b(3,3);
  Mat2b(3,1) = 4.0;
  Mat2b(3,2) = 5.0;
  Mat2b(3,3) = 6.0;
  SparseMatrix::calcCSR(IA2, JA2, 2, 3, Mat2b.getValues());

  EXPECT_EQ(IA1.size(), IA2.size());
  for (size_t i = 0; i < IA1.size(); ++i)
    EXPECT_EQ(IA1[i], IA2[i]);

  EXPECT_EQ(JA1.size(), JA2.size());
  for (size_t j = 0; j < JA1.size(); ++j)
    EXPECT_EQ(JA1[j], JA2[j]);
}
