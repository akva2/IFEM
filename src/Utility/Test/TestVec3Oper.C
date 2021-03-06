//==============================================================================
//!
//! \file TestVec3Oper.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for global algebraic operations involving Vec3 objects.
//!
//==============================================================================

#include "Vec3.h"
#include "Vec3Oper.h"

#include "gtest/gtest.h"

TEST(TestVec3Oper, GetAndSet)
{
  Vec3 a(1.0,2.0,3.0);
  EXPECT_TRUE(a[0] == 1.0);
  EXPECT_TRUE(a[1] == 2.0);
  EXPECT_TRUE(a[2] == 3.0);
  EXPECT_TRUE(a(1) == 1.0);
  EXPECT_TRUE(a(2) == 2.0);
  EXPECT_TRUE(a(3) == 3.0);
  EXPECT_TRUE(a.x  == 1.0);
  EXPECT_TRUE(a.y  == 2.0);
  EXPECT_TRUE(a.z  == 3.0);
  a[0] = 4.0;
  a(2) = 5.0;
  a.z  = 6.0;
  EXPECT_TRUE(a[0] == 4.0);
  EXPECT_TRUE(a[1] == 5.0);
  EXPECT_TRUE(a[2] == 6.0);
}

TEST(TestVec3Oper, MxV)
{
  utl::matrix<Real> A(3,3);
  A.fill(std::vector<Real>({1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0}).data());
  std::vector<Real> x({1.0,2.0,3.0});
  Vec3 result = A * x;

  EXPECT_TRUE(result.x == 30.0);
  EXPECT_TRUE(result.y == 36.0);
  EXPECT_TRUE(result.z == 42.0);
}

TEST(TestVec3Oper, MultiplyScalar)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 result = a*10.0;
  Vec3 result2 = 10.0*a;

  EXPECT_TRUE(result == result2);
  EXPECT_TRUE(result.x == 10.0);
  EXPECT_TRUE(result.y == 20.0);
  EXPECT_TRUE(result.z == 30.0);
}

TEST(TestVec3Oper, DivideScalar)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 result = a/10.0;

  EXPECT_TRUE(result.x == 0.1);
  EXPECT_TRUE(result.y == 0.2);
  EXPECT_TRUE(result.z == 0.3);
}

TEST(TestVec3Oper, Dot)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Real result = a * b;

  EXPECT_TRUE(result == 32.0);
}

TEST(TestVec3Oper, Addition)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 result = a + b;

  EXPECT_TRUE(result.x == 5.0);
  EXPECT_TRUE(result.y == 7.0);
  EXPECT_TRUE(result.z == 9.0);
}

TEST(TestVec3Oper, Subtraction)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 result = a - b;

  EXPECT_TRUE(result.x == -3.0);
  EXPECT_TRUE(result.y == -3.0);
  EXPECT_TRUE(result.z == -3.0);
}

TEST(TestVec3Oper, Length)
{
  Vec3 a(1.0,1.0,1.0);
  double sqrt3 = sqrt(3.0);

  EXPECT_FLOAT_EQ(a.length(), sqrt3);
}

TEST(TestVec3Oper, Equality)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 c(1.0,2.0,3.0);

  EXPECT_TRUE (a == c);
  EXPECT_FALSE(a == b);
}

TEST(TestVec3Oper, InEquality)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);
  Vec3 c(1.0,2.0,3.0);

  EXPECT_TRUE (a != b);
  EXPECT_FALSE(a != c);
}

TEST(TestVec3Oper, Less)
{
  Vec3 a(1.0,2.0,3.0);
  Vec3 b(4.0,5.0,6.0);

  EXPECT_TRUE(a < b);
  EXPECT_FALSE(b < a);
}

TEST(TestVec3Oper, StreamOut)
{
  Vec3 a(1.0,2.0,3.0);
  std::stringstream str;
  str << a;

  EXPECT_STREQ(str.str().c_str(), "1 2 3");
}

TEST(TestVec3Oper, StreamIn)
{
  std::stringstream str;
  str << "1.0 2.0 3.0";
  Vec3 result;
  str >> result;

  EXPECT_TRUE(result.x == 1.0);
  EXPECT_TRUE(result.y == 2.0);
  EXPECT_TRUE(result.z == 3.0);
}
