//==============================================================================
//!
//! \file MatrixTests.h
//!
//! \date Apr 11 2016
//!
//! \author Eivind Fonn / SINTEF
//!
//! \brief Templates implementing unit tests for matrix.
//!
//==============================================================================

#include <numeric>

#include "gtest/gtest.h"


namespace {

template<template<class T> class Vector, template<class T> class Matrix>
void multiplyTest()
{
  Vector<double> u(14), v(9), x, y;
  Matrix<double> A(3,5);

  std::iota(u.begin(),u.end(),1.0);
  std::iota(v.begin(),v.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);

  ASSERT_TRUE(A.multiply(u,x,1.0,0.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,0.0,true,4,2));

  EXPECT_FLOAT_EQ(x(3),370.0);
  EXPECT_FLOAT_EQ(x(7),410.0);
  EXPECT_FLOAT_EQ(x(11),450.0);
  EXPECT_FLOAT_EQ(y(1),38.0);
  EXPECT_FLOAT_EQ(y(3),83.0);
  EXPECT_FLOAT_EQ(y(5),128.0);
  EXPECT_FLOAT_EQ(y(7),173.0);
  EXPECT_FLOAT_EQ(y(9),218.0);

  ASSERT_TRUE(A.multiply(u,x,1.0,-1.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,-1.0,true,4,2));

  EXPECT_FLOAT_EQ(x.sum(),0.0);
  EXPECT_FLOAT_EQ(y.sum(),0.0);

  u.resize(5,utl::RETAIN);
  ASSERT_TRUE(A.multiply(u,v));
  v *= 0.5;

  EXPECT_FLOAT_EQ(v(1),67.5);
  EXPECT_FLOAT_EQ(v(2),75.0);
  EXPECT_FLOAT_EQ(v(3),82.5);
}


template<template<class T> class Matrix>
void normTest()
{
  Matrix<double> a(4,5);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  EXPECT_FLOAT_EQ(a.sum(),210.0);
  EXPECT_FLOAT_EQ(a.sum(5),34.0);
  EXPECT_FLOAT_EQ(a.asum(5),34.0);
  EXPECT_FLOAT_EQ(a.trace(),34.0);
  EXPECT_NEAR(a.norm2(5),sqrt(414.0),1.0e-15);
}


/*
TEST(TestMatrix, ExtractBlock)
{
  utl::matrix<int> a(3,3), b(2,2);

  std::iota(a.begin(), a.end(), 1);

  a.extractBlock(b,1,1);
  EXPECT_EQ(b(1,1), 1);
  EXPECT_EQ(b(2,1), 2);
  EXPECT_EQ(b(1,2), 4);
  EXPECT_EQ(b(2,2), 5);

  a.extractBlock(b,2,2,true);
  EXPECT_EQ(b(1,1), 1+5);
  EXPECT_EQ(b(2,1), 2+6);
  EXPECT_EQ(b(1,2), 4+8);
  EXPECT_EQ(b(2,2), 5+9);
}


TEST(TestMatrix, AddRows)
{
  utl::matrix<int> a(3,5);
  std::iota(a.begin(),a.end(),1);
  std::cout <<"A:"<< a;

  a.expandRows(1);
  std::cout <<"B:"<< a;
  int fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++)
  {
    for (size_t i = 1; i <= 3; i++, fasit++)
      EXPECT_EQ(a(i,j), fasit);
    EXPECT_EQ(a(4,j), 0);
  }

  a.expandRows(-2);
  std::cout <<"C:"<< a;
  fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++, fasit++)
    for (size_t i = 1; i <= 2; i++, fasit++)
      EXPECT_EQ(a(i,j), fasit);
}


TEST(TestMatrix, AugmentRows)
{
  utl::matrix<int> a(5,3), b(4,3), c(3,2);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),16);
  std::cout <<"A:"<< a;
  std::cout <<"B:"<< b;
  size_t nA = a.size();
  size_t na = a.rows();
  size_t nb = b.rows();
  ASSERT_TRUE(a.augmentRows(b));
  ASSERT_FALSE(a.augmentRows(c));
  std::cout <<"C:"<< a;
  for (size_t j = 1; j <= a.cols(); j++)
    for (size_t i = 1; i <= a.rows(); i++)
      if (i <= na)
        EXPECT_EQ(a(i,j), i+na*(j-1));
      else
        EXPECT_EQ(a(i,j), nA-na+i+nb*(j-1));
}


TEST(TestMatrix, AugmentCols)
{
  utl::matrix<int> a(3,5), b(3,4), c(2,3);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),16);
  std::cout <<"A:"<< a;
  std::cout <<"B:"<< b;
  ASSERT_TRUE(a.augmentCols(b));
  ASSERT_FALSE(a.augmentCols(c));
  std::cout <<"C:"<< a;
  int fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++)
    for (size_t i = 1; i <= a.rows(); i++, fasit++)
      EXPECT_EQ(a(i,j), fasit);
}


TEST(TestMatrix, SumCols)
{
  utl::matrix<int> a(5,3);
  std::iota(a.begin(),a.end(),1);
  std::cout <<"A:"<< a;
  EXPECT_EQ(a.sum(-1), 15);
  EXPECT_EQ(a.sum(-2), 40);
  EXPECT_EQ(a.sum(-3), 65);
}


TEST(TestMatrix, Multiply)
{
  utl::vector<double> u(14), v(9), x, y;
  utl::matrix<double> A(3,5);

  std::iota(u.begin(),u.end(),1.0);
  std::iota(v.begin(),v.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);

  ASSERT_TRUE(A.multiply(u,x,1.0,0.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,0.0,true,4,2));

  EXPECT_FLOAT_EQ(x(3),370.0);
  EXPECT_FLOAT_EQ(x(7),410.0);
  EXPECT_FLOAT_EQ(x(11),450.0);
  EXPECT_FLOAT_EQ(y(1),38.0);
  EXPECT_FLOAT_EQ(y(3),83.0);
  EXPECT_FLOAT_EQ(y(5),128.0);
  EXPECT_FLOAT_EQ(y(7),173.0);
  EXPECT_FLOAT_EQ(y(9),218.0);

  ASSERT_TRUE(A.multiply(u,x,1.0,-1.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,-1.0,true,4,2));

  EXPECT_FLOAT_EQ(x.sum(),0.0);
  EXPECT_FLOAT_EQ(y.sum(),0.0);

  u.resize(5,utl::RETAIN);
  ASSERT_TRUE(A.multiply(u,v));
  v *= 0.5;

  EXPECT_FLOAT_EQ(v(1),67.5);
  EXPECT_FLOAT_EQ(v(2),75.0);
  EXPECT_FLOAT_EQ(v(3),82.5);
}


TEST(TestMatrix, Norm)
{
  utl::matrix<double> a(4,5);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  EXPECT_FLOAT_EQ(a.sum(),210.0);
  EXPECT_FLOAT_EQ(a.sum(5),34.0);
  EXPECT_FLOAT_EQ(a.asum(5),34.0);
  EXPECT_FLOAT_EQ(a.trace(),34.0);
  EXPECT_NEAR(a.norm2(5),sqrt(414.0),1.0e-15);
}


TEST(TestMatrix, Read)
{
  utl::vector<double> a(26);
  utl::matrix<double> A(2,3);
  std::iota(a.begin(),a.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);
  std::cout <<"a:"<< a;
  std::cout <<"A:"<< A;

  auto&& checkVector = [&a](const char* fname)
  {
    std::ifstream is(fname,std::ios::in);
    utl::vector<double> b;
    is >> b;
    std::cout <<"b:"<< b;
    ASSERT_EQ(a.size(),b.size());
    for (size_t i = 1; i <= a.size(); i++)
      EXPECT_NEAR(a(i), b(i), 1.0e-13);
  };

  auto&& checkMatrix = [&A](const char* fname)
  {
    std::ifstream is(fname,std::ios::in);
    utl::matrix<double> B;
    is >> B;
    std::cout <<"B:"<< B;
    ASSERT_EQ(A.rows(),B.rows());
    ASSERT_EQ(A.cols(),B.cols());
    for (size_t i = 1; i <= A.rows(); i++)
      for (size_t j = 1; j <= A.cols(); j++)
        EXPECT_NEAR(A(i,j), B(i,j), 1.0e-13);
  };

  const char* fname0 = "/tmp/testVector.dat";
  const char* fname1 = "/tmp/testMatrix1.dat";
  const char* fname2 = "/tmp/testMatrix2.dat";
  const char* fname3 = "/tmp/testMatrix3.dat";
  const char* fname4 = "/tmp/testMatrix4.dat";

  std::ofstream os(fname0);
  os << a.size() << a;
  os.close();

  os.open(fname1);
  os << A.rows() <<' '<< A.cols() << A;
  os.close();

  os.open(fname2);
  os << A.rows() <<' '<< A.cols();
  for (size_t i = 1; i <= A.rows(); i++)
    for (size_t j = 1; j <= A.cols(); j++)
      os << (j == 1 ? '\n' : ' ') << A(i,j);
  os <<'\n';
  os.close();

  checkVector(fname0);
  checkMatrix(fname1);
  checkMatrix(fname2);

  double value = 0.0;
  A.resize(6,6);
  for (size_t i = 1; i <= A.rows(); i++)
    for (size_t j = i; j <= A.cols(); j++)
      A(i,j) = A(j,i) = ++value;
  std::cout <<"Symmetric A:"<< A;

  os.open(fname3);
  os <<"Symmetric: "<< A.rows() << A;
  os.close();

  checkMatrix(fname3);

  std::iota(A.begin(),A.end(),1.0);
  std::cout <<"Non-symmetric A:"<< A;
  os.open(fname4);
  os <<"Column-oriented: "<< A.rows() <<" "<< A.cols();
  for (double v : A) os <<"\n"<< v;
  os <<"\n";
  os.close();

  checkMatrix(fname4);
}


TEST(TestMatrix3D, Trace)
{
  utl::matrix3d<double> a(4,3,3);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  for (size_t i = 1; i <= 4; i++)
    EXPECT_FLOAT_EQ(a.trace(i),3.0*i+48.0);
}


TEST(TestMatrix3D, GetColumn)
{
  utl::matrix3d<int> A(4,3,2);
  std::iota(A.begin(),A.end(),1);
  std::cout <<"A:"<< A;

  int value = 1;
  for (size_t c = 1; c <= A.dim(3); c++)
    for (size_t r = 1; r <= A.dim(2); r++)
    {
      utl::vector<int> column = A.getColumn(r,c);
      EXPECT_EQ(column.size(), A.dim(1));
      for (size_t i = 0; i < column.size(); i++, value++)
        EXPECT_EQ(value, column[i]);
    }

  EXPECT_EQ(value, 1+(int)A.size());
}


TEST(TestMatrix3D, DumpRead)
{
  int i = 0;
  utl::matrix3d<double> A(2,3,4);
  for (double& v : A) v = 3.14159*(++i);

  const char* fname = "/tmp/testMatrix3D.dat";
  std::ofstream os(fname);
  os << std::setprecision(16) << A;
  std::ifstream is(fname,std::ios::in);
  utl::matrix3d<double> B(is);
  B -= A;
  ASSERT_NEAR(B.norm2(), 0.0, 1.0e-13);
}


TEST(TestMatrix3D, Multiply)
{
  std::vector<double> a(10);
  utl::matrix<double> A(2,5);
  utl::matrix3d<double> B(5,4,3), C, D;

  std::iota(a.begin(),a.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);
  std::iota(B.begin(),B.end(),1.0);

  C.multiply(A,B);
  ASSERT_TRUE(D.multiplyMat(a,B));

  std::vector<double>::const_iterator c = C.begin();
  for (double d : D)
    EXPECT_EQ(d,*(c++));
}
*/

}
