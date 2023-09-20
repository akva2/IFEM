//==============================================================================
//!
//! \file TestFunctions.C
//!
//! \date Apr 28 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for parsing of functions.
//!
//==============================================================================

#include "ExprFunctions.h"
#include "Functions.h"
#include <cstdlib>
#include <cmath>

#include "gtest/gtest.h"


TEST(TestScalarFunc, ParseDerivative)
{
  const char* func1 = "sin(1.5*t)*t";
  const char* func2 = "sin(1.5*t)*t:1.5*cos(1.5*t)*t+sin(1.5*t)";

  ScalarFunc* f1 = utl::parseTimeFunc(func1,"expression");
  ScalarFunc* f2 = utl::parseTimeFunc(func2,"expression");

  ASSERT_TRUE(f1 != nullptr);
  ASSERT_TRUE(f2 != nullptr);

  EXPECT_FALSE(f1->isConstant());
  EXPECT_FALSE(f2->isConstant());

  double t = 0.0;
  for (int i = 0; i < 20; i++)
  {
    t += 0.314*(double)random()/(double)RAND_MAX;
    std::cout <<"f("<< t <<") = "<< (*f1)(t)
	      <<"  f'("<< t <<") = "<< f1->deriv(t) << std::endl;
    EXPECT_FLOAT_EQ((*f1)(t),sin(1.5*t)*t);
    EXPECT_FLOAT_EQ((*f2)(t),sin(1.5*t)*t);
    EXPECT_FLOAT_EQ(f1->deriv(t),1.5*cos(1.5*t)*t+sin(1.5*t));
    EXPECT_FLOAT_EQ(f2->deriv(t),1.5*cos(1.5*t)*t+sin(1.5*t));
  }
}


TEST(TestRealFunc, GradientFD)
{
  const char* f1 = "sin(x)*sin(y)*sin(z)";

  const double eps = 1e-6;

  EvalFunction f(f1, eps);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 r((sin(x+0.5*eps)-sin(x-0.5*eps))*sin(y)*sin(z) / eps,
                     sin(x)*(sin(y+0.5*eps)-sin(y-0.5*eps))*sin(z) / eps,
                     sin(x)*sin(y)*(sin(z+0.5*eps)-sin(z-0.5*eps)) / eps);
        EXPECT_NEAR(f.deriv(X, 1), r[0], 1e-8);
        EXPECT_NEAR(f.deriv(X, 2), r[1], 1e-8);
        EXPECT_NEAR(f.deriv(X, 3), r[2], 1e-8);
        const Vec3 grad = f.gradient(X);
        EXPECT_NEAR(grad[0], r[0], 1e-8);
        EXPECT_NEAR(grad[1], r[1], 1e-8);
        EXPECT_NEAR(grad[2], r[2], 1e-8);
      }
}


TEST(TestRealFunc, Gradient)
{
  const char* f1 = "sin(x)*sin(y)*sin(z)";
  const char* d1 = "cos(x)*sin(y)*sin(z)";
  const char* d2 = "sin(x)*cos(y)*sin(z)";
  const char* d3 = "sin(x)*sin(y)*cos(z)";

  EvalFunction f(f1);
  f.addDerivative(d1, "", 1);
  f.addDerivative(d2, "", 2);
  f.addDerivative(d3, "", 3);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 r(cos(x)*sin(y)*sin(z),
                     sin(x)*cos(y)*sin(z),
                     sin(x)*sin(y)*cos(z));
        EXPECT_DOUBLE_EQ(f.deriv(X, 1), r[0]);
        EXPECT_DOUBLE_EQ(f.deriv(X, 2), r[1]);
        EXPECT_DOUBLE_EQ(f.deriv(X, 3), r[2]);
        const Vec3 grad = f.gradient(X);
        EXPECT_DOUBLE_EQ(grad[0], r[0]);
        EXPECT_DOUBLE_EQ(grad[1], r[1]);
        EXPECT_DOUBLE_EQ(grad[2], r[2]);
      }
}


TEST(TestEvalFunction, ExtraParam)
{
  EvalFunction f("x*foo");
  f.setParam("foo", 2.0);
  Vec3 X(1.0,0.0,0.0);
  EXPECT_FLOAT_EQ(f(X), 2.0);
  X.x = 0.5;
  f.setParam("foo", 4.0);
  EXPECT_FLOAT_EQ(f(X), 2.0);
}


TEST(TestEvalFunction, isConstant)
{
  EXPECT_TRUE (EvalFunction("2.0*x*y").isConstant());
  EXPECT_FALSE(EvalFunction("2.0*x*t").isConstant());
  EXPECT_TRUE (EvalFunction("1.8*tan(x)*x").isConstant());
  EXPECT_FALSE(EvalFunction("2.0*x*tan(t*3)+y").isConstant());
}
