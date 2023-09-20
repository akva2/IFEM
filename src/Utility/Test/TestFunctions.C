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
#include "autodiff/reverse/var/var.hpp"
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


TEST(TestRealFunc, Hessian)
{
  const char* f1 = "sin(x)*sin(y)*sin(z)";
  const char* d11 = "-sin(x)*sin(y)*sin(z)";
  const char* d22 = "-sin(x)*sin(y)*sin(z)";
  const char* d33 = "-sin(x)*sin(y)*sin(z)";
  const char* d12 = "cos(x)*cos(y)*sin(z)";
  const char* d13 = "cos(x)*sin(y)*cos(z)";
  const char* d23 = "sin(x)*cos(y)*cos(z)";

  EvalFunction f(f1);
  f.addDerivative(d11,"",1,1);
  f.addDerivative(d12,"",1,2);
  f.addDerivative(d13,"",1,3);
  f.addDerivative(d22,"",2,2);
  f.addDerivative(d23,"",2,3);
  f.addDerivative(d33,"",3,3);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1})//, 0.2, 0.3})
    for (double y : {0.1})//, 0.6, 0.7})
      for (double z : {0.1}){//, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const std::array<double,6> r{-sin(x)*sin(y)*sin(z),
                                      -sin(x)*sin(y)*sin(z)
                                      -sin(x)*sin(y)*sin(z),
                                      cos(x)*cos(y)*sin(z),
                                      cos(x)*sin(y)*cos(z),
                                      sin(x)*cos(y)*cos(z)};
        EXPECT_DOUBLE_EQ(f.dderiv(X,1,1), r[0]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,2,2), r[1]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,3,3), r[2]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,1,2), r[3]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,2,1), r[3]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,1,3), r[4]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,3,1), r[4]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,2,3), r[5]);
        EXPECT_DOUBLE_EQ(f.dderiv(X,3,2), r[5]);
        const SymmTensor grad = f.hessian(X);
        EXPECT_DOUBLE_EQ(grad(1,1), r[0]);
        EXPECT_DOUBLE_EQ(grad(2,2), r[1]);
        EXPECT_DOUBLE_EQ(grad(3,3), r[2]);
        EXPECT_DOUBLE_EQ(grad(1,2), r[3]);
        EXPECT_DOUBLE_EQ(grad(1,3), r[4]);
        EXPECT_DOUBLE_EQ(grad(2,3), r[5]);
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

TEST(TestEvalFunc, AutoDiff)
{
  EvalFuncImpl<autodiff::var> func("sin(x)");
  EXPECT_DOUBLE_EQ(func(M_PI / 2.0), 1.0);
  EXPECT_DOUBLE_EQ(func.deriv(0.0), 1.0);
}


TEST(TestEvalFunction, AutoDiff)
{
  EvalFunctionImpl<autodiff::var> func("sin(x)*sin(y)*sin(z)*sin(2*t)");
  Vec4 X;
  X.x = 1.0;
  X.y = 2.0;
  X.z = 3.0;
  X.t = 4.0;

  EXPECT_DOUBLE_EQ(func(X),          sin(1.0) * sin(2.0) * sin(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.deriv(X, 1), cos(1.0) * sin(2.0) * sin(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.deriv(X, 2), sin(1.0) * cos(2.0) * sin(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.deriv(X, 3), sin(1.0) * sin(2.0) * cos(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.deriv(X, 4), sin(1.0) * sin(2.0) * sin(3.0) * 2.0 * cos(2.0 * 4.0));

  EXPECT_DOUBLE_EQ(func.dderiv(X, 1, 1), -sin(1.0) * sin(2.0) * sin(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.dderiv(X, 1, 2), cos(1.0) * cos(2.0) * sin(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.dderiv(X, 1, 3), cos(1.0) * sin(2.0) * cos(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.dderiv(X, 2, 2), sin(1.0) * -sin(2.0) * sin(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.dderiv(X, 2, 3), sin(1.0) * cos(2.0) * cos(3.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(func.dderiv(X, 3, 3), sin(1.0) * sin(2.0) * -sin(3.0) * sin(2.0 * 4.0));
}


TEST(TestVecFuncExpr, AutoDiff2D)
{
  EvalMultiFunction<VecFunc,Vec3,autodiff::var> func("sin(2*x)*sin(y)*sin(t) |"
                                                     "sin(x)*sin(2*y)*sin(t)");
  Vec4 X;
  X.x = 1.0;
  X.y = 2.0;
  X.t = 4.0;

  const auto val = func(X);

  EXPECT_DOUBLE_EQ(val[0],          sin(2.0 * 1.0) * sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val[1],          sin(1.0) * sin(2.0 * 2.0) * sin(4.0));

  const auto grad = func.gradient(X);
  EXPECT_DOUBLE_EQ(grad(1,1), 2.0 * cos(2.0 * 1.0) * sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,1), cos(1.0) * sin(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,2), sin(2.0 * 1.0) * cos(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,2), sin(1.0) * 2.0 * cos(2.0 * 2.0) * sin(4.0));

  const auto ut = func.tgradient(X);
  EXPECT_DOUBLE_EQ(ut[0], sin(2.0 * 1.0) * sin(2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut[1], sin(1.0) * sin(2.0 * 2.0) * cos(4.0));

  const auto hess = func.hessian(X);
  EXPECT_DOUBLE_EQ(hess(1,1,1), -4.0 * sin(2.0 * 1.0) * sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,1), -sin(1.0) * sin(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1), 2.0 * cos(2.0 * 1.0) * cos(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1), hess(1,1,2));
  EXPECT_DOUBLE_EQ(hess(2,1,2), cos(1.0) * 2.0 * cos(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,1), hess(2,1,2));
  EXPECT_DOUBLE_EQ(hess(1,2,2), sin(2.0 * 1.0) * -sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,2), sin(1.0) * -4.0 * sin(2.0 * 2.0) * sin(4.0));
}


TEST(TestVecFuncExpr, AutoDiff3D)
{
  EvalMultiFunction<VecFunc,Vec3,autodiff::var> func("sin(2*x)*sin(y)*sin(z)*sin(t) |"
                                                     "sin(x)*sin(2*y)*sin(z)*sin(t) |"
                                                     "sin(x)*sin(y)*sin(2*z)*sin(t)");
  Vec4 X;
  X.x = 1.0;
  X.y = 2.0;
  X.z = 3.0;
  X.t = 4.0;

  const auto val = func(X);

  EXPECT_DOUBLE_EQ(val[0],          sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val[1],          sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val[2],          sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));

  const auto grad = func.gradient(X);
  EXPECT_DOUBLE_EQ(grad(1,1), 2.0 * cos(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,1), cos(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,1), cos(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,2), sin(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,2), sin(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,2), sin(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,3), sin(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,3), sin(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,3), sin(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));

  const auto ut = func.tgradient(X);
  EXPECT_DOUBLE_EQ(ut[0], sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut[1], sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut[2], sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * cos(4.0));

  const auto hess = func.hessian(X);
  EXPECT_DOUBLE_EQ(hess(1,1,1), -4.0 * sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,1), -sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,1,1), -sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1), 2.0 * cos(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1), hess(1,1,2));
  EXPECT_DOUBLE_EQ(hess(2,1,2), cos(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,1), hess(2,1,2));
  EXPECT_DOUBLE_EQ(hess(3,1,2), cos(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,1), hess(3,1,2));
  EXPECT_DOUBLE_EQ(hess(1,1,3), 2.0 * cos(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,1,3), hess(1,3,1));
  EXPECT_DOUBLE_EQ(hess(2,1,3), cos(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,3), hess(2,3,1));
  EXPECT_DOUBLE_EQ(hess(3,1,3), cos(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,1,3), hess(3,3,1));
  EXPECT_DOUBLE_EQ(hess(1,2,2), sin(2.0 * 1.0) * -sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,2), sin(1.0) * -4.0 * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,2), sin(1.0) * -sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,3), sin(2.0 * 1.0) * cos(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,3), hess(1,3,2));
  EXPECT_DOUBLE_EQ(hess(2,2,3), sin(1.0) * 2.0 * cos(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,3), hess(2,3,2));
  EXPECT_DOUBLE_EQ(hess(3,2,3), sin(1.0) * cos(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,3), hess(3,3,2));
  EXPECT_DOUBLE_EQ(hess(1,3,3), sin(2.0 * 1.0) * sin(2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,3,3), sin(1.0) * sin(2.0 * 2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,3,3), sin(1.0) * sin(2.0) * -4.0 * sin(2.0 * 3.0) * sin(4.0));
}


TEST(TestTensorFuncExpr, AutoDiff2D)
{
  EvalMultiFunction<TensorFunc,Tensor,autodiff::var> func("sin(2*x)*sin(y)*sin(t) |"
                                                          "sin(x)*sin(2*y)*sin(t) |"
                                                          "sin(x)*sin(2*y)*cos(t) |"
                                                          "sin(x)*sin(y)*sin(2*t)");
  Vec4 X;
  X.x = 1.0;
  X.y = 2.0;
  X.z = 0.0;
  X.t = 4.0;

  const auto val = func(X);
  EXPECT_DOUBLE_EQ(val(1,1), sin(2.0 * 1.0) * sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(1,2), sin(1.0) * sin(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(2,1), sin(1.0) * sin(2.0 * 2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(val(2,2), sin(1.0) * sin(2.0) * sin(2.0 * 4.0));

  const auto grad = func.gradient(X);
  EXPECT_DOUBLE_EQ(grad(1,1,1), 2.0 * cos(2.0 * 1.0) * sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,1,1), cos(1.0) * sin(2.0 * 2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(grad(1,2,1), cos(1.0) * sin(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,2,1), cos(1.0) * sin(2.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(grad(1,1,2), sin(2.0 * 1.0) * cos(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,1,2), sin(1.0) * 2.0 * cos(2.0 * 2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(grad(1,2,2), sin(1.0) * 2.0 * cos(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,2,2), sin(1.0) * cos(2.0) * sin(2.0 * 4.0));

  const auto ut = func.tgradient(X);
  EXPECT_DOUBLE_EQ(ut(1,1), sin(2.0 * 1.0) * sin(2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(1,2), sin(1.0) * sin(2.0 * 2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(2,1), sin(1.0) * sin(2.0 * 2.0) * -sin(4.0));
  EXPECT_DOUBLE_EQ(ut(2,2), sin(1.0) * sin(2.0) * 2.0 * cos(2.0 * 4.0));

  const auto hess = func.hessian(X);
  EXPECT_DOUBLE_EQ(hess(1,1,1,1), -4.0 * sin(2.0 * 1.0) * sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1,1), -sin(1.0) * sin(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,1,1), -sin(1.0) * sin(2.0 * 2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,1,1), -sin(1.0) * sin(2.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(hess(1,1,1,2), 2.0 * cos(2.0 * 1.0) * cos(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,1,1,2), hess(1,1,2,1));
  EXPECT_DOUBLE_EQ(hess(1,2,1,2), cos(1.0) * 2.0 * cos(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1,2), hess(1,2,2,1));
  EXPECT_DOUBLE_EQ(hess(2,1,1,2), hess(2,1,2,1));
  EXPECT_DOUBLE_EQ(hess(2,2,1,2), cos(1.0) * cos(2.0) * sin(2.0 * 4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,1,2), hess(2,2,2,1));
  EXPECT_DOUBLE_EQ(hess(1,1,2,2), sin(2.0 * 1.0) * -sin(2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,2,2), sin(1.0) * -4.0 * sin(2.0 * 2.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,2,2), sin(1.0) * -4.0 * sin(2.0 * 2.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,2,2), sin(1.0) * -sin(2.0) * sin(2.0 * 4.0));
}


TEST(TestTensorFuncExpr, AutoDiff3D)
{
  EvalMultiFunction<TensorFunc,Tensor,autodiff::var> func("sin(2*x)*sin(y)*sin(z)*sin(t) |"
                                                          "sin(x)*sin(2*y)*sin(z)*sin(t) |"
                                                          "sin(x)*sin(y)*sin(2*z)*sin(t) |"
                                                          "sin(2*x)*sin(y)*sin(z)*sin(t) |"
                                                          "sin(x)*sin(2*y)*sin(z)*sin(t) |"
                                                          "sin(x)*sin(y)*sin(2*z)*sin(t) |"
                                                          "sin(2*x)*sin(y)*sin(z)*sin(t) |"
                                                          "sin(x)*sin(2*y)*sin(z)*sin(t) |"
                                                          "sin(x)*sin(y)*sin(2*z)*sin(t)");
  Vec4 X;
  X.x = 1.0;
  X.y = 2.0;
  X.z = 3.0;
  X.t = 4.0;

  const auto val = func(X);
  EXPECT_DOUBLE_EQ(val(1,1), sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(1,2), sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(1,3), sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(2,1), sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(2,2), sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(2,3), sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(3,1), sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(3,2), sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(val(3,3), sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));

  const auto grad = func.gradient(X);
  EXPECT_DOUBLE_EQ(grad(1,1,1), 2.0 * cos(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,2,1), cos(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,3,1), cos(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,1,1), 2.0 * cos(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,2,1), cos(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,3,1), cos(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,1,1), 2.0 * cos(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,2,1), cos(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,3,1), cos(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));

  EXPECT_DOUBLE_EQ(grad(1,1,2), sin(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,2,2), sin(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,3,2), sin(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,1,2), sin(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,2,2), sin(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,3,2), sin(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,1,2), sin(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,2,2), sin(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,3,2), sin(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));

  EXPECT_DOUBLE_EQ(grad(1,1,3), sin(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,2,3), sin(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(1,3,3), sin(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,1,3), sin(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,2,3), sin(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(2,3,3), sin(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,1,3), sin(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,2,3), sin(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(grad(3,3,3), sin(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));

  const auto ut = func.tgradient(X);
  EXPECT_DOUBLE_EQ(ut(1,1), sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(1,2), sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(1,3), sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(2,1), sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(2,2), sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(2,3), sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(3,1), sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(3,2), sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * cos(4.0));
  EXPECT_DOUBLE_EQ(ut(3,3), sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * cos(4.0));

  const auto hess = func.hessian(X);
  EXPECT_DOUBLE_EQ(hess(1,1,1,1), -4.0 * sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1,1), -sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,3,1,1), -sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,1,1), -4.0 * sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,1,1), -sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,3,1,1), -sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,1,1,1), -4.0 * sin(2.0 * 1.0) * sin(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,1,1), -sin(1.0) * sin(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,3,1,1), -sin(1.0) * sin(2.0) * sin(2.0 * 3.0) * sin(4.0));

  EXPECT_DOUBLE_EQ(hess(1,1,1,2), 2.0 * cos(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,1,1,2), hess(1,1,2,1));
  EXPECT_DOUBLE_EQ(hess(1,2,1,2), cos(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1,2), hess(1,2,2,1));
  EXPECT_DOUBLE_EQ(hess(1,3,1,2), cos(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,3,1,2), hess(1,3,2,1));
  EXPECT_DOUBLE_EQ(hess(2,1,1,2), 2.0 * cos(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,1,2), hess(2,1,2,1));
  EXPECT_DOUBLE_EQ(hess(2,2,1,2), cos(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,1,2), hess(2,2,2,1));
  EXPECT_DOUBLE_EQ(hess(2,3,1,2), cos(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,3,1,2), hess(2,3,2,1));
  EXPECT_DOUBLE_EQ(hess(3,1,1,2), 2.0 * cos(2.0 * 1.0) * cos(2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,1,1,2), hess(3,1,2,1));
  EXPECT_DOUBLE_EQ(hess(3,2,1,2), cos(1.0) * 2.0 * cos(2.0 * 2.0) * sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,1,2), hess(3,2,2,1));
  EXPECT_DOUBLE_EQ(hess(3,3,1,2), cos(1.0) * cos(2.0) * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,3,1,2), hess(3,3,2,1));

  EXPECT_DOUBLE_EQ(hess(1,1,1,3), 2.0 * cos(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,1,1,3), hess(1,1,3,1));
  EXPECT_DOUBLE_EQ(hess(1,2,1,3), cos(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,1,3), hess(1,2,3,1));
  EXPECT_DOUBLE_EQ(hess(1,3,1,3), cos(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,3,1,3), hess(1,3,3,1));
  EXPECT_DOUBLE_EQ(hess(2,1,1,3), 2.0 * cos(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,1,3), hess(2,1,3,1));
  EXPECT_DOUBLE_EQ(hess(2,2,1,3), cos(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,1,3), hess(2,2,3,1));
  EXPECT_DOUBLE_EQ(hess(2,3,1,3), cos(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,3,1,3), hess(2,3,3,1));
  EXPECT_DOUBLE_EQ(hess(3,1,1,3), 2.0 * cos(2.0 * 1.0) * sin(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,1,1,3), hess(3,1,3,1));
  EXPECT_DOUBLE_EQ(hess(3,2,1,3), cos(1.0) * sin(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,1,3), hess(3,2,3,1));
  EXPECT_DOUBLE_EQ(hess(3,3,1,3), cos(1.0) * sin(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,3,1,3), hess(3,3,3,1));

  EXPECT_DOUBLE_EQ(hess(1,1,2,3), sin(2.0 * 1.0) * cos(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,1,2,3), hess(1,1,3,2));
  EXPECT_DOUBLE_EQ(hess(1,2,2,3), sin(1.0) * 2.0 * cos(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,2,3), hess(1,2,3,2));
  EXPECT_DOUBLE_EQ(hess(1,3,2,3), sin(1.0) * cos(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,3,2,3), hess(1,3,3,2));
  EXPECT_DOUBLE_EQ(hess(2,1,2,3), sin(2.0 * 1.0) * cos(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,2,3), hess(2,1,3,2));
  EXPECT_DOUBLE_EQ(hess(2,2,2,3), sin(1.0) * 2.0 * cos(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,2,3), hess(2,2,3,2));
  EXPECT_DOUBLE_EQ(hess(2,3,2,3), sin(1.0) * cos(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,3,2,3), hess(2,3,3,2));
  EXPECT_DOUBLE_EQ(hess(3,1,2,3), sin(2.0 * 1.0) * cos(2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,1,2,3), hess(3,1,3,2));
  EXPECT_DOUBLE_EQ(hess(3,2,2,3), sin(1.0) * 2.0 * cos(2.0 * 2.0) * cos(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,2,3), hess(3,2,3,2));
  EXPECT_DOUBLE_EQ(hess(3,3,2,3), sin(1.0) * cos(2.0) * 2.0 * cos(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,3,2,3), hess(3,3,3,2));

  EXPECT_DOUBLE_EQ(hess(1,1,3,3), sin(2.0 * 1.0) * sin(2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,2,3,3), sin(1.0) * sin(2.0 * 2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(1,3,3,3), sin(1.0) * sin(2.0) * -4.0 * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,1,3,3), sin(2.0 * 1.0) * sin(2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,2,3,3), sin(1.0) * sin(2.0 * 2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(2,3,3,3), sin(1.0) * sin(2.0) * -4.0 * sin(2.0 * 3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,1,3,3), sin(2.0 * 1.0) * sin(2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,2,3,3), sin(1.0) * sin(2.0 * 2.0) * -sin(3.0) * sin(4.0));
  EXPECT_DOUBLE_EQ(hess(3,3,3,3), sin(1.0) * sin(2.0) * -4.0 * sin(2.0 * 3.0) * sin(4.0));
}
