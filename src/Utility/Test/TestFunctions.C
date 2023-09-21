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
#include "TensorFunction.h"
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
        const Vec3 Xp(x+0.5*eps, y+0.5*eps, z+0.5*eps);
        const Vec3 Xm(x-0.5*eps, y-0.5*eps, z-0.5*eps);
        const Vec3 r((sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z) / eps,
                     sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z) / eps,
                     sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)) / eps);

        const Vec3 grad = f.gradient(X);
        for (size_t i = 1; i <= 3; ++i) {
          EXPECT_NEAR(f.deriv(X, i), r[i-1], 1e-8);
          EXPECT_NEAR(grad[i-1], r[i-1], 1e-8);
        }
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

        const Vec3 grad = f.gradient(X);
        for (size_t i = 1; i <= 3; ++i) {
          EXPECT_DOUBLE_EQ(f.deriv(X, i), r[i-1]);
          EXPECT_DOUBLE_EQ(grad[i-1], r[i-1]);
        }
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

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const SymmTensor r({-sin(x)*sin(y)*sin(z),
                            -sin(x)*sin(y)*sin(z),
                            -sin(x)*sin(y)*sin(z),
                             cos(x)*cos(y)*sin(z),
                             sin(x)*cos(y)*cos(z),
                             cos(x)*sin(y)*cos(z)});

        const SymmTensor grad = f.hessian(X);
        for (size_t i = 1; i <= 3; ++i)
          for (size_t j = 1; j <= 3; ++j) {
            EXPECT_DOUBLE_EQ(f.dderiv(X,i,j), r(i,j));
            EXPECT_DOUBLE_EQ(grad(i,j), r(i,j));
          }
      }
}


TEST(TestVecFunction, Gradient2D)
{
  const char* g   = "sin(x)*sin(y) | x*x*y*y";
  const char* g_x = "cos(x)*sin(y) | 2*x*y*y";
  const char* g_y = "sin(x)*cos(y) | 2*x*x*y";

  VecFuncExpr f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Tensor r({cos(x)*sin(y), 2*x*y*y, sin(x)*cos(y), 2*x*x*y});

      const Tensor grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const Vec3 dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i) {
          EXPECT_DOUBLE_EQ(dx[i-1], r(i,d));
          EXPECT_DOUBLE_EQ(grad(i,d), r(i,d));
        }
      }
    }
}


TEST(TestVecFunction, Gradient2DFD)
{
  const char* g = "sin(x)*sin(y) | x*x*y*y";

  const double eps = 1e-6;

  VecFuncExpr f(g,"",eps);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Vec3 Xp(x + 0.5*eps, y + 0.5*eps);
      const Vec3 Xm(x - 0.5*eps, y - 0.5*eps);
      const Tensor r({(sin(Xp.x) - sin(Xm.x))*sin(y) / eps,
                      (Xp.x*Xp.x - Xm.x*Xm.x)*y*y / eps,
                      sin(x)*(sin(Xp.y) - sin(Xm.y)) / eps,
                      x*x*(Xp.y*Xp.y - Xm.y*Xm.y) / eps});

      const Tensor grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const Vec3 dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i) {
          EXPECT_NEAR(dx[i-1], r(i,d), 1e-8);
          EXPECT_NEAR(grad(i,d), r(i,d), 1e-8);
        }
      }
    }
}


TEST(TestVecFunction, Gradient3D)
{
  const char* g   = "sin(x)*sin(y)*sin(z) | x*x*y*y*z*z*z | exp(-x)*exp(2*y)*exp(z*z)";
  const char* g_x = "cos(x)*sin(y)*sin(z) | 2*x*y*y*z*z*z | -exp(-x)*exp(2*y)*exp(z*z)";
  const char* g_y = "sin(x)*cos(y)*sin(z) | 2*x*x*y*z*z*z | 2*exp(-x)*exp(2*y)*exp(z*z)";
  const char* g_z = "sin(x)*sin(y)*cos(z) | 3*x*x*y*y*z*z | 2*z*exp(-x)*exp(2*y)*exp(z*z)";

  VecFuncExpr f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);
  f.addDerivative(g_z,"",3);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Tensor r({cos(x)*sin(y)*sin(z), 2*x*y*y*z*z*z, -exp(-x)*exp(2*y)*exp(z*z),
                        sin(x)*cos(y)*sin(z), 2*x*x*y*z*z*z, 2*exp(-x)*exp(2*y)*exp(z*z),
                        sin(x)*sin(y)*cos(z), 3*x*x*y*y*z*z, 2*z*exp(-x)*exp(2*y)*exp(z*z)});

        const Tensor grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const Vec3 dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i) {
            EXPECT_DOUBLE_EQ(dx[i-1], r(i,d));
            EXPECT_DOUBLE_EQ(grad(i,d), r(i,d));
          }
        }
      }
}


TEST(TestVecFunction, Gradient3DFD)
{
  const char* g   = "sin(x)*sin(y)*sin(z) | x*x*y*y*z*z*z | exp(-x)*exp(2*y)*exp(z*z)";

  const double eps = 1e-6;
  VecFuncExpr f(g,"",eps);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 Xp(x + 0.5*eps, y + 0.5*eps, z + 0.5*eps);
        const Vec3 Xm(x - 0.5*eps, y - 0.5*eps, z - 0.5*eps);
        const Tensor r({(sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z) / eps,
                        (Xp.x*Xp.x - Xm.x*Xm.x)*X.y*X.y*X.z*X.z*X.z / eps,
                        (exp(-Xp.x) - exp(-Xm.x))*exp(2*y)*exp(z*z) / eps,
                        sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z) / eps,
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y)*z*z*z / eps,
                        exp(-x)*(exp(2*Xp.y) - exp(2*Xm.y))*exp(z*z) / eps,
                        sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)) / eps,
                        x*x*y*y*(Xp.z*Xp.z*Xp.z - Xm.z*Xm.z*Xm.z) / eps,
                        exp(-x)*exp(2*y)*(exp(Xp.z*Xp.z) - exp(Xm.z*Xm.z)) / eps});

        const Tensor grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const Vec3 dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i) {
            EXPECT_NEAR(dx[i-1], r(i,d), 1e-8);
            EXPECT_NEAR(grad(i,d), r(i,d), 1e-8);
          }
        }
      }
}


TEST(TestTensorFunction, Gradient2D)
{
  const char* g   = "sin(x)*sin(y) | x*x*y*y | exp(x)*exp(2*y)   | exp(-2*x)*exp(y)";
  const char* g_x = "cos(x)*sin(y) | 2*x*y*y | exp(x)*exp(2*y)   | -2*exp(-2*x)*exp(y)";
  const char* g_y = "sin(x)*cos(y) | 2*x*x*y | 2*exp(x)*exp(2*y) | exp(-2*x)*exp(y)";

  TensorFuncExpr f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{cos(x)*sin(y), 2*x*y*y, exp(x)*exp(2*y), -2*exp(-2*x)*exp(y),
                        sin(x)*cos(y), 2*x*x*y, 2*exp(x)*exp(2*y), exp(-2*x)*exp(y)}.data());

      const utl::matrix3d<Real> grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const Tensor dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j) {
            EXPECT_DOUBLE_EQ(dx(i,j), r(i,j,d));
            EXPECT_DOUBLE_EQ(grad(i,j,d), r(i,j,d));
          }
      }
    }
}


TEST(TestTensorFunction, Gradient2DFD)
{
  const char* g   = "sin(x)*sin(y) | x*x*y*y | exp(x)*exp(2*y) | exp(-2*x)*exp(y)";

  const double eps = 1e-6;
  TensorFuncExpr f(g,"",eps);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Vec3 Xp(x + 0.5*eps, y + 0.5*eps);
      const Vec3 Xm(x - 0.5*eps, y - 0.5*eps);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y),
                        (Xp.x*Xp.x - Xm.x*Xm.x)*y*y,
                        (exp(Xp.x) - exp(Xm.x))*exp(2*y),
                        (exp(-2*Xp.x) - exp(-2*Xm.x))*exp(y),
                        sin(x)*(sin(Xp.y) - sin(Xm.y)),
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y),
                        exp(x)*(exp(2*Xp.y) - exp(2*Xm.y)),
                        exp(-2*x)*(exp(Xp.y) - exp(Xm.y))}.data());
      r.multiply(1.0 / eps);

      const utl::matrix3d<Real> grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const Tensor dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j) {
            EXPECT_NEAR(dx(i,j), r(i,j,d), 1e-8);
            EXPECT_NEAR(grad(i,j,d), r(i,j,d), 1e-8);
          }
      }
    }
}


TEST(TestTensorFunction, Gradient3D)
{
  const char* g   = "sin(x)*sin(y)*sin(z)  | x*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                    "exp(-2*x)*exp(y)*z    | x*y*z     | x*y*z*z |"
                    "x                     | y         | z";
  const char* g_x = "cos(x)*sin(y)*sin(z)  | 2*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                    "-2*exp(-2*x)*exp(y)*z | y*z       | y*z*z |"
                    "1.0                   | 0.0       | 0.0";
  const char* g_y = "sin(x)*cos(y)*sin(z)  | 2*x*x*y*z | 2*exp(x)*exp(2*y)*z*z |"
                    "exp(-2*x)*exp(y)*z    | x*z       | x*z*z |"
                    "0.0                   | 1.0       | 0.0";
  const char* g_z = "sin(x)*sin(y)*cos(z)  | x*x*y*y   | exp(x)*exp(2*y)*2*z |"
                    "exp(-2*x)*exp(y)      | y*z       | x*y*2*z |"
                    "0.0                   | 0.0       | 1.0";

  TensorFuncExpr f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);
  f.addDerivative(g_z,"",3);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{cos(x)*sin(y)*sin(z),  2*x*y*y*z, exp(x)*exp(2*y)*z*z,
                          -2*exp(-2*x)*exp(y)*z, y*z,       y*z*z,
                          1.0,                   0.0,       0.0,
                          sin(x)*cos(y)*sin(z),  2*x*x*y*z, 2*exp(x)*exp(2*y)*z*z,
                          exp(-2*x)*exp(y)*z,    x*z,       x*z*z,
                          0.0,                   1.0,       0.0,
                          sin(x)*sin(y)*cos(z),  x*x*y*y,   exp(x)*exp(2*y)*2*z,
                          exp(-2*x)*exp(y),      y*z,       x*y*2*z,
                          0.0,                   0.0,       1.0}.data());

        const utl::matrix3d<Real> grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const Tensor dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j) {
              EXPECT_DOUBLE_EQ(dx(i,j), r(i,j,d));
              EXPECT_DOUBLE_EQ(grad(i,j,d), r(i,j,d));
            }
        }
    }
}


TEST(TestTensorFunction, Gradient3DFD)
{
  const char* g   = "sin(x)*sin(y)*sin(z)  | x*x*y*y*z | exp(x)*exp(2*y)*z*z |"
                    "exp(-2*x)*exp(y)*z    | x*y*z     | x*y*z*z |"
                    "x                     | y         | z";

  const double eps = 1e-6;
  TensorFuncExpr f(g,"",eps);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 Xp(x + 0.5*eps, y + 0.5*eps, z + 0.5*eps);
        const Vec3 Xm(x - 0.5*eps, y - 0.5*eps, z - 0.5*eps);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z),
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*y*z,
                          (exp(Xp.x) - exp(Xm.x))*exp(2*y)*z*z,
                          (exp(-2*Xp.x) - exp(-2*Xm.x))*exp(y)*z,
                          (Xp.x - Xm.x)*y*z,
                          (Xp.x - Xm.x)*y*z*z,
                          Xp.x -Xm.x,
                          0.0,
                          0.0,
                          sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z),
                          x*x*(Xp.y*Xp.y - Xm.y*Xm.y)*z,
                          exp(x)*(exp(2*Xp.y) - exp(2*Xm.y))*z*z,
                          exp(-2*x)*(exp(Xp.y) - exp(Xm.y))*z,
                          x*(Xp.y - Xm.y)*z,
                          x*(Xp.y - Xm.y)*z*z,
                          0.0,
                          Xp.y - Xm.y,
                          0.0,
                          sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)),
                          x*x*y*y*(Xp.z - Xm.z),
                          exp(x)*exp(2*y)*(Xp.z*Xp.z - Xm.z*Xm.z),
                          exp(-2*x)*exp(y)*(Xp.z - Xm.z),
                          x*y*(Xp.z - Xm.z),
                          x*y*(Xp.z*Xp.z - Xm.z*Xm.z),
                          0.0,
                          0.0,
                          Xp.z - Xm.z}.data());
        r.multiply(1.0 / eps);

        const utl::matrix3d<Real> grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const Tensor dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j) {
              EXPECT_NEAR(dx(i,j), r(i,j,d), 1e-8);
              EXPECT_NEAR(grad(i,j,d), r(i,j,d), 1e-8);
            }
        }
    }
}


TEST(TestSTensorFunction, Gradient2D)
{
  const char* g   = "sin(x)*sin(y) | exp(x)*exp(2*y) | x*x*y*y";
  const char* g_x = "cos(x)*sin(y) | exp(x)*exp(2*y) | 2*x*y*y";
  const char* g_y = "sin(x)*cos(y) | 2*exp(x)*exp(2*y) | 2*x*x*y";

  STensorFuncExpr f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{cos(x)*sin(y), 2*x*y*y, 2*x*y*y, exp(x)*exp(2*y),
                        sin(x)*cos(y), 2*x*x*y, 2*x*x*y, 2*exp(x)*exp(2*y)}.data());

      const utl::matrix3d<Real> grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const SymmTensor dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j) {
            EXPECT_DOUBLE_EQ(dx(i,j), r(i,j,d));
            EXPECT_DOUBLE_EQ(grad(i,j,d), r(i,j,d));
          }
      }
    }
}


TEST(TestSTensorFunction, Gradient2DFD)
{
  const char* g   = "sin(x)*sin(y) | exp(x)*exp(2*y) | x*x*y*y";

  const double eps = 1e-6;
  STensorFuncExpr f(g,"",eps);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7}) {
      const Vec3 X(x,y);
      const Vec3 Xp(x + 0.5*eps, y + 0.5*eps);
      const Vec3 Xm(x - 0.5*eps, y - 0.5*eps);
      utl::matrix3d<Real> r(2,2,2);
      r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y),
                        (Xp.x*Xp.x - Xm.x*Xm.x)*y*y,
                        (Xp.x*Xp.x - Xm.x*Xm.x)*y*y,
                        (exp(Xp.x) - exp(Xm.x))*exp(2*y),
                        sin(x)*(sin(Xp.y) - sin(Xm.y)),
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y),
                        x*x*(Xp.y*Xp.y - Xm.y*Xm.y),
                        exp(x)*(exp(2*Xp.y) - exp(2*Xm.y))}.data());
      r.multiply(1.0 / eps);

      const utl::matrix3d<Real> grad = f.gradient(X);
      for (size_t d = 1; d <= 2; ++d) {
        const SymmTensor dx = f.deriv(X,d);
        for (size_t i = 1; i <= 2; ++i)
          for (size_t j = 1; j <= 2; ++j) {
            EXPECT_NEAR(dx(i,j), r(i,j,d), 1e-8);
            EXPECT_NEAR(grad(i,j,d), r(i,j,d), 1e-8);
          }
      }
    }
}


TEST(TestSTensorFunction, Gradient3D)
{
  const char* g   = "sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*z*z |"
                    "x*y*z | x*x*y*z | z";
  const char* g_x = "cos(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | 2*x*y*y*z*z |"
                    "y*z | 2*x*y*z | 0";
  const char* g_y = "sin(x)*cos(y)*sin(z) | exp(x)*2*exp(2*y)*exp(z) | x*x*2*y*z*z |"
                    "x*z | x*x*z | 0";
  const char* g_z = "sin(x)*sin(y)*cos(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*2*z |"
                    "x*y | x*x*y | 1";

  STensorFuncExpr f(g);
  f.addDerivative(g_x,"",1);
  f.addDerivative(g_y,"",2);
  f.addDerivative(g_z,"",3);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{cos(x)*sin(y)*sin(z),
                          y*z,
                          0.0,
                          y*z,
                          exp(x)*exp(2*y)*exp(z),
                          2*x*y*z,
                          0.0,
                          2*x*y*z,
                          2*x*y*y*z*z,

                          sin(x)*cos(y)*sin(z),
                          x*z,
                          0.0,
                          x*z,
                          exp(x)*2*exp(2*y)*exp(z),
                          x*x*z,
                          0.0,
                          x*x*z,
                          x*x*2*y*z*z,

                          sin(x)*sin(y)*cos(z),
                          x*y,
                          1.0,
                          x*y,
                          exp(x)*exp(2*y)*exp(z),
                          x*x*y,
                          1.0,
                          x*x*y,
                          x*x*y*y*2*z}.data());

        const utl::matrix3d<Real> grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const SymmTensor dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j) {
              EXPECT_DOUBLE_EQ(dx(i,j), r(i,j,d));
              EXPECT_DOUBLE_EQ(grad(i,j,d), r(i,j,d));
            }
        }
      }
}


TEST(TestSTensorFunction, Gradient3DFD)
{
  const char* g   = "sin(x)*sin(y)*sin(z) | exp(x)*exp(2*y)*exp(z) | x*x*y*y*z*z |"
                    "x*y*z | x*x*y*z | z";

  const double eps = 1e-6;
  STensorFuncExpr f(g,"",eps);

  EXPECT_TRUE(f.isConstant());

  for (double x : {0.1, 0.2, 0.3})
    for (double y : {0.5, 0.6, 0.7})
      for (double z : {0.8, 0.9, 1.0}) {
        const Vec3 X(x,y,z);
        const Vec3 Xp(x + 0.5*eps, y + 0.5*eps, z + 0.5*eps);
        const Vec3 Xm(x - 0.5*eps, y - 0.5*eps, z - 0.5*eps);
        utl::matrix3d<Real> r(3,3,3);
        r.fill(std::array{(sin(Xp.x) - sin(Xm.x))*sin(y)*sin(z),
                          (Xp.x - Xm.x)*y*z,
                          0.0,
                          (Xp.x - Xm.x)*y*z,
                          (exp(Xp.x) - exp(Xm.x))*exp(2*y)*exp(z),
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*z,
                          0.0,
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*z,
                          (Xp.x*Xp.x - Xm.x*Xm.x)*y*y*z*z,

                          sin(x)*(sin(Xp.y) - sin(Xm.y))*sin(z),
                          x*(Xp.y - Xm.y)*z,
                          0.0,
                          x*(Xp.y - Xm.y)*z,
                          exp(x)*(exp(2*Xp.y) - exp(2*Xm.y))*exp(z),
                          x*x*(Xp.y - Xm.y)*z,
                          0.0,
                          x*x*(Xp.y - Xm.y)*z,
                          x*x*(Xp.y*Xp.y - Xm.y*Xm.y)*z*z,

                          sin(x)*sin(y)*(sin(Xp.z) - sin(Xm.z)),
                          x*y*(Xp.z - Xm.z),
                          Xp.z - Xm.z,
                          x*y*(Xp.z - Xm.z),
                          exp(x)*exp(2*y)*(exp(Xp.z) -  exp(Xm.z)),
                          x*x*y*(Xp.z - Xm.z),
                          Xp.z - Xm.z,
                          x*x*y*(Xp.z - Xm.z),
                          x*x*y*y*(Xp.z*Xp.z - Xm.z*Xm.z)}.data());
        r.multiply(1.0 / eps);

        const utl::matrix3d<Real> grad = f.gradient(X);
        for (size_t d = 1; d <= 3; ++d) {
          const SymmTensor dx = f.deriv(X,d);
          for (size_t i = 1; i <= 3; ++i)
            for (size_t j = 1; j <= 3; ++j) {
              EXPECT_NEAR(dx(i,j), r(i,j,d), 1e-8);
              EXPECT_NEAR(grad(i,j,d), r(i,j,d), 1e-8);
            }
        }
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
