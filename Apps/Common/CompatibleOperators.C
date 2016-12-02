//==============================================================================
//!
//! \file CompatibleOperators.C
//!
//! \date Oct 9 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, discrete div-compatible operators.
//!
//==============================================================================

#include "CompatibleOperators.h"
#include "EqualOrderOperators.h"
#include "FiniteElement.h"
#include "Vec3.h"


void CompatibleOperators::Weak::Advection(std::vector<Matrix>& EM,
                                          const FiniteElement& fe,
                                          const Vec3& AC, double scale)
{
  size_t nsd = fe.grad(1).cols();
  Vector c(AC.ptr(), nsd);
  for (size_t n = 1; n <= nsd; ++n)
    EM[n].outer_product(fe.basis(n), fe.grad(n)*c, true, scale*fe.detJxW);
}


void CompatibleOperators::Weak::Convection(std::vector<Matrix>& EM,
                                           const FiniteElement& fe,
                                           const Vec3& U,
                                           const Tensor& dUdX,
                                           double scale,
                                           bool conservative)
{
  // Convection
  static const double vidx[3][3] = {{1, 6, 7},
                                    {10, 2, 11},
                                    {14, 15, 3}};
  size_t nsd = fe.grad(1).cols();
  Vector u(U.ptr(), nsd);
  for (size_t m=1; m <= nsd; ++m)
    for (size_t n = 1; n <= nsd; ++n) {
      EM[vidx[m-1][n-1]].outer_product(fe.basis(m), fe.basis(n),
                                       true, dUdX(m,n)*scale*fe.detJxW);
      if (m == n)
        EM[vidx[m-1][n-1]].outer_product(fe.basis(m), fe.grad(n)*u,
                                         true, scale*fe.detJxW);
    }
}


void CompatibleOperators::Weak::Gradient(std::vector<Matrix>& EM,
                                        const FiniteElement& fe,
                                        double scale)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    for (size_t i=1; i <= fe.basis(n).size(); ++i)
      for (size_t j=1; j <= fe.basis(nsd+1).size(); ++j)
        EM[8+4*(n-1)](i,j) += -scale*fe.grad(n)(i,n)*fe.basis(nsd+1)(j)*fe.detJxW;
}


void CompatibleOperators::Weak::Gradient(Vectors& EV, const FiniteElement& fe,
                                         double scale)
{
  for (size_t n = 1; n <= fe.grad(n).cols(); ++n)
    for (size_t i = 1; i <= fe.basis(n).size(); ++i)
      EV[n](i) += scale*fe.grad(n)(i,n)*fe.detJxW;
}


void CompatibleOperators::Weak::Laplacian(std::vector<Matrix>& EM,
                                          const FiniteElement& fe,
                                          double scale, bool stress)
{
  size_t nsd = fe.grad(1).cols();
  for (size_t n = 1; n <= nsd; ++n)
    EqualOrderOperators::Weak::Laplacian(EM[n], fe, scale, false, n);

  for (size_t m = 1; m <= nsd && stress; m++)
    for (size_t n = m; n <= nsd; n++) {
      int idx = m == n ? m : (m == 1 ? 5+n-m : 10+n-m);
      EM[idx].multiply(fe.grad(m), fe.grad(n), false, true,
                       true, scale*fe.detJxW);
    }
}


void CompatibleOperators::Weak::Mass(std::vector<Matrix>& EM,
                                     const FiniteElement& fe, double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Mass(EM[k], fe, scale, k);
}

void CompatibleOperators::Weak::Source(Vectors& EV,
                                       const FiniteElement& fe,
                                       const Vec3& f, double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Source(EV[k], fe, scale*f[k-1], 1, k);
}


void CompatibleOperators::Weak::Source(Vectors& EV,
                                       const FiniteElement& fe,
                                       double scale)
{
  for (size_t k = 1; k <= fe.grad(1).cols(); ++k)
    EqualOrderOperators::Weak::Source(EV[k], fe, scale, 1, k);
}


void CompatibleOperators::Residual::Convection(Vectors& EV, const FiniteElement& fe,
                                               const Vec3& U, const Tensor& dUdX,
                                               const Vec3& UC, double scale,
                                               bool conservative)
{
  size_t nsd = fe.grad(1).cols();
  Vec3 udu = dUdX*U;
  for (size_t n=1; n <= nsd; ++n)
    EV[n].add(fe.basis(n), -scale*udu[n-1]*fe.detJxW);
}


void CompatibleOperators::Residual::Laplacian(Vectors& EV,
                                              const FiniteElement& fe,
                                              const Tensor& dUdX,
                                              double scale, bool stress)
{
  size_t nsd = fe.grad(1).cols();
  Matrix dUdXt(nsd, nsd);
  dUdXt = dUdX;
  for (size_t k = 1; k <= nsd; ++k)
    for (size_t i = 1; i <= fe.basis(k).size(); ++i) {
      double diff = 0.0;
      for (size_t m = 1; m <= nsd; ++m)
        diff += fe.grad(k)(i,m)*dUdX(k,m);
      EV[k](i) += scale*diff*fe.detJxW;
    }
}
