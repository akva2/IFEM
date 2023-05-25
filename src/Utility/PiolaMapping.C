// $Id$
//==============================================================================
//!
//! \file PiolaMapping.C
//!
//! \date Apr 13 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utilities for Piola mapping transformations.
//!
//==============================================================================

#include "PiolaMapping.h"

#include "FiniteElement.h"

namespace utl
{

bool piolaMapping (FiniteElement& fe,
                   const matrix<Real>& J,
                   const matrix<Real>& Ji,
                   const std::vector<matrix<Real>>& dNdu,
                   const matrix3d<Real>& H)
{
  // Only implemented in 2D for now
  if (fe.grad(1).cols() != 2 || fe.getNoBasis() != 3)
    return false;

  // Apply piola map to basis functions
  double Jt1 = H(1,1,1)*J(2,2) + J(1,1)*H(2,1,2) - H(1,1,2)*J(2,1) - J(1,2)*H(2,1,1);
  double Jt2 = H(1,1,2)*J(2,2) + J(1,1)*H(2,2,2) - H(1,2,2)*J(2,1) - J(1,2)*H(2,1,2);

  double detJ_x = Ji(1,1)*Jt1 + Ji(2,1)*Jt2;
  double detJ_y = Ji(1,2)*Jt1 + Ji(2,2)*Jt2;

  Matrix Jx(2,2);
  Jx(1,1) = Ji(1,1)*H(1,1,1) + Ji(2,1)*H(1,1,2);
  Jx(1,2) = Ji(1,1)*H(1,1,2) + Ji(2,1)*H(1,2,2);
  Jx(2,1) = Ji(1,1)*H(2,1,1) + Ji(2,1)*H(2,1,2);
  Jx(2,2) = Ji(1,1)*H(2,2,1) + Ji(2,1)*H(2,2,2);
  Jx *= 1.0 / fe.detJxW;

  Matrix Jy(2,2);
  Jy(1,1) = Ji(1,2)*H(1,1,1) + Ji(2,2)*H(1,1,2);
  Jy(1,2) = Ji(1,2)*H(1,1,2) + Ji(2,2)*H(1,2,2);
  Jy(2,1) = Ji(1,2)*H(2,1,1) + Ji(2,2)*H(2,1,2);
  Jy(2,2) = Ji(1,2)*H(2,1,2) + Ji(2,2)*H(2,2,2);
  Jy *= 1.0 / fe.detJxW;

  Matrix Px = Jx;
  Px.add(J, -detJ_x/(fe.detJxW*fe.detJxW));
  Matrix Py = Jy;
  Py.add(J, -detJ_y/(fe.detJxW*fe.detJxW));

  size_t NP = fe.basis(1).size() + fe.basis(2).size();
  fe.dPndX.resize(4, NP);

  for (size_t i = 1; i <= fe.basis(1).size(); ++i) {
    Vector bf(2);
    bf(1) = fe.basis(1)(i);

    Vector du(2);
    du(1) = dNdu[0](i,1)*Ji(1,1) + dNdu[0](i,2)*Ji(2,1);
    du *= 1.0 / fe.detJxW;
    Vector dudx = Px*bf;
    J.multiply(du, dudx, false, 1);
    fe.dPndX(1,i) = dudx(1);
    fe.dPndX(2,i) = dudx(2);

    Vector dv(2);
    dv(1) = dNdu[0](i,1)*Ji(1,2) + dNdu[0](i,2)*Ji(2,2);
    dv *= 1.0 / fe.detJxW;
    Vector dudy = Py*bf;
    J.multiply(dv, dudy, false, 1);
    fe.dPndX(3,i) = dudy(1);
    fe.dPndX(4,i) = dudy(2);
  }

  const size_t N1 = fe.basis(1).size();
  for (size_t j = 1; j <= fe.basis(2).size(); ++j) {
    Vector bf(2);
    bf(2) = fe.basis(2)(j);

    Vector du(2);
    du(2) = dNdu[1](j,1)*Ji(1,1) + dNdu[1](j,2)*Ji(2,1);
    du *= 1.0 / fe.detJxW;
    Vector dudx = Px*bf;
    J.multiply(du, dudx, false, 1);
    fe.dPndX(1,j + N1) = dudx(1);
    fe.dPndX(2,j + N1) = dudx(2);

    Vector dv(2);
    dv(2) = dNdu[1](j,1)*Ji(1,2) + dNdu[1](j,2)*Ji(2,2);
    dv *= 1.0 / fe.detJxW;
    Vector dudy = Py*bf;
    J.multiply(dv, dudy, false, 1);
    fe.dPndX(3,j + N1) = dudy(1);
    fe.dPndX(4,j + N1) = dudy(2);
  }

  Matrix Jt(2,1);
  Jt(1,1) = Jt1;
  Jt(2,1) = Jt2;
  fe.grad(3).multiply(Jt, fe.basis(3), 1.0 / fe.detJxW, -1.0);

  // Vector transform for velocity basis
  utl::piolaBasis(fe, J);

  return true;
}


void piolaBasis (FiniteElement& fe,
                 const utl::matrix<Real>& J)
{
  // Vector transform for velocity basis
  size_t NP = fe.basis(1).size() + fe.basis(2).size();
  Matrix N(2, NP);
  for (size_t i = 1; i <= fe.basis(1).size(); ++i)
    N(1,i) = fe.basis(1)(i);
  for (size_t i = 1; i <= fe.basis(2).size(); ++i)
    N(2,fe.basis(1).size() + i) = fe.basis(2)(i);
  fe.PN.multiply(J,N,false,false,false,1.0 / fe.detJxW);

  fe.basis(3) *= 1.0 / fe.detJxW;
}

}
