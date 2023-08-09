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

void piolaMapping (FiniteElement& fe,
                   const double detJ,
                   const matrix<Real>& J,
                   const matrix<Real>& Ji,
                   const std::vector<matrix<Real>>& dNdu,
                   const matrix3d<Real>& H)
{
  piolaBasis(fe, J);
  piolaGradient(fe, detJ, J, Ji, dNdu, H);
}


std::vector<matrix<Real>> jacobianGradient (const matrix<Real>& dudX,
                                            const matrix3d<Real>& d2Xdu2)
{
  const size_t dim = dudX.rows();

  Matrices dJdU(dim);
  for (size_t k = 1; k <= dim; ++k) {
    Matrix& dJ = dJdU[k-1];
    dJ.resize(dim, dim);
    for (size_t i = 1; i <= dim; ++i)
      for (size_t j = 1; j <= dim; ++j)
        dJ(i,j) = d2Xdu2(i,k,j);
  }

  Matrices dJdX(dim);
  for (size_t j = 1; j <= dim; ++j) {
    Matrix& dJ = dJdX[j-1];
    dJ.resize(dim, dim);
    for (size_t i = 0; i < dim; ++i)
      dJ.add(dJdU[i], dudX(i+1,j));
  }

  return dJdX;
}


vector<Real> determinantGradient (const matrix<Real>& J,
                                  const matrix<Real>& Ji,
                                  const matrix3d<Real>& H)
{
  const size_t dim = J.rows();
  vector<Real> ddet(dim), ddetu(dim);

  if (dim == 2) {
    ddetu(1) =   H(1,1,1)*J(2,2) + J(1,1)*H(2,2,1)
               - H(1,2,1)*J(2,1) - J(1,2)*H(2,1,1);

    ddetu(2) =   H(1,1,2)*J(2,2) + J(1,1)*H(2,2,2)
               - H(1,2,2)*J(2,1) - J(1,2)*H(2,1,2);
  } else {
    auto det = [&J](size_t i)
    {
        const size_t i1 = (i == 1 ? 2 : 1);
        const size_t i2 = (i == 3 ? 2 : 3);

        return J(2,i1) * J(3,i2) - J(2,i2) * J(3,i1);
    };

    auto detD = [&J,&H](size_t i, size_t d)
    {
        const size_t i1 = (i == 1 ? 2 : 1);
        const size_t i2 = (i == 3 ? 2 : 3);

        return   H(2,i1,d) * J(3,i2) + J(2,i1) * H(3,i2,d)
               - H(2,i2,d) * J(3,i1) - J(2,i2) * H(3,i1,d);
    };

    for (size_t i = 1; i <= 3; ++i)
      ddetu(i) =    H(1,1,i) * det(1) + J(1,1) * detD(1, i)
                 - (H(1,2,i) * det(2) + J(1,2) * detD(2, i))
                 +  H(1,3,i) * det(3) + J(1,3) * detD(3, i);
  }

  for (size_t i = 1; i <= dim; ++i)
    for (size_t j = 1; j <= dim; ++j)
      ddet(i) += ddetu(j) * Ji(j,i);

  return ddet;
}

void piolaBasis (FiniteElement& fe,
                 const utl::matrix<Real>& J)
{
  const size_t dim = J.rows();
  size_t NP = 0;
  for (size_t b = 1; b <= dim; ++b)
    NP += fe.basis(b).size();

  Matrix N(dim, NP);
  size_t k = 1;
  for (size_t b = 1; b <= dim; ++b)
    for (size_t i = 1; i <= fe.basis(b).size(); ++i, ++k)
      N(b,k) = fe.basis(b)(i);

  fe.piola.N.multiply(J,N,false,false,false,1.0 / fe.detJxW);
}


void piolaGradient (FiniteElement& fe,
                    const Real detJ,
                    const matrix<Real>& J,
                    const matrix<Real>& Ji,
                    const std::vector<matrix<Real>>& dNdu,
                    const matrix3d<Real>& H)
{
  const size_t dim = J.rows();
  size_t NP = 0;
  for (size_t b = 1; b <= dim; ++b)
    NP += fe.basis(b).size();
  fe.piola.dNdX.resize(dim*dim, NP);

  const Matrices dJdX = jacobianGradient(Ji, H);
  const auto det = determinantGradient(J, Ji, H);

  std::vector<Matrix> P(dim);
  for (size_t i = 0; i < dim; ++i) {
    P[i] = dJdX[i];
    P[i] *= 1.0 / detJ;
    P[i].add(J, -det[i] / (detJ * detJ));
  }

  size_t k = 0;
  for (size_t b = 1; b <= dim; k += fe.basis(b).size(), ++b)
    for (size_t i = 1; i <= fe.basis(b).size(); ++i) {
      Vector bf(dim);
      bf(b) = fe.basis(b)(i);
      Matrix dBf(dim,dim);
      for (size_t j = 1; j <= dim; ++j)
        dBf(b,j) = dNdu[b-1](i,j);
      for (size_t d = 1; d <= dim; ++d) {
        Vector du(dim);
        for (size_t j = 1; j <= dim; ++j)
          du(j) = Ji(j,d);
        Vector tmp(dim), res(dim);
        dBf.multiply(du, tmp);
        J.multiply(tmp, res, 1.0 / detJ);
        P[d-1].multiply(bf, res, false, 1);
        for (size_t j = 1; j <= dim; ++j)
          fe.piola.dNdX((d-1) * dim + j,k + i) = res(j);
      }
    }
}

}
