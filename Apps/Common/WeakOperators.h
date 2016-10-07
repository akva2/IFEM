//==============================================================================
//!
//! \file EqualOrderOperators.h
//!
//! \date Jul 22 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various equal-ordered discrete operators.
//!
//==============================================================================

#ifndef WEAKOPERATORS_H_
#define WEAKOPERATORS_H_

class Vec3;

#include "FiniteElement.h"
#include "MatVec.h"

/*! \brief Common discrete operators using equal-ordered discretizations.
 */

class EqualOrderOperators
{
public:
  //! \brief Common weak operators using equal-ordered discretizations.
  class Weak {
  public:
    //! \brief Compute an advection term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] AC Advecting field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Advection(Matrix& EM, const FiniteElement& fe,
                          const Vec3& AC, double scale=1.0, int basis=1);

    //! \brief Wrapper checking whether blocks are in use or not.
    static void Advection(std::vector<Matrix>& EM, const FiniteElement& fe,
                          const Vec3& AC, double scale=1.0, int basis=1)
    {
      Advection(EM[EM.size() > 1 ? 1 : 0], fe, AC, scale, basis);
    }

    //! \brief Compute a (nonlinear) convection term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U  Advecting field
    //! \param[in] conservative True to use the conservative formulation
    //! \param[in] basis Basis to use
    template<class T>
    static void Convection(Matrix& EM, const FiniteElement& fe,
                           const Vec3& U, const T& dUdX, double scale,
                           bool conservative=false, int basis=1)
    {
      size_t cmp = EM.rows() / fe.basis(basis).size();
      if (conservative) {
        Advection(EM, fe, U, -scale, basis);
        for (size_t i = 1;i <= fe.basis(basis).size();i++)
          for (size_t j = 1;j <= fe.basis(basis).size();j++) {
            for (size_t k = 1;k <= cmp;k++) {
              for (size_t l = 1;l <= cmp;l++)
                EM((j-1)*cmp+l,(i-1)*cmp+k) -= scale*U[l-1]*fe.basis(basis)(i)*fe.grad(basis)(j,k)*fe.detJxW;
            }
          }
      }
      else {
        Advection(EM, fe, U, scale, basis);
        for (size_t i = 1;i <= fe.basis(basis).size();i++)
          for (size_t j = 1;j <= fe.basis(basis).size();j++) {
            for (size_t k = 1;k <= cmp;k++) {
              for (size_t l = 1;l <= cmp;l++)
                EM((j-1)*cmp+l,(i-1)*cmp+k) += scale*dUdX(l,k)*fe.basis(basis)(j)*fe.basis(basis)(i)*fe.detJxW;
            }
          }
      }
    }

    //! \brief Compute a divergence term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis for field
    //! \param[in] tbasis Test function basis
    static void Divergence(Matrix& EM, const FiniteElement& fe,
                           double scale=1.0, int basis=1, int tbasis=1);

    //! \brief Compute a divergence term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] D Divergence of field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Test function basis
    static void Divergence(Vector& EV, const FiniteElement& fe,
                           const Vec3& D, double scale=1.0, int basis=1);

    //! \brief Compute a gradient term for a (potentially mixed) vector/scalar field.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis for field
    //! \param[in] tbasis Test function basis
    static void Gradient(Matrix& EM, const FiniteElement& fe,
                         double scale=1.0, int basis=1, int tbasis=1);

    //! \brief Compute a gradient term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] tbasis Test function basis
    static void Gradient(Vector& EV, const FiniteElement& fe,
                         double scale=1.0, int basis=1);

    //! \brief Wrapper checking whether blocks are in use or not.
    static void Gradient(Vectors& EV, const FiniteElement& fe,
                         double scale=1.0, int basis=1)
    {
      Gradient(EV[EV.size() > 1 ? 1 : 0], fe, scale, basis);
    }

    //! \brief Compute a laplacian.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    //! \param[in] basis Basis to use
    static void Laplacian(Matrix& EM, const FiniteElement& fe,
                          double scale=1.0, bool stress=false, int basis=1);

    //! \brief Wrapper checking whether blocks are in use or not.
    static void Laplacian(std::vector<Matrix>& EM, const FiniteElement& fe,
                          double scale=1.0, bool stress=false, int basis=1)
    {
      Laplacian(EM[EM.size() > 1 ? 1 : 0], fe, scale, stress, basis);
    }

    //! \brief Compute a heteregenous coefficient laplacian.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[out] K The coefficient matrix
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    static void LaplacianCoeff(Matrix& EM, const Matrix& K, const FiniteElement& fe,
                               double scale=1.0, int basis=1);

    //! \brief Compute a mass term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Mass(Matrix& EM, const FiniteElement& fe,
                     double scale=1.0, int basis=1);

    //! \brief Wrapper checking whether blocks are in use or not.
    static void Mass(std::vector<Matrix>& EM, const FiniteElement& fe,
                     double scale=1.0, int basis=1)
    {
      Mass(EM[EM.size() > 1 ? 1 : 0], fe, scale, basis);
    }

    //! \brief Compute a source term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    //! \param[in] cmp Component to add (0 for all)
    static void Source(Vector& EV, const FiniteElement& fe,
                       double scale=1.0, int cmp=1, int basis=1);

    //! \brief Wrapper checking whether blocks are in use or not.
    static void Source(Vectors& EV, const FiniteElement& fe,
                       double scale=1.0, int cmp=1, int basis=1)
    {
      Source(EV[EV.size() > 1 ? 1 : 0], fe, scale, cmp, basis);
    }

    //! \brief Compute a vector-source term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Vector with contributions
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Source(Vector& EV, const FiniteElement& fe,
                       const Vec3& f, double scale=1.0, int basis=1);

    //! \brief Wrapper checking whether blocks are in use or not.
    static void Source(Vectors& EV, const FiniteElement& fe,
                       const Vec3& f, double scale=1.0, int basis=1)
    {
      Source(EV[EV.size() > 1 ? 1 : 0], fe, f, scale, basis);
    }
  };

  //! \brief Common weak residual operators using equal-ordered discretizations.
  class Residual {
  public:
    //! \brief Compute a convection term in a residual vector.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U  Advected field
    //! \param[in] dUdX Advected field gradient
    //! \param[in] UC Advecting field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] conservative True to use the conservative formulation
    //! \param[in] basis Basis to use
    template<class T>
    static void Convection(Vector& EV, const FiniteElement& fe,
                           const Vec3& U, const T& dUdX, const Vec3& UC,
                           double scale, bool conservative=false, int basis=1)
    {
      size_t cmp = EV.size() / fe.basis(basis).size();
      if (conservative) {
        for (size_t i = 1;i <= fe.basis(basis).size();i++)
          for (size_t k = 1;k <= cmp;k++)
            for (size_t l = 1;l <= cmp;l++)
              // Convection
              EV((i-1)*cmp + k) += scale*U[k-1]*UC[l-1]*fe.grad(basis)(i,l)*fe.detJxW;
      }
      else {
        for (size_t k = 1;k <= cmp;k++) {
          // Convection
          double conv = 0.0;
          for (size_t l = 1;l <= cmp;l++)
            conv += UC[l-1]*dUdX(k,l);
          conv *= scale*fe.detJxW;

          for (size_t i = 1;i <= fe.basis(basis).size();i++)
            EV((i-1)*cmp + k) -= conv*fe.basis(basis)(i);
        }
      }
    }

    //! \brief Compute a divergence term in a residual vector.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Gradient of field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    template<class T>
    static void Divergence(Vector& EV, const FiniteElement& fe,
                           const T& dUdX, double scale=1.0,
                           size_t basis=1)
    {
      for (size_t i = 1; i <= fe.basis(basis).size(); ++i) {
        double div=0.0;
        for (size_t k = 1; k <= fe.grad(basis).cols(); ++k)
          div += dUdX(k,k);
        EV(i) += scale*div*fe.basis(basis)(i)*fe.detJxW;
      }
    }

    //! \brief Compute a laplacian term in a residual vector.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Current solution gradient
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    //! \param[in] basis Basis to use
    template<class T>
    static void Laplacian(Vector& EV, const FiniteElement& fe,
                          const T& dUdX, double scale=1.0,
                          bool stress=false, int basis=1)
    {
      size_t cmp = EV.size() / fe.basis(basis).size();
      for (size_t i = 1;i <= fe.basis(basis).size();i++) {
        for (size_t k = 1;k <= cmp;k++) {
          double diff = 0.0;
          for (size_t l = 1;l <= cmp;l++)
            diff += dUdX(k,l)*fe.grad(basis)(i,l);
          diff *= scale*fe.detJxW;

          // Add negative residual to rhs of momentum equation
          EV((i-1)*cmp + k) -= diff;
        }
      }

      // Use stress formulation
      if (stress) {
        for (size_t i = 1;i <= fe.basis(basis).size();i++)
          for (size_t k = 1;k <= cmp;k++)
            for (size_t l = 1;l <= cmp;l++)
              // Diffusion
              EV((i-1)*cmp + k) -= scale*dUdX(l,k)*fe.grad(basis)(i,l)*fe.detJxW;
      }
    }
  };
};

#endif
