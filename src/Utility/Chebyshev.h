// $Id$
//==============================================================================
//!
//! \file Chebyshev.h
//!
//! \date Jul 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Evaluation of Chebyshev polynomials.
//!
//==============================================================================

#ifndef _CHEBYSHEV_H
#define _CHEBYSHEV_H

#include "Function.h"
#include "TensorFunction.h"

#include <array>
#include <memory>


/*!
  \brief Evaluation of Chebyshev polynomials.
*/

namespace Chebyshev
{
  //! \brief Evaluates a 1D Chebyshev polynomial of first kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  double evalPol1(int polnum, double xi);
  //! \brief Evaluates a 1D Chebyshev polynomial of second kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  double evalPol2(int polnum, double xi);

  //! \brief Evaluates the first derivative of a 1D Chebyshev polynomial of first kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  double evalDer1(int polnum, double xi);
  //! \brief Evaluates the first derivative of a 1D Chebyshev polynomial of second kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  double evalDer2(int polnum, double xi);

  //! \brief Evaluates the second derivative of a 1D Chebyshev polynomial of first kind.
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  double eval2Der1(int polnum, double xi);
}


/*!
  \brief A scalar-valued spatial function, chebyshev polynomials.
*/

class ChebyshevFunc : public RealFunc
{
  std::vector<Real> coefs; //!< Function coefficients
  std::array<int,3> n;     //!< Number of coefficients

public:
  //! \brief The constructor initializes the function parameters from a file.
  //! \param[in] file Name of file to read \ref coefs from
  ChebyshevFunc(const char* file);

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return n[0] == n[1] == n[2] == 0; }

  //! \brief Returns a const reference to the coefficients.
  const std::vector<Real>& getCoefs() const { return coefs; }
  //! \brief Returns the number of polynomials in each parameter direction.
  const std::array<int,3>& getSize() const { return n; }

protected:
  //! \brief Evaluates the function at point \a X.
  virtual Real evaluate(const Vec3& X) const;
};


/*!
  \brief A vector-valued spatial function, chebyshev polynomials.
*/

class ChebyshevVecFunc : public VecFunc
{
  std::array<std::unique_ptr<ChebyshevFunc>,3> f; //!< Functions
  bool secondDer; //!< True to take second derivatives

public:
  //! \brief The constructor initializes the function parameters from a file.
  //! \param[in] file Name of files to read coefs from
  //! \param[in] second True to take second derivatives
  ChebyshevVecFunc(const std::vector<const char*>& file, bool second = false);

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return f[0]->isZero(); }

protected:
  //! \brief Evaluates the function at point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief A tensor-valued spatial function, chebyshev polynomials.

  \details If 2 or 3 functions: Take the derivative of the interpolants.
           If 4 or 9 functions: Each component has their own interpolant.
*/

class ChebyshevTensorFunc : public TensorFunc
{
  std::array<std::unique_ptr<ChebyshevVecFunc>,3> f; //!< Array of vector components

public:
  //! \brief The constructor initializes the function parameters from files.
  //! \param[in] file Name of files to read from
  //! \param[in] second True to take second derivatives
  ChebyshevTensorFunc(const std::vector<const char*>& file, bool second);

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const { return !(f[0] || f[1] || f[2]); }

protected:
  //! \brief Evaluates the function at point \a X.
  virtual Tensor evaluate(const Vec3& X) const;
};

#endif
