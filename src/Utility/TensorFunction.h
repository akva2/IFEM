// $Id$
//==============================================================================
//!
//! \file TensorFunction.h
//!
//! \date May 10 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Spatial tensor-valued functions.
//!
//==============================================================================

#ifndef _TENSOR_FUNCTION_H
#define _TENSOR_FUNCTION_H

#include "Function.h"
#include "matrixnd.h"
#include "Tensor.h"


/*!
  \brief Tensor-valued unary function of a spatial point.
*/

class TensorFunc : public utl::Function<Vec3,Tensor>, public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  explicit TensorFunc(size_t n = 0)
  {
    ncmp = n*n;
  }

public:
  //! \brief Empty destructor.
  virtual ~TensorFunc() {}

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const { return 3; }

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3& X) const
  {
    return this->evaluate(X);
  }

  //! \brief Returns a representative scalar equivalent of the function value.
  virtual Real getScalarValue(const Vec3& X) const
  {
    return this->evaluate(X).trace();
  }

  //! \brief Evaluates first derivatives of the function.
  utl::matrix3d<Real> gradient(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    utl::matrix3d<Real> result(nsd, nsd, nsd);
    result.fill(this->evalGradient(X).data());

    return result;
  }

  //! \brief Evaluates time derivatives of the function.
  Tensor timeDerivative(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    Tensor result(nsd);
    result = this->evalTimeDerivative(X);
    return result;
  }

  //! \brief Evaluates second derivatives of the function.
  utl::matrix4d<Real> hessian(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    utl::matrix4d<Real> result(nsd, nsd, nsd, nsd);
    result.fill(this->evalHessian(X).data());

    return result;
  }

protected:
  //! \brief Returns the gradient of the function as a 1D array.
  virtual std::vector<Real> evalGradient(const Vec3& X) const
  {
    return std::vector<Real>(ncmp);
  }

  //! \brief Returns the second derivatives of the function as a 1D array.
  virtual std::vector<Real> evalHessian(const Vec3& X) const
  {
    return {};
  }

  //! \brief Returns the time derivatives of the function as a 1D array.
  virtual std::vector<Real> evalTimeDerivative(const Vec3& X) const
  {
    return std::vector<Real>(ncmp);
  }
};


/*!
  \brief Symmetric tensor-valued unary function of a spatial point.
*/

class STensorFunc : public utl::Function<Vec3,SymmTensor>,
                    public FunctionBase
{
  size_t index (size_t nsd, size_t i, size_t j) const
  {
    if (i == j)
      return i-1; // diagonal term
    else if (nsd == 2)
      return ncmp-1; // off-diagonal term (2D)

    if (i == j+1 || i+2 == j)
      std::swap(i,j);
    return i+2; // upper triangular term (3D)
  }

protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  STensorFunc(size_t n = 0, bool with33 = false)
  {
    ncmp = SymmTensor(n,with33).size();
  }

public:
  //! \brief Empty destructor.
  virtual ~STensorFunc() {}

  //! \brief Returns the function type flag.
  virtual unsigned char getType() const { return 3; }

  //! \brief Returns the function value as an array.
  virtual std::vector<Real> getValue(const Vec3& X) const
  {
    return this->evaluate(X);
  }

  //! \brief Returns a representative scalar equivalent of the function value.
  virtual Real getScalarValue(const Vec3& X) const
  {
    return this->evaluate(X).trace();
  }

  //! \brief Evaluates first derivatives of the function.
  utl::matrix3d<Real> gradient(const Vec3& X) const
  {
    const size_t nsd = ncmp > 5 ? 3 : (ncmp > 2 ? 2 : 1);
    utl::matrix3d<Real> result(nsd,nsd,nsd);
    const std::vector<Real> temp = this->evalGradient(X);

    for (size_t d = 1; d <= nsd; ++d)
      for (size_t i = 1; i <= nsd; ++i)
        for (size_t j = 1; j <= nsd; ++j)
          result(i,j,d) = temp[index(nsd,i,j) + (d-1)*ncmp];

    return result;
  }

  //! \brief Evaluates second derivatives of the function.
  utl::matrix4d<Real> hessian(const Vec3& X) const
  {
    const size_t nsd = ncmp > 5 ? 3 : (ncmp > 2 ? 2 : 1);
    utl::matrix4d<Real> result(nsd, nsd, nsd, nsd);
    const std::vector<Real> temp = this->evalHessian(X);

    // de-voigt the blocks
    auto it_out = result.begin();
    for (size_t j = 1; j <= nsd; ++j)
      for (size_t i = 1; i <= nsd; ++i, it_out += nsd*nsd) {
        const size_t idx = index(nsd,i,j);
        std::copy(temp.begin() + idx*nsd*nsd,
                  temp.begin() + (idx+1)*nsd*nsd, it_out);
      }

    return result;
  }

  //! \brief Evaluates time derivatives of the function.
  SymmTensor timeDerivative(const Vec3& X) const
  {
    const size_t nsd = ncmp > 5 ? 3 : (ncmp > 2 ? 2 : 1);
    SymmTensor result(nsd, ncmp == 4);
    result = this->evalTimeDerivative(X);
    return result;
  }

protected:
  //! \brief Returns the gradient of the function as a 1D array.
  virtual std::vector<Real> evalGradient(const Vec3& X) const
  {
    return std::vector<Real>(ncmp*ncmp);
  }

  //! \brief Returns the second derivatives of the function as a 1D array.
  virtual std::vector<Real> evalHessian(const Vec3& X) const
  {
    return {};
  }

  //! \brief Returns the time derivatives of the function as a 1D array.
  virtual std::vector<Real> evalTimeDerivative(const Vec3& X) const
  {
    return std::vector<Real>(ncmp);
  }
};

#endif
