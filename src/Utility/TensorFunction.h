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
#include "Tensor.h"


/*!
  \brief Tensor-valued unary function of a spatial point.
*/

class TensorFunc : public utl::SpatialFunction<Tensor>, public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  explicit TensorFunc(size_t n = 0) : utl::SpatialFunction<Tensor>(Tensor(n))
  {
    ncmp = zero.size();
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

  utl::matrix3d<Real> gradient(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    utl::matrix3d<Real> result(nsd, nsd, nsd);
    result.fill(this->evalGradient(X).data());

    return result;
  }

  Tensor tgradient(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    Tensor result(nsd);
    result = this->evalTGradient(X);

    return result;
  }

  utl::matrix4d<Real> hessian(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    utl::matrix4d<Real> result(nsd, nsd, nsd, nsd);
    result.fill(this->evalHessian(X).data());

    return result;
  }

protected:
  virtual std::vector<Real> evalGradient(const Vec3& X) const
  {
    return std::vector<Real>(ncmp*this->dim()*this->dim());
  }

  virtual std::vector<Real> evalTGradient(const Vec3& X) const
  {
    return std::vector<Real>(this->dim()*this->dim());
  }

  virtual std::vector<Real> evalHessian(const Vec3& X) const
  {
    return std::vector<Real>(ncmp*this->dim()*this->dim()*this->dim());
  }
};


/*!
  \brief Symmetric tensor-valued unary function of a spatial point.
*/

class STensorFunc : public utl::SpatialFunction<SymmTensor>,
                    public FunctionBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  STensorFunc(size_t n = 0, bool with33 = false)
    : utl::SpatialFunction<SymmTensor>(SymmTensor(n,with33))
  {
    ncmp = zero.size();
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

  utl::matrix3d<Real> gradient(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    utl::matrix3d<Real> result(nsd, nsd, nsd);
    result.fill(this->evalGradient(X).data());

    return result;
  }

  SymmTensor tgradient(const Vec3& X) const
  {
    return SymmTensor(this->evalTGradient(X));
  }

  utl::matrix4d<Real> hessian(const Vec3& X) const
  {
    const size_t nsd = sqrt(this->dim());
    utl::matrix4d<Real> result(nsd, nsd, nsd, nsd);
    result.fill(this->evalHessian(X).data());

    return result;
  }

protected:
  virtual std::vector<Real> evalGradient(const Vec3& X) const
  {
    return std::vector<Real>(ncmp*this->dim()*this->dim());
  }

  virtual std::vector<Real> evalTGradient(const Vec3& X) const
  {
    return std::vector<Real>(this->dim()*this->dim());
  }

  virtual std::vector<Real> evalHessian(const Vec3& X) const
  {
    return std::vector<Real>(ncmp*this->dim()*this->dim()*this->dim());
  }
};

#endif
