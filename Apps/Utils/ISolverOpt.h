//==============================================================================
//!
//! \file ISolverOpt.h
//!
//! \date Sep 19 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Abstract simulator interface for adjoint-based optimal control.
//!
//==============================================================================

#ifndef _I_SOLVER_OPT_H_
#define _I_SOLVER_OPT_H_

#include "MatVec.h"


/*!
  \brief Abstract interface for PDE-constrained optimal control simulators.

  \details A simulator implementing this interface can be driven by the
  SIMSolverOpt template. The interface describes the reduced-space formulation
  of the optimal control problem

  \f[ \min_{f} J(f,u(f)) \quad\mbox{subject to}\quad e(u,f) = 0 \f]

  where \b f is the control, \b u is the state and \f$e(u,f)=0\f$ is the
  discretized state (forward) equation. All vectors exchanged with the driver
  are in DOF-order, and the control and the gradient live in the same space.

  The driver never forms the reduced Hessian or the mass matrix itself.
  It only needs the four operators solveForward(), solveAdjoint(),
  computeGradient() and applyHessian(), plus the metric applyMass() in which
  norms and inner products of controls and gradients are measured.

  \note The interface is deliberately not enforced through inheritance, since
  the simulator templates in IFEM are combined statically. This class therefore
  only serves as documentation of the methods SIMSolverOpt expects to find.
*/

class ISolverOpt
{
public:
  //! \brief Initializes the optimal control formulation.
  //! \details This is invoked once before the optimization loop is entered.
  //! It is the place to assemble and factorize the quantities that stay
  //! constant throughout the optimization, e.g. the state coefficient matrix,
  //! the control space mass matrix and the data of the objective function.
  virtual bool initOptimalControl() = 0;

  //! \brief Returns the initial control field.
  //! \param[out] F Control vector in DOF-order
  virtual bool getControl(Vector& F) const = 0;

  //! \brief Solves the forward (state) problem for a given control.
  //! \param[in] F Control vector
  //! \param[out] U Resulting state vector
  virtual bool solveForward(const Vector& F, Vector& U) = 0;

  //! \brief Evaluates the objective function.
  //! \param[in] F Control vector
  //! \param[in] U State vector associated with \a F
  virtual double computeObjective(const Vector& F, const Vector& U) const = 0;

  //! \brief Solves the adjoint problem.
  //! \param[in] U Current state vector
  //! \param[in] F Current control vector
  //! \param[out] P Resulting adjoint vector
  virtual bool solveAdjoint(const Vector& U, const Vector& F, Vector& P) = 0;

  //! \brief Evaluates the reduced gradient of the objective function.
  //! \param[in] F Current control vector
  //! \param[in] P Adjoint vector associated with \a F
  //! \param[out] G Riesz representation of the gradient, in the metric
  //! defined by applyMass()
  virtual bool computeGradient(const Vector& F, const Vector& P,
                               Vector& G) const = 0;

  //! \brief Applies the reduced Hessian to a direction.
  //! \param[in] D Direction in control space
  //! \param[out] HD The Hessian-vector product, in the same metric as the
  //! gradient returned by computeGradient()
  //!
  //! \details Only needed by the Newton-CG method. Simulators that do not
  //! provide second order information should return \e false, which makes
  //! the driver truncate the inner iteration, falling back to a steepest
  //! descent step if this happens already on the first one.
  virtual bool applyHessian(const Vector& D, Vector& HD) = 0;

  //! \brief Applies the control space metric (mass matrix) to a vector.
  //! \param[in] X Vector in control space
  //! \param[out] MX The product \b M \b X
  //!
  //! \details Returning \e false makes the driver fall back to the Euclidean
  //! inner product of the DOF values, which is mesh dependent.
  virtual bool applyMass(const Vector& X, Vector& MX) const = 0;

  //! \brief Saves the converged optimal solution for postprocessing.
  //! \param[in] F Optimal control vector
  //! \param[in] U State vector associated with \a F
  //! \param[in] P Adjoint vector associated with \a F
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  virtual bool saveOptimalResults(const Vector& F, const Vector& U,
                                  const Vector& P, char* fileName,
                                  int& geoBlk, int& nBlock) = 0;

  //! \brief Prints errors with respect to an analytical solution, if any.
  //! \param[in] F Optimal control vector
  //! \param[in] U State vector associated with \a F
  virtual void printErrors(const Vector& F, const Vector& U,
                           const Vector& P) = 0;
};

#endif
