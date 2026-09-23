//==============================================================================
//!
//! \file OptControlOps.h
//!
//! \date Sep 19 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Reduced-space operators for linear-quadratic optimal control.
//!
//==============================================================================

#ifndef _OPT_CONTROL_OPS_H_
#define _OPT_CONTROL_OPS_H_

#include "MatVec.h"

#include <memory>
#include <vector>

class SIMbase;
class SparseMatrix;


/*!
  \brief Reduced-space operators of a linear-quadratic optimal control problem.

  \details This class holds the algebra that is common to adjoint-based optimal
  control of any linear state equation with a tracking-type objective

  \f[ J(f) = \frac12\|u(f)-d\|^2 + \frac\alpha2\|f\|^2 \f]

  where the control \b f enters the state equation as a distributed source. It
  owns the control space mass matrix \b M, the load vector of the desired state
  \f${\bf b}_d\f$, and the constant part \f${\bf b}_0\f$ of the state load
  vector, and provides the four operators the SIMSolverOpt driver needs.

  The state coefficient matrix \b K is the one held by the simulator, assembled
  and factorized once. Writing \b R for the restriction from nodal to equation
  ordering, the operators are

  \f[ {\bf K}{\bf u} = {\bf R}({\bf M}{\bf f}+{\bf b}_0) \qquad
      {\bf K}{\bf p} = {\bf R}({\bf M}{\bf u}-{\bf b}_d) \qquad
      \nabla J = \alpha{\bf f}+{\bf p} \f]

  and the reduced Hessian is applied by repeating the two solves with
  \f${\bf b}_0={\bf b}_d={\bf 0}\f$.

  The adjoint equation uses the transpose of \b K. For a symmetric state
  operator that is \b K itself. It is also \b K for a saddle-point operator
  assembled with a skew-symmetric coupling block, as IFEM does for Stokes:
  there \f${\bf K}^T={\bf S}{\bf K}{\bf S}\f$ with \b S negating the
  constraint degrees of freedom, and since the adjoint right-hand-side has
  entries only in the equations the control acts on, the solve reduces to an
  ordinary one whose constraint part comes out with the opposite sign.

  The control need not live on all degrees of freedom of the model. For a mixed
  formulation it typically lives on the velocity degrees of freedom only, which
  is expressed by restricting it to the degrees of freedom of one nodal type.
*/

class OptControlOps
{
public:
  //! \brief The constructor initializes the reference to the simulator.
  //! \param sim The simulator holding the state equation
  explicit OptControlOps(SIMbase& sim);
  //! \brief The destructor deletes the control space mass matrix.
  ~OptControlOps();

  //! \brief Toggles use of a transposed solve for the adjoint equation.
  //! \details Needed when the state operator is not symmetric. It is left off
  //! by default, as a symmetric operator is its own transpose and the
  //! untransposed solve is then both cheaper and better tested.
  bool setAdjointTranspose(bool yes);

  //! \brief Defines the Tikhonov regularization parameter.
  void setRegularization(double a) { alpha = a; }
  //! \brief Returns the Tikhonov regularization parameter.
  double getRegularization() const { return alpha; }

  //! \brief Allocates the equation system of the state problem.
  bool initSystem();

  //! \brief Captures the constant part of the state load vector.
  //! \details To be invoked after the state equation system has been assembled.
  //! The captured vector carries any body load of the model, as well as the
  //! contributions from inhomogeneous Dirichlet conditions.
  bool captureFixedLoad();

  //! \brief Allocates and returns the control space mass matrix.
  SparseMatrix* initMassMatrix();

  //! \brief Allocates the coupling blocks of the control-to-load operator.
  //! \param[out] C The coupling block, to be filled by the simulator
  //! \param[out] Ct Where the simulator should assemble its transpose
  //!
  //! \details The control enters the state equation through the load operator
  //! \b L. When that is a multiple of the mass matrix, which is the case for
  //! an unstabilized formulation, nothing needs to be allocated here and the
  //! scaling alone describes it. A consistently stabilized formulation, on the
  //! other hand, also lets the control enter the equations it stabilizes with,
  //! and the resulting coupling block is what this allocates:
  //! \f${\bf L}=c_u{\bf M}+c_c{\bf C}\f$.
  //! \param[out] Mr Where the simulator should assemble a second copy of the
  //! mass matrix, which this class factorizes for the Riesz map
  void initLoadOperator(SparseMatrix*& C, SparseMatrix*& Ct,
                        SparseMatrix*& Mr);

  //! \brief Defines the scaling of the two parts of the load operator.
  //! \param[in] cu Scaling of the mass matrix part
  //! \param[in] cc Scaling of the coupling block
  void setLoadScaling(double cu, double cc = 1.0)
  {
    massScale = cu;
    coupScale = cc;
  }
  //! \brief Returns the control space mass matrix.
  SparseMatrix* getMassMatrix() const { return myMass.get(); }
  //! \brief Returns the load vector of the desired state.
  Vector& getTargetLoad() { return myTarget; }
  //! \brief Defines the squared \f$L_2\f$-norm of the desired state.
  void setTargetNorm(double d2) { normTarget2 = d2; }
  //! \brief Returns the \f$L_2\f$-norm of the desired state.
  double getTargetNorm() const;

  //! \brief Restricts the control to a subset of the degrees of freedom.
  //! \param[in] nodeTypes Nodal DOF classifications the control lives on,
  //! empty for all nodes in the model
  //! \param[in] nComp If nonzero, only the first \a nComp degrees of freedom
  //! of each selected node belong to the control
  //!
  //! \details The two criteria cover the layouts the mixed and non-mixed
  //! formulations give rise to. A field discretized on one basis of a mixed
  //! model is selected by the nodal type of that basis, whereas a field that
  //! shares its nodes with other fields, as in an equal-order formulation, is
  //! selected by the number of leading components.
  bool setControlDofs(const std::vector<char>& nodeTypes = {},
                      unsigned char nComp = 0);
  //! \brief Returns the number of degrees of freedom the control lives on.
  size_t getNoControls() const;
  //! \brief Zeroes the entries of a vector outside the control space.
  void maskControl(Vector& X) const;

  //! \brief Applies the control space mass matrix to a vector.
  bool applyMass(const Vector& X, Vector& MX) const;

  //! \brief Applies the control-to-load operator to a control field.
  //! \param[in] F The control field
  //! \param[out] LF The resulting load vector, in nodal ordering
  bool applyLoad(const Vector& F, Vector& LF) const;

  //! \brief Maps a field in state space back onto the control space.
  //! \param[in] Y The field to map, typically an adjoint solution
  //! \param[out] G The result, as a field in the control space
  //!
  //! \details This is \f${\bf M}^{-1}{\bf L}^T\f$, that is the transpose of
  //! the load operator followed by the Riesz map of the control space. The
  //! mass matrix part of it needs no solve, as the two cancel; only the
  //! coupling block of a stabilized formulation leaves an \f$L_2\f$
  //! projection to be performed.
  bool applyLoadTranspose(const Vector& Y, Vector& G) const;

  //! \brief Solves the state problem for a given control field.
  bool solveForward(const Vector& F, Vector& U);
  //! \brief Solves the state equation with homogeneous boundary conditions.
  //! \param[in] nodalRHS Right-hand-side vector in nodal ordering
  //! \param[out] sol The resulting solution, in nodal ordering
  //! \param[in] adjoint If \e true, solve the adjoint of the state equation
  bool solveHomogeneous(const Vector& nodalRHS, Vector& sol,
                        bool adjoint = false);

private:
  //! \brief Solves the already assembled equation system.
  //! \param[out] sol The resulting solution, in nodal ordering
  //! \param[in] scaleSD Scaling factor for the prescribed values
  //!
  //! \details The state coefficient matrix is factorized by the first solve
  //! and reused by all the subsequent ones. The expansion of the solution is
  //! done here rather than through SIMbase::solveEqSystem(), since the adjoint
  //! and Hessian solves need the homogeneous version of the constraints,
  //! which that method only applies to a second right-hand-side vector, and
  //! allocating one of those makes the assembly loop try to assemble the
  //! additional element vectors of a block integrand into it.
  bool solveSystem(Vector& sol, double scaleSD, bool transposed = false);

  //! \brief Applies the Riesz map of the control space to a load vector.
  //! \details This is a mass matrix solve, and is only needed when the load
  //! operator has a coupling block. The factorization is computed on the
  //! first call and reused by all the later ones.
  bool applyRiesz(Vector& X) const;

public:
  //! \brief Solves the adjoint problem for a given state.
  bool solveAdjoint(const Vector& U, Vector& P);
  //! \brief Applies the reduced Hessian to a direction in control space.
  bool applyHessian(const Vector& D, Vector& HD);

  //! \brief Evaluates the objective function.
  double computeObjective(const Vector& F, const Vector& U) const;
  //! \brief Evaluates the reduced gradient of the objective function.
  bool computeGradient(const Vector& F, const Vector& P, Vector& G) const;

private:
  SIMbase& mySim; //!< The simulator holding the state equation

  double alpha = 1.0e-6; //!< Tikhonov regularization parameter

  std::unique_ptr<SparseMatrix> myMass; //!< Control space mass matrix
  std::unique_ptr<SparseMatrix> myCoup;  //!< Coupling block of the load operator
  std::unique_ptr<SparseMatrix> myCoupT; //!< Transpose of the coupling block
  std::unique_ptr<SparseMatrix> myRiesz; //!< Factorizable copy of \ref myMass
  mutable bool rieszReady = false; //!< If \e true, \ref myRiesz is prepared

  double massScale = 1.0; //!< Scaling of the mass matrix part of the load
  double coupScale = 1.0; //!< Scaling of the coupling block of the load
  Vector myTarget;          //!< Load vector of the desired state
  double normTarget2 = 0.0; //!< Squared L2-norm of the desired state
  Vector myFixedLoad;       //!< Constant part of the state load vector

  //! Marks the degrees of freedom the control lives on, empty if all of them
  std::vector<bool> ctrlDof;

  bool adjTranspose = false; //!< If \e true, the adjoint solve is transposed
};

#endif
