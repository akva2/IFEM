// $Id$
//==============================================================================
//!
//! \file SIMSolverOptAdap.h
//!
//! \brief Adaptive solver class template for optimal control problems.
//!
//==============================================================================

#ifndef _SIM_SOLVER_OPT_ADAP_H_
#define _SIM_SOLVER_OPT_ADAP_H_

#include "OptimizationDriver.h"
#include "SIMSolverAdap.h"

#include "tinyxml2.h"

#include <cstring>


/*!
  \brief Adaptive driver for PDE-constrained optimal control problems.

  \details This combines the reduced-space optimizer with the adaptive
  simulation driver, so that the control problem is solved on a sequence of
  meshes that are refined towards the objective functional rather than towards
  a norm of the solution.

  The two fit together without any new error estimation machinery. Each cycle
  solves the optimal control problem to convergence and hands back the state
  \b and the adjoint as the two solution fields. The adaptive driver already
  forms its refinement indicator as the product of the recovery-based energy
  error of the first field and that of the second, which for this pair is
  precisely the dual-weighted-residual estimate of the error in the objective:
  the state residual weighted by the adjoint, plus the adjoint residual
  weighted by the state. Since the optimizer cannot have converged without
  computing both fields, the estimate costs no additional solves.

  The control is carried across meshes alongside the state and the adjoint, so
  every cycle after the first is warm started from the previous optimum and
  typically needs only one or two iterations.

  This is the discrete analogue of the error representation in Becker, Kapp and
  Rannacher, "Adaptive finite element methods for optimal control of partial
  differential equations: Basic concept", SIAM J. Control Optim. 39 (2000).
*/

template<class T1>
class AdaptiveOptSolver : public AdaptiveSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param sim The simulator to solve the control problem for
  //! \param[in] sa If \e true, this driver owns the input file parsing
  AdaptiveOptSolver(T1& sim, bool sa)
    : AdaptiveSIM(sim,sa), model(sim), optDrv(sim) {}

  //! \brief Empty destructor.
  virtual ~AdaptiveOptSolver() = default;

  using AdaptiveSIM::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    if (!strcasecmp(elem->Value(),"optimalcontrol") ||
        !strcasecmp(elem->Value(),"adjoint"))
      return optDrv.parse(elem);

    return this->AdaptiveSIM::parse(elem);
  }

  //! \brief Gives access to the optimizer, for command-line overrides.
  OptimizationDriver<T1>& getOptimizer() { return optDrv; }

  //! \brief Refines the mesh based on the goal-oriented indicator.
  //! \param[in] iStep Refinement step counter
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  //! \details The control is mapped onto the refined mesh alongside the state
  //! and the adjoint, so that the next cycle can be warm started from it.
  bool adaptMesh(int iStep, std::streamsize outPrec = 0)
  {
    if (iStep > 1)
    {
      transfered = solution;
      transfered.push_back(control);
    }

    return this->AdaptiveSIM::adaptMesh(iStep,transfered,outPrec);
  }

protected:
  //! \brief Solves the optimal control problem on the current mesh.
  //! \details The state and the adjoint are returned as the two solution
  //! fields, which is what lets the adaptive driver form the goal-oriented
  //! indicator from them without knowing that it is a control problem.
  bool assembleAndSolveSystem() override
  {
    // On a refined mesh the control transferred from the previous cycle is
    // the starting point, so the optimizer is warm started
    if (transfered.size() > 2)
      control = transfered.back();

    if (!model.initOptimalControl())
      return false;

    if (control.empty() && !model.getControl(control))
      return false;

    Vector U, P;
    int status = optDrv.optimize(control,U,P);
    if (status != 0 && status != 12)
      return false;

    model.printErrors(control,U,P);

    solution.resize(2);
    solution[0] = U;
    solution[1] = P;

    return true;
  }

  T1& model;          //!< Reference to the control problem simulator
  Vectors transfered; //!< Fields mapped from the previous mesh
  Vector control;     //!< Control vector, carried across the refinements

private:
  OptimizationDriver<T1> optDrv; //!< The reduced-space optimization algorithm
};


//! \brief Convenience alias template for the adaptive control driver.
template<class T1>
using SIMSolverOptAdap = SIMSolverAdapImpl<T1,AdaptiveOptSolver<T1>>;

#endif
