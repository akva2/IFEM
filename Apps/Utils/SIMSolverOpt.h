// $Id$
//==============================================================================
//!
//! \file SIMSolverOpt.h
//!
//! \brief Stationary solver class template for optimal control problems.
//!
//==============================================================================

#ifndef _SIM_SOLVER_OPT_H_
#define _SIM_SOLVER_OPT_H_

#include "OptimizationDriver.h"
#include "SIMSolver.h"

#include "tinyxml2.h"

#include <cstring>
#include <iostream>


/*!
  \brief Template class for adjoint-based optimal control simulator drivers.

  \details This driver solves a PDE-constrained optimization problem on a
  single, fixed discretization. The optimization algorithm itself lives in
  OptimizationDriver, which is shared with the adaptive driver in
  SIMSolverOptAdap.h.
*/

template<class T1>
class SIMSolverOpt : public SIMSolverStat<T1>
{
public:
  //! \brief Optimization method, see OptimizationDriver::Method.
  using Method = typename OptimizationDriver<T1>::Method;

  //! \brief The constructor initializes the reference to the simulator.
  explicit SIMSolverOpt(T1& s1)
    : SIMSolverStat<T1>(s1,"Adjoint / Optimal Control Driver"), optDrv(s1) {}

  //! \brief Empty destructor.
  virtual ~SIMSolverOpt() = default;

  //! \brief Defines the optimization method to use.
  void setMethod(Method m) { optDrv.setMethod(m); }
  //! \brief Defines the absolute gradient norm tolerance.
  void setTolerance(double t) { optDrv.setTolerance(t); }
  //! \brief Defines the maximum number of optimization iterations.
  void setMaxIterations(int it) { optDrv.setMaxIterations(it); }
  //! \brief Defines box constraints on the control.
  void setBounds(double low, double high) { optDrv.setBounds(low,high); }

  //! \brief Reads solver parameters from the specified input file.
  bool read(const char* file) override { return this->SIMadmin::read(file); }

  //! \brief Solves the optimal control problem.
  //! \param[in] infile File name used to construct the VTF-file name from
  //! \param[in] heading Additional heading to print before the iterations
  //! \return Zero on success, a positive error code if one of the simulator
  //! operations failed, and 12 if the optimization ran out of iterations or
  //! the line search gave up before the convergence criterion was met
  int solveProblem(char* infile, const char* heading = nullptr) override
  {
    int geoBlk = 0, nBlock = 0;
    if (!this->S1.saveModel(infile,geoBlk,nBlock))
      return 1;

    this->printHeading(heading ? heading :
                       "Solving PDE-constrained optimal control problem");

    if (!this->S1.initOptimalControl())
      return this->error("Failed to initialize the optimal control problem",2);

    Vector F;
    if (!this->S1.getControl(F))
      return this->error("Failed to obtain the initial control",3);

    Vector U, P;
    int status = optDrv.optimize(F,U,P);
    if (status != 0 && status != 12)
      return status;

    this->S1.printErrors(F,U,P);

    if (!this->S1.saveOptimalResults(F,U,P,infile,geoBlk,nBlock))
      return this->error("Failed to save the optimal solution",10);

    if (this->exporter && !this->exporter->dumpTimeLevel())
      return 11;

    return status;
  }

protected:
  //! \brief Parses the optimization parameters from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    if (strcasecmp(elem->Value(),"optimalcontrol") &&
        strcasecmp(elem->Value(),"adjoint"))
      return this->SIMSolverStat<T1>::parse(elem);

    return optDrv.parse(elem);
  }

  //! \brief Prints an error message and returns the given error code.
  int error(const char* msg, int code) const
  {
    std::cerr <<" *** SIMSolverOpt: "<< msg <<"."<< std::endl;
    return code;
  }

  OptimizationDriver<T1> optDrv; //!< The reduced-space optimization algorithm
};


//! \brief Convenience alias for adjoint simulations.
template<class T1> using SIMSolverAdjoint = SIMSolverOpt<T1>;

#endif
