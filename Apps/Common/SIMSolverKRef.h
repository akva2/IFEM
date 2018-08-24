// $Id$
//==============================================================================
//!
//! \file SIMSolverKRef.h
//!
//! \date Feb 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Stationary solver class template for K-refinement.
//!
//==============================================================================

#ifndef _SIM_SOLVER_KREF_H_
#define _SIM_SOLVER_KREF_H_

#include "SIMSolver.h"
#include "SIMMxV.h"
#include "SIMPC.h"
#include "Profiler.h"


/*!
  \brief Template class for stationary K-refinement simulator drivers.
  \details This template can be instanciated over any type implementing the
  ISolver interface. It provides a solver loop with data output.
*/

template<class T1> class SIMSolverKRef : public SIMSolver<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMSolverKRef(T1& s1, T1& s2) : 
    SIMSolver<T1>(s1), mxv(s1), S2(s2), pc(s1, s2)
  {
    this->S1.opt.solver = SystemMatrix::PETSC;
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverKRef() {}

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, 
                           const char* heading = nullptr, bool = false)
  {
    int geoBlk= 0, nBlock=0;
    if (!this->saveState(geoBlk,nBlock,true,infile,false))
      return 2;

    this->printHeading(heading);

    mxv.setPC(&pc);
    pc.setupPC();
    S2.assembleStep(this->tp);
    Vector csol(S2.getNoDOFs());
    S2.solveSystem(csol);
    StdVector tmp(csol.size());
    for (size_t i = 1; i <= csol.size(); ++i)
      tmp(i) = csol(i);
    StdVector tmp2(pc.N.rows());
    pc.N.multiply(tmp, tmp2);
    pc.B.solve(tmp2);
    PETScVector tmp3(this->S1.getProcessAdm(),
                     this->S1.getSAM()->getNoEquations());
    const int* meqn = this->S1.getSAM()->getMEQN();
    for (size_t i = 1; i <= this->S1.getPatch(1)->getNoNodes(); ++i) {
      int eq = meqn[i-1];
      if (eq > 0)
        tmp3(eq) = tmp2(i);
    }
    static_cast<PETScMatrix*>(this->S1.getAlgEqSystem()->getMatrix(0))->hack = &tmp3;
    tmp3.beginAssembly();
    tmp3.endAssembly();

    if (!mxv.solveStep(this->tp))
      return 1;
    else if (!this->saveState(geoBlk, nBlock,false,infile,true))
      return 2;

    return 0;
  }

  int initSystems(const char* infile)
  {
    this->S1.addRef = true;
    if (!this->S1.read(infile) || !this->S2.read(infile) || !this->read(infile))
      return 1;

    if (!this->S1.preprocess() ||
        !this->S2.preprocess())
      return 2;

    this->S1.setQuadratureRule(this->S1.opt.nGauss[0],true);
    this->S2.setQuadratureRule(this->S2.opt.nGauss[0],true);
    const_cast<LinSolParams*>(this->S1.getSolParams())->addValue("matrixfree","1");
    return this->S1.initSystem(this->S1.opt.solver,1,1,0,true) &&
           this->S2.initSystem(this->S2.opt.solver,1,1,0,true) ? 0 : 3;
  }

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"krefinement"))
      return this->SIMSolver<T1>::parse(elem);

    utl::getAttribute(elem,"cycles",cycles);
    return true;
  }
  
  size_t cycles;   //!< Number of k-refinement cycles
  SIMMxV<T1> mxv; //!< Higher-order solver wrapped as a matrix-vector product.
  T1& S2;
  SIMKCyclePC<T1> pc;   //!< Lower-order solver wrapped as a preconditioner
};


/*!
  \brief Template class for specializing solveSystem for the inner solver.
  \details This template can be instantiated over any type implementing the
  SIMbase interface. It reimplements the solveSystem method to suppress some
  output.
*/

template<class T1> class SIMInnerKCycle : public T1
{
  using SetupProps = typename T1::SetupProps;
public:
  SIMInnerKCycle(const SetupProps& props) : T1(props) {}

  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[out] rCond Reciprocal condition number
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] newLHS If \e false, reuse the LHS-matrix from previous call.
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  virtual bool solveSystem(Vector& solution, int printSol, double* rCond,
                           const char* compName = "displacement",
                           bool newLHS = true, size_t idxRHS = 0) override
  {
    SystemMatrix* A = this->myEqSys->getMatrix();
    SystemVector* b = this->myEqSys->getVector(idxRHS);
    if (!A) std::cerr <<" *** SIMKCycleInner::solveSystem: No LHS matrix."<< std::endl;
    if (!b) std::cerr <<" *** SIMKCycleInner::solveSystem: No RHS vector."<< std::endl;
    if (!A || !b) return false;

    // Solve the linear system of equations
    utl::profiler->start("Equation solving");
    bool status = A->solve(*b,newLHS,nullptr);
    utl::profiler->stop("Equation solving");

    // Expand solution vector from equation ordering to DOF-ordering
    if (status)
      status = this->mySam->expandSolution(*b,solution);

    return status;
  }
};


#endif
