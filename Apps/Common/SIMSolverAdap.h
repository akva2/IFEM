// $Id$
//==============================================================================
//!
//! \file SIMSolverAdap.h
//!
//! \date Feb 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Stationary, adaptive SIM solver class template.
//!
//==============================================================================

#ifndef _SIM_SOLVER_ADAP_H_
#define _SIM_SOLVER_ADAP_H_

#include "SIMSolver.h"
#include "AdaptiveSIM.h"
#include "TimeStep.h"

/*!
  \brief Adaptive simulator driver using the ISolver interface.
*/

template<class T1>
class AdaptiveISolver : public AdaptiveSIM
{
public:
  //! \brief The constructor initializes default adaptation parameters.
  //! \param sim The FE model
  //! \param sa If \e true, this is a stand-alone driver
  AdaptiveISolver(T1& sim, bool sa=true) :
    AdaptiveSIM(sim, sa), model(sim) {}

  //! \brief Assemble and solve the equation system.
  virtual bool solveSystem()
  {
    TimeStep dummy;
    model.init(dummy);
    bool result = model.solveStep(dummy);
    if (result)
      solution = model.getSolutions();

    return result;
  }

protected:
  T1& model; //!< Reference to the actual sim
};


/*!
  \brief Template class for stationary adaptive simulator drivers.
  \details This template can be instanciated over any type implementing the
  ISolver interface. It provides an adaptive loop with data output.
*/

template<class T1, class AdapSim>
class SIMSolverAdapImpl : public SIMSolver<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMSolverAdapImpl(T1& s1) : SIMSolver<T1>(s1), aSim(s1,false)
  {
    this->S1.setSol(&aSim.getSolution());
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverAdapImpl() {}

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, const char* heading, bool = false)
  {
    if (!aSim.initAdaptor())
      return 1;

    this->printHeading(heading);

    for (int iStep = 1; aSim.adaptMesh(iStep); iStep++)
      if (!aSim.solveStep(infile,iStep))
        return 1;
      else if (!aSim.writeGlv(infile,iStep))
        return 2;
      else if (SIMSolver<T1>::exporter)
        SIMSolver<T1>::exporter->dumpTimeLevel(nullptr,true);

    return 0;
  }

protected:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyw, std::istream& is)
  {
    return this->SIMSolver<T1>::parse(keyw,is) && aSim.parse(keyw,is);
  }
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    return this->SIMSolver<T1>::parse(elem) && aSim.parse(elem);
  }

  AdapSim aSim; //!< Adaptive simulation driver
};

template<class T1>
using SIMSolverAdap = SIMSolverAdapImpl<T1, AdaptiveSIM>; //!< Convenience typedef

#endif
