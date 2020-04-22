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
  //! \brief The constructor forwards to the parent class constructor.
  AdaptiveISolver(T1& sim, bool sa) : AdaptiveSIM(sim,sa), model(sim) {}

protected:
  //! \brief Assembles and solves the linearized FE equation system.
  virtual bool assembleAndSolveSystem()
  {
    TimeStep dummy;
    model.init(dummy);
    if (model.solveStep(dummy))
      solution = model.getSolutions();
    else
      return false;

    return true;
  }

  T1& model; //!< Reference to the actual sim
};


/*!
  \brief Template class for stationary adaptive simulator drivers.
  \details This template can be instantiated over any type implementing the
  ISolver interface. It provides an adaptive loop with data output.
*/

template<class T1, class AdapSim>
class SIMSolverAdapImpl : public SIMSolverStat<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMSolverAdapImpl(T1& s1) : SIMSolverStat<T1>(s1), aSim(s1,false)
  {
    this->S1.setSol(&aSim.getSolution());
    this->S1.opt.saveNorms = true;
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverAdapImpl() {}

  //! \brief Reads solver data from the specified input file.
  virtual bool read(const char* file) { return this->SIMadmin::read(file); }

  //! \brief Solves the problem on a sequence of adaptively refined meshes.
  virtual int solveProblem(char* infile, const char* = nullptr, bool = false)
  {
    if (!aSim.initAdaptor())
      return 1;

    if (SIMSolverStat<T1>::exporter) {
      SIMSolverStat<T1>::exporter->setFieldValue(exporterName, &this->S1,
                                                 &aSim.getSolution(),
                                                 &aSim.getProjections(),
                                                 &aSim.getEnorm());
      if (!this->S1.opt.project.empty()) {
        std::vector<std::string> pref;
        for (const auto& it : this->S1.opt.project)
          pref.push_back(it.second);
        SIMSolverStat<T1>::exporter->setNormPrefixes(pref);
      }
    }

    for (int iStep = 1; aSim.adaptMesh(iStep); iStep++)
      if (!aSim.solveStep(infile,iStep))
        return 1;
      else if (!aSim.writeGlv(infile,iStep))
        return 2;
      else if (SIMSolverStat<T1>::exporter)
        SIMSolverStat<T1>::exporter->dumpTimeLevel(nullptr,true);

    return 0;
  }

  //! \brief Set name of data exporter registration to use.
  void setExporterName(const std::string& name) { exporterName = name; }

protected:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* kyw, std::istream& is) { return aSim.parse(kyw,is); }
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) { return aSim.parse(elem); }

  AdapSim aSim; //!< Adaptive simulation driver
  std::string exporterName = "u"; //!< Name for data exporter registration to use
};


//! Convenience alias template
template<class T1>
using SIMSolverAdap = SIMSolverAdapImpl<T1,AdaptiveSIM>;


/*!
  \brief Template class for stationary adaptive, staggered simulator drivers.
*/

template<class T1>
class SIMSolverStaggeredAdap : public SIMSolverStat<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMSolverStaggeredAdap(T1& s1) : SIMSolverStat<T1>(s1)
  {
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverStaggeredAdap() {}

  //! \brief Solves the problem on a sequence of adaptively refined meshes.
  virtual int solveProblem(char* infile, const char* = nullptr, bool = false)
  {
    // Save FE model to VTF and HDF5 for visualization
    // Optionally save the initial configuration also
    int geoBlk = 0, nBlock = 0;
    if (!this->saveState(geoBlk,nBlock,true,infile,false))
      return 2;

    for (int iStep = 1; this->S1.adaptMesh(iStep); iStep++)
      if (!this->S1.solveStep(infile,iStep))
        return 1;
      else if (!this->saveState(geoBlk,nBlock))
        return 2;

    return 0;
  }

protected:
  //! \brief Saves geometry and results to VTF and HDF5 for current time step.
  bool saveState(int& geoBlk, int& nBlock, bool newMesh = false,
                 char* infile = nullptr, bool saveRes = true)
  {
    return true;
//    if (newMesh && !this->S1.saveModel(infile,geoBlk,nBlock))
//      return false;

//    if (saveRes && !this->S1.saveStep(tp,nBlock))
//      return false;

//    if (saveRes && SIMSolverStat<T1>::exporter) {
//      HDF5Restart::SerializeData data;
//      if (restartAdm && restartAdm->dumpStep(tp) && this->serialize(data))
//        if (!restartAdm->writeData(tp,data))
//          return false;

//      return SIMSolverStat<T1>::exporter->dumpTimeLevel(&tp,newMesh);
//    }

//    return true;
  }
};

#endif
