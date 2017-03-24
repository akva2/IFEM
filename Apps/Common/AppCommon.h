// $Id$
//==============================================================================
//!
//! \file AppCommon.h
//!
//! \date Nov 06 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Common helper templates for applications.
//!
//==============================================================================

#ifndef _APP_COMMON_H_
#define _APP_COMMON_H_

#include "XMLInputBase.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "IFEM.h"


namespace SIM
{
  //! \brief Base class for input file pre-parsing in applications.
  class AppXMLInputBase : public XMLInputBase
  {
  public:
    //! \brief Default constructor.
    AppXMLInputBase() : dim(3) {}

  protected:
    //! \brief Parses a data section from an XML element.
    virtual bool parse(const TiXmlElement* elem);

  public:
    int dim; //!< Dimensionality of simulation
    std::string restartfile; //!< File to read restart data from
    int restartlevel = -1; //!< Restart level
  };

  //! \brief Handles application restarts.
  //! \param[in] simulator The top SIMbase instance of your application
  //! \param[in] solver The SIMSolver instance of your application
  //! \param[in] restartfile The file to read from
  //! \param[in] level The level to read from
  template<class Simulator, class Solver>
  void handleRestart(Simulator& simulator, Solver& solver,
                     const std::string& restartfile, int level=-1)
  {
    HDF5Writer hdf(restartfile,solver.getProcessAdm(),true);
    DataExporter::SerializeData data;
    int alevel = hdf.readRestartData(data, level);
    if (alevel >= 0) {
      IFEM::cout << "\n=== Restarting from a serialized state ==="
                 << "\n  file  = " << restartfile
                 << "\n  level = " << alevel << std::endl;
      solver.deSerialize(data);
      simulator.deSerialize(data);
    }
  }

  //! \brief Handles application data output.
  //! \param[in] simulator The top SIMbase instance of your application
  //! \param[in] solver The SIMSolver instance of your application
  //! \param[in] hdf5file The file to save to
  //! \param[in] append Whether or not to append to file
  //! \param[in] interval The stride in the output file
  template<class Simulator, class Solver>
  DataExporter* handleDataOutput(Simulator& simulator, Solver& solver,
                                 const std::string& hdf5file,
                                 bool append = false,
                                 int interval = 1,
                                 int restartInterval = 0)
  {
    DataExporter* writer = new DataExporter(true,interval,restartInterval);
    XMLWriter* xml = new XMLWriter(hdf5file,solver.getProcessAdm());
    HDF5Writer* hdf = new HDF5Writer(hdf5file,solver.getProcessAdm(),append);
    writer->registerWriter(xml);
    writer->registerWriter(hdf);
    simulator.registerFields(*writer);
    IFEM::registerCallback(*writer);
    return writer;
  }
}

#endif
