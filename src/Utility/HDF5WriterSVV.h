// $Id$
//==============================================================================
//!
//! \file HDF5WriterSVV.h
//!
//! \date Jul 7 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of model and results to HDF5 file for SVV.
//!
//==============================================================================

#ifndef _HDF5_WRITER_SVV_H
#define _HDF5_WRITER_SVV_H

#include "DataExporter.h"
#include "HDF5Writer.h"

class SIMbase;


/*!
  \brief Write data to a HDF5 file for SVV.

  \details The HDF5 SVV writer writes data to a HDF5 file. It supports parallel I/O.
*/

class HDF5WriterSVV : public HDF5Writer
{
public:
  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (without extension) of the data file
  //! \param[in] adm The process administrator
  HDF5WriterSVV(const std::string& name, const ProcessAdm& adm);

  //! \brief Destructor.
  virtual ~HDF5WriterSVV();

  //! \brief Returns the last time level stored in the HDF5 file.
  virtual int getLastTimeLevel();

  //! \brief Opens the file at a given time level.
  //! \param[in] level The requested time level
  virtual void openFile(int level);

  //! \brief Closes the file.
  //! \param[in] level Level we just wrote to the file
  virtual void closeFile(int level);

  //! \brief Writes data from a SIM to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the data to write
  //! \param[in] geometryUpdated Whether or not geometries should be written
  //! \param[in] prefix Field name prefix
  //!
  //! \details If prefix is non-empty and we are asked to dump secondary
  //! solutions, we assume they come from different projections
  //! \sa SIMbase::project
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool geometryUpdated, const std::string& prefix);

protected:
  bool firstStep = true;
};

#endif
