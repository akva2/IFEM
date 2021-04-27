// $Id$
//==============================================================================
//!
//! \file HDF5WriterSVV.C
//!
//! \date Apr 27 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of model and results to HDF5 file for SVV.
//!
//==============================================================================

#include "HDF5WriterSVV.h"

#include "ASMbase.h"
#include "IntegrandBase.h"
#include "MatVec.h"
#include "ProcessAdm.h"
#include "SIMbase.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif


HDF5WriterSVV::HDF5WriterSVV (const std::string& name, const ProcessAdm& adm)
  : HDF5Writer(name+"_svv",adm,false)
{
}

HDF5WriterSVV::~HDF5WriterSVV()
{
  if (m_file != -1)
    H5Fclose(m_file);
}


int HDF5WriterSVV::getLastTimeLevel ()
{
  return 0;
}


void HDF5WriterSVV::openFile(int level)
{
#ifdef HAS_HDF5
  if (m_file != -1)
    return;

  if (!HDF5Base::openFile(m_flag, true))
    return;

  if (!checkGroupExistence(m_file,"/0"))
    H5Gclose(H5Gcreate2(m_file,"/0",0,H5P_DEFAULT,H5P_DEFAULT));
#endif
}


void HDF5WriterSVV::closeFile(int level)
{
#ifdef HAS_HDF5
  if (m_file != -1)
    H5Fflush(m_file,H5F_SCOPE_GLOBAL);
#endif
}


void HDF5WriterSVV::writeSIM (int level, const DataEntry& entry,
                              bool geometryUpdated, const std::string& prefix)
{
  if (firstStep) {
    this->writeSIM2(0, entry, geometryUpdated, prefix, true);
    firstStep = false;
  }

  auto status = H5Fstart_swmr_write(m_file);

  return this->HDF5Writer::writeSIM2(0, entry, geometryUpdated, prefix, false);
}
