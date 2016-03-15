// $Id$
//==============================================================================
//!
//! \file PETScSolParams.h
//!
//! \date Mar 10 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//! \details Includes linear solver method, preconditioner
//! and convergence criteria.
//!
//==============================================================================

#ifndef _PETSCSOLPARAMS_H
#define _PETSCSOLPARAMS_H

#if HAS_PETSC

#include "LinSolParams.h"
#include "PETScSupport.h"

#include <array>
#include <iostream>
#include <string>
#include <vector>

class ProcessAdm;


typedef std::vector<int>         IntVec;       //!< Integer vector
typedef std::vector<IntVec>      IntMat;       //!< Integer matrix
typedef std::vector<std::string> StringVec;    //!< String vector
typedef std::vector<StringVec>   StringMat;    //!< String matrix
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix

//! \brief Schur preconditioner methods
enum SchurPrec { SIMPLE, MSIMPLER, PCD };


/*!
  \brief Class for PETSc solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class PETScSolParams
{
public:
  //! \brief Set linear solver parameters for KSP object
  static void setParams(const LinSolParams& params,
                        const ProcessAdm& adm,
                        KSP& ksp, PetscRealVec& coords,
                        ISMat& dirIndexSet, int nsd);

  //! \brief Set directional smoother
  //! \param[in] PC The preconditioner to add smoother for
  //! \param[in] P The preconditioner matrix
  //! \param[in] params The block parameters
  //! \param[in] iBlock The index of the block to add smoother to
  //! \param[in] dirIndexSet The index set for the smoother
  static bool addDirSmoother(PC pc, Mat P, const LinSolParams::BlockParams& block,
                             int iBlock, ISMat& dirIndexSet);

  //! \brief Set ML options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The map of settings to use
  static void setMLOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Set GAMG options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] block The index of the block to set parameters for
  static void setGAMGOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Set Hypre options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The settings to apply
  static void setHypreOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Setup the coarse solver in a multigrid
  //! \param[in] PC The preconditioner to set coarse solver for
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The settings to apply
  static void setupCoarseSolver(PC& pc, const std::string& prefix, const SettingMap& map);

  //! \brief Setup the smoothers in a multigrid
  //! \param[in] PC The preconditioner to set coarse solver for
  //! \param[in] params The linear solver parameters
  //! \param[in] iBlock The index of the  block to set parameters for
  //! \param[in] dirIndexSet The index set for direction smoothers
  //! \param[in] locSubdDofs Local subdomain DOFs for ASM preconditioners
  //! \param[in] subdDofs Subdomain DOFs for ASM preconditioners
  static void setupSmoothers(PC& pc, const LinSolParams& params, size_t iBlock,
                             ISMat& dirIndexSet, const ProcessAdm& adm);

  //! \brief Setup an additive Schwarz preconditioner
  //! \param[in] PC The preconditioner to set coarse solver for
  //! \param[in] overlap The overlap
  //! \param[in] nx Subdomains per patch in x
  //! \param[in] ny Subdomains per patch in y
  //! \param[in] nz Subdomains per patch in z
  //! \param[in] asmlu True to use LU subdomain solvers
  //! \param[in] adm The process administrator (including the domain decomposition)
  //! \param[in] locSubdDofs Local subdomain DOFs for ASM preconditioners
  //! \param[in] subdDofs Subdomain DOFs for ASM preconditioners
  //! \param[in] smoother True if this is a smoother in multigrid
  static void setupAdditiveSchwarz(PC& pc, int overlap, 
                                   int nx, int ny, int nz,
                                   bool asmlu, bool smoother,
                                   const ProcessAdm& adm);
};

#endif
#endif
