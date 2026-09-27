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

#ifndef _PETSC_SOL_PARAMS_H
#define _PETSC_SOL_PARAMS_H

#include "LinSolParams.h"
#include "PETScSupport.h"

#include <set>
#include <string>
#include <vector>

class ProcessAdm;
class SettingMap;

typedef std::vector<int>         IntVec;       //!< Integer vector
typedef std::vector<IntVec>      IntMat;       //!< Integer matrix
typedef std::vector<std::string> StringVec;    //!< String vector
typedef std::vector<StringVec>   StringMat;    //!< String matrix
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix

//! \brief Schur preconditioner methods
enum SchurPrec { SIMPLE, MSIMPLER, PCD };


/*!
  \brief A geometric multigrid hierarchy for one block of the linear system.
  \details The levels are ordered with the coarsest first and the level the
  linear system itself is posed on last. The operators are optional; without
  them PETSc forms the coarse ones as the Galerkin products of the finest.
*/

struct PETScMGLevels
{
  std::vector<Mat> A; //!< Operator on each level, the finest one excluded
  std::vector<Mat> P; //!< Prolongation from level \a i to level \a i+1

  //! \brief Equations of the subdomains of each level, coarsest level first.
  //!
  //! \details A smoother given these solves the subdomains rather than the
  //! points, which is what an anisotropic mesh needs: the error a point
  //! smoother leaves behind is smooth along the strong coupling and
  //! oscillatory across it, and a subdomain spanning that coupling removes
  //! it in one go. Empty when the levels carry no subdomains of their own.
  std::vector<std::vector<std::vector<int>>> subdomains;

  //! \brief How many of the subdomains of each level are mesh lines.
  //!
  //! \details The lines of a level are followed by one subdomain holding
  //! whatever they did not reach, where they did not reach everything, and
  //! the two are not solved with the same method by default: a line is small
  //! and the remainder need not be. Empty when every subdomain is a line.
  std::vector<size_t> nLines;

  //! \brief Returns the number of levels in the hierarchy.
  size_t size() const { return P.empty() ? 0 : P.size()+1; }
};


/*!
  \brief Class for PETSc solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class PETScSolParams
{
public:
  //! \brief Default constructor.
  //! \param[in] spar The base linear solver parameters
  //! \param[in] padm The process administrator
  PETScSolParams(const LinSolParams& spar, const ProcessAdm& padm) :
    params(spar), adm(padm)
  {}

  //! \brief Destructor.
  ~PETScSolParams()
  {
  }

  //! \brief Set up preconditioner parameters for PC object
  //! \param pc Preconditioner to configure
  //! \param block Block this preconditioner applies to
  //! \param prefix PETsc param prefix for block
  //! \param blockEqs The local equations belonging to block
  //! \param setup True to setup preconditioner
  //! \param[in] mg Geometric multigrid hierarchy for this block, if any
  void setupPC(PC& pc, size_t block,
               const std::string& prefix,
               const std::set<int>& blockEqs,
               bool setup,
               const PETScMGLevels* mg = nullptr);

  //! \brief Sets up a geometric multigrid preconditioner.
  //! \param pc The preconditioner to configure
  //! \param[in] mg The hierarchy of transfer operators and level operators
  //! \param[in] map The settings to apply
  //!
  //! \details The operator of the finest level is left alone, since the KSP
  //! the preconditioner belongs to already has it. Only the coarse levels are
  //! taken from \a mg. This is public so that a preconditioner assembled
  //! outside this class, such as the inner solve of PETScSchurPC, can use a
  //! hierarchy as well.
  bool setupGeometricMG(PC& pc, const PETScMGLevels& mg,
                        const SettingMap& map);

  //! \brief Sets up an additive Schwarz smoother over the given subdomains.
  //! \param pc The smoother of one multigrid level
  //! \param[in] subdomains The equations of each subdomain on that level
  //! \param[in] iBlock Matrix block the smoother belongs to
  //! \param[in] nLines How many of the subdomains are mesh lines, the rest
  //! of them, if any, being the one holding what the lines did not reach
  void setupSubdomainSmoother(PC& pc,
                              const std::vector<std::vector<int>>& subdomains,
                              size_t iBlock, size_t nLines);

  //! \brief Obtain number of blocks
  size_t getNoBlocks() const { return params.getNoBlocks(); }

  //! \brief Obtain settings for a given block
  const LinSolParams::BlockParams& getBlock(size_t i) const { return params.getBlock(i); }

  //! \brief Get integer setting
  int getIntValue(const std::string& key) const { return params.getIntValue(key); }

  //! \brief Get integer setting
  double getDoubleValue(const std::string& key) const { return params.getDoubleValue(key); }

  //! \brief Get string setting
  std::string getStringValue(const std::string& key) const { return params.getStringValue(key); }

  //! \brief Get integer setting
  bool hasValue(const std::string& key) const { return params.hasValue(key); }

  //! \brief Returns the linear system type.
  LinAlg::LinearSystemType getLinSysType() const { return params.getLinSysType(); }

protected:
  //! \brief Set directional smoother
  //! \param[in] pc The preconditioner to add smoother for
  //! \param[in] P The preconditioner matrix
  //! \param[in] iBlock The index of the block to add smoother to
  //! \param[in] dirIndexSet The index set for the smoother
  bool addDirSmoother(PC pc, const Mat& P,
                      int iBlock, const ISMat& dirIndexSet);

  //! \brief Set ML options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The map of settings to use
  void setMLOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Set GAMG options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The map of settings to use
  void setGAMGOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Set Hypre options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The settings to apply
  void setHypreOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Setup the coarse solver in a multigrid
  //! \param[in] pc The preconditioner to set coarse solver for
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The settings to apply
  void setupCoarseSolver(PC& pc, const std::string& prefix, const SettingMap& map);

  //! \brief Setup the smoothers in a multigrid
  //! \param[in] pc The preconditioner to set coarse solver for
  //! \param[in] iBlock The index of the  block to set parameters for
  //! \param[in] dirIndexSet The index set for direction smoothers
  //! \param blockEqs The local equations belonging to block
  //! \param setup True to setup preconditioner
  //! \param[in] mg Hierarchy whose subdomains the smoothers solve over,
  //! if it brought any. Levels without subdomains are set up from the
  //! settings as usual.
  void setupSmoothers(PC& pc, size_t iBlock,
                      const ISMat& dirIndexSet,
                      const std::set<int>& blockEqs,
                      bool setup,
                      const PETScMGLevels* mg = nullptr);

  //! \brief Setup an additive Schwarz preconditioner
  //! \param pc The preconditioner to set coarse solver for
  //! \param[in] block The block the preconditioner belongs to
  //! \param[in] asmlu True to use LU subdomain solvers
  //! \param[in] smoother True if this is a smoother in multigrid
  //! \param blockEqs The local equations belonging to block
  //! \param setup True to setup preconditioner
  void setupAdditiveSchwarz(PC& pc, size_t block,
                            bool asmlu, bool smoother,
                            const std::set<int>& blockEqs,
                            bool setup);

  const LinSolParams& params; //!< Reference to linear solver parameters.
  const ProcessAdm& adm;      //!< Reference to process administrator.
};

#endif
