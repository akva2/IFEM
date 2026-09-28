// $Id$
//==============================================================================
//!
//! \file PETScSchurPC.h
//!
//! \date Jun 4 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Schur-complement preconditioner using PETSc.
//!
//==============================================================================

#include "PETScMatrix.h"
#include "SchurComplementPrecond.h"

#include <memory>


/*!
  \brief Class implementing a Schur-complement preconditioner for 2-block systems.
 */


class PETScSchurPC {
public:
  //! \brief The constructor sets up the preconditioner.
  //! \param pc_init The PETSc PC to set up
  //! \param blocks The matrix blocks
  //! \param params Linear solver parameters
  //! \param adm Processes administrator
  //! \param verbosity Print the solvers that are set up if this is positive
  //! \param[in] schur What the simulator has set for this Schur complement:
  //! the preconditioner of the solve against it, and whether the block has a
  //! diagonal block C
  PETScSchurPC(PC& pc_init, const std::vector<Mat>& blocks,
               const LinSolParams::BlockParams& params, const ProcessAdm& adm,
               int verbosity = 0,
               const SchurSetup& schur = SchurSetup());

  //! \brief The destructor frees the PETSc structures.
  ~PETScSchurPC();

  //! \brief PETSc compatible function applying the Schur complement matrix.
  //! \param A Shell matrix to apply
  //! \param x Vector to apply matrix to
  //! \param y Result of matrix-vector product
  static PetscErrorCode Apply_Schur(Mat A, Vec x, Vec y);

  //! \brief PETSc compatible function applying the preconditioner of the
  //! solve against the Schur complement.
  //! \param pc Shell preconditioner to apply
  //! \param x Vector to apply preconditioner to
  //! \param y Result of preconditioner evaluation
  //!
  //! \details The work is the simulator's, this only finds it and hands the
  //! two vectors over.
  static PetscErrorCode Apply_Precond(PC pc, Vec x, Vec y);

  //! \brief PETSc compatible function applying the Schur complement preconditioner.
  //! \param pc Shell preconditioner to apply
  //! \param x Vector to apply preconditioner to
  //! \param y Result of preconditioner evaluation
  static PetscErrorCode Apply_Outer(PC pc, Vec x, Vec y);

  //! \brief Destroy a Schur preconditioner.
  //! \param pc Preconditioner to destroy
  static PetscErrorCode Destroy(PC pc);

protected:
  KSP inner_ksp; //!< The KSP for the approximation of inner matrix inverse
  KSP outer_ksp; //!< The KSP for the approximated Schur complement
  Mat outer_mat; //!< Matrix shell describing the Schur complement
  //! What the solve against that shell is preconditioned with, which is an
  //! approximation of the Schur complement with its sign turned round
  Mat prec = nullptr;

  //! What the solve against the Schur complement is preconditioned with,
  //! or nothing. Owned by whoever set it on the matrix.
  SchurComplementPrecond* precond = nullptr;
  Vec tmp;  //!< Temporary vector, laid out like the leading block
  Vec ptmp = nullptr; //!< Temporary vector, laid out like the block
  const std::vector<Mat>* m_blocks; //!< Matrix blocks in system
  //! Whether the system has a diagonal block C for the block being
  //! preconditioned, which settles whether it belongs in the complement
  bool hasC;
};
