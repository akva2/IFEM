// $Id$
//==============================================================================
//!
//! \file SchurComplementPrecond.h
//!
//! \date Sep 28 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Interface for preconditioning a solve against a Schur complement.
//!
//==============================================================================

#ifndef _SCHUR_COMPLEMENT_PRECOND_H
#define _SCHUR_COMPLEMENT_PRECOND_H

#include "PETScSupport.h"

#include <vector>

class ProcessAdm;


/*!
  \brief Interface for an operator preconditioning a solve against a Schur
  complement.
*/

class SchurComplementPrecond
{
public:
  //! \brief The constructor says whether the block has a diagonal block C.
  //! \param[in] withC Whether the discretization contributes a diagonal block C
  //!
  //! \details For a two by two system
  //! \code
  //!   [ A  B^T ]
  //!   [ B  -C  ]
  //! \endcode
  explicit SchurComplementPrecond(bool withC = false) : C(withC) {}
  //! \brief Empty destructor.
  virtual ~SchurComplementPrecond() {}

  //! \brief Returns whether the block has a diagonal block C.
  bool hasC() const { return C; }

  //! \brief Sets up whatever the preconditioner needs to be applied.
  //! \param[in] adm The process administrator the solve lives on
  //! \param[in] blocks The blocks of the two by two system, in the order
  //! \code
  //!   [ 0  1 ]
  //!   [ 2  3 ]
  //! \endcode
  //! so that 0 is the leading block whose inverse the complement is taken
  //! over, 1 and 2 are the couplings and 3 is the diagonal block of the one
  //! being preconditioned.
  //!
  //! \details Called once for each linear solver the preconditioner is
  //! handed to, which may be several over the life of one of these, so it
  //! has to be able to run again and let go of what a previous call set up.
  virtual bool setUp(const ProcessAdm& adm,
                     const std::vector<Mat>& blocks) = 0;

  //! \brief Applies the preconditioner to a vector.
  //! \param[in] x The vector to apply it to
  //! \param[out] y The result
  virtual PetscErrorCode apply(Vec x, Vec y) = 0;

  //! \brief Returns the name PETSc is to report the preconditioner under.
  virtual const char* getName() const = 0;

private:
  bool C; //!< Whether the block has a diagonal block C
};

#endif
