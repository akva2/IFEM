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

class ProcessAdm;


/*!
  \brief Interface for an operator preconditioning a solve against a Schur
  complement.

  \details Applying the Schur complement of a block system means a Krylov
  solve against it, and every iteration of that solve costs an application of
  the inverse of the leading block. Unpreconditioned it takes hundreds of
  them, so what preconditions it decides whether the whole arrangement is
  worth having.

  What makes a good one cannot be said here. It follows from the equations the
  system discretizes, and a linear algebra layer serving any saddle point
  problem knows nothing about those. So it is supplied from the outside: the
  simulator which knows what its Schur complement is like builds one of these,
  sets it on the matrix for the block it belongs to, and the Schur
  preconditioner applies it in place of the approximation it would otherwise
  factorize.

  The object is owned by whoever sets it, and has to outlive the linear
  solver it is handed to, which the simulator that assembled its operators
  does.
*/

class SchurComplementPrecond
{
public:
  //! \brief Empty destructor.
  virtual ~SchurComplementPrecond() {}

  //! \brief Sets up whatever the preconditioner needs to be applied.
  //! \param[in] adm The process administrator the solve lives on
  //!
  //! \details Called once for each linear solver the preconditioner is
  //! handed to, which may be several over the life of one of these, so it
  //! has to be able to run again and let go of what a previous call set up.
  virtual bool setUp(const ProcessAdm& adm) = 0;

  //! \brief Applies the preconditioner to a vector.
  //! \param[in] x The vector to apply it to
  //! \param[out] y The result
  virtual PetscErrorCode apply(Vec x, Vec y) = 0;

  //! \brief Returns the name PETSc is to report the preconditioner under.
  virtual const char* getName() const = 0;
};

#endif
