// $Id$
//==============================================================================
//!
//! \file MultigridProvider.h
//!
//! \date Sep 23 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulator interface for geometric multigrid.
//!
//==============================================================================

#ifndef _MULTIGRID_PROVIDER_H
#define _MULTIGRID_PROVIDER_H

#include "MultigridTransfer.h"

#include <memory>
#include <vector>

class SystemMatrix;


/*!
  \brief Interface for simulators which can drive a geometric multigrid solver.

  \details A simulator implementing this interface tells the multigrid driver
  which operators it wants a hierarchy built for, and knows how to assemble
  each of them on whatever mesh it currently holds. The driver takes care of
  the meshes, the transfer operators and the plumbing into PETSc.

  The operators are deliberately not identified with blocks of the system
  matrix. A Poisson problem does want a hierarchy for the system matrix
  itself, but a Schur complement preconditioner for Stokes wants one for an
  auxiliary pressure operator, a pressure Laplacian or mass matrix, which is
  never assembled as part of the system. Letting the simulator assemble the
  operator covers both, and lets it use a different discretization on the
  coarse levels than on the fine one if it wants to.
*/

class MultigridProvider
{
public:
  //! \brief Empty destructor.
  virtual ~MultigridProvider() {}

  //! \brief Returns the operators to build multigrid hierarchies for.
  //!
  //! \details An empty list disables the geometric multigrid machinery, which
  //! is what a simulator that does not use it should return.
  virtual std::vector<MG::Operator> getMGOperators() const { return {}; }

  //! \brief Assembles one of the operators on the mesh currently held.
  //! \param[in] op The operator to assemble
  //! \return The assembled operator, or null on failure
  //!
  //! \details The operator is handed over rather than copied out: a copy of a
  //! PETSc matrix loses the block structure, and the operator of a block
  //! hierarchy is a block of the level matrix. Handing it over is also what
  //! lets the driver drop the simulator once it has taken everything else it
  //! needs off the mesh, the operator being the one thing which has to stay
  //! for as long as the hierarchy is in use.
  //!
  //! Returning null makes the driver fall back on letting PETSc form the
  //! coarse operators as Galerkin products of the finest one, which needs
  //! only the transfer operators. That is an option for a hierarchy of an
  //! actual system matrix block, but not for an auxiliary operator, which has
  //! no fine level counterpart to coarsen.
  virtual std::unique_ptr<SystemMatrix>
  assembleMGOperator(const MG::Operator& op) = 0;

  //! \brief Creates an empty simulator configured like this one.
  //!
  //! \details The driver uses this to snapshot a level of the adaptive
  //! sequence: it creates a simulator, has it read the same input file, gives
  //! it a copy of the mesh as it stands, and preprocesses it. The returned
  //! object must therefore be of the same type as this one, and constructed
  //! with the same arguments, but must not be read or preprocessed yet.
  virtual std::unique_ptr<MultigridProvider> createMGLevel() const = 0;
};

#endif
