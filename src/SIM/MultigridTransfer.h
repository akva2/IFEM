// $Id$
//==============================================================================
//!
//! \file MultigridTransfer.h
//!
//! \date Sep 23 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Grid transfer operators for geometric multigrid.
//!
//==============================================================================

#ifndef _MULTIGRID_TRANSFER_H
#define _MULTIGRID_TRANSFER_H

#include <cstddef>
#include <memory>
#include <string>

class SIMbase;
class SparseMatrix;


namespace MG //! Utilities for geometric multigrid.
{
  /*!
    \brief Method used to compute the grid transfer operators.
  */

  enum class Transfer
  {
    CHANGE_OF_BASIS, //!< Element-local change of basis
    L2_PROJECTION    //!< Global L2 projection
  };

  /*!
    \brief Identifies an operator to build a multigrid hierarchy for.
    \details A simulator announces the operators it wants a hierarchy for
    through the MultigridProvider interface. The operator is not necessarily a
    block of the system matrix; for a Schur complement preconditioner it is
    typically an auxiliary operator which only exists on the pressure basis.
  */

  struct Operator
  {
    std::string name;      //!< Name, referred to from the linear solver input
    int         basis = 1; //!< Index of the basis the operator is defined on
    size_t      comps = 0; //!< Components of that basis, 0 means all of them
    size_t      block = 0; //!< Linear solver block the hierarchy applies to
  };

  //! \brief Builds the prolongation operator between two FE spaces.
  //! \param[in] coarse The simulator holding the coarse mesh
  //! \param[in] fine The simulator holding the fine mesh
  //! \param[in] op The operator to build the prolongation for
  //! \param[in] method The method used to compute the operator
  //! \param[out] rowsOwned Rows of the operator this process owns
  //! \param[out] colsOwned Columns of the operator this process owns
  //!
  //! \details The returned matrix maps the free DOFs of \a op on the coarse
  //! mesh to those on the fine mesh, with both numbered consecutively in
  //! order of increasing global equation number. Constrained DOFs are left
  //! out on both sides, which is the transfer a multigrid cycle needs, since
  //! it moves corrections satisfying homogeneous constraints between levels.
  //!
  //! Both simulators must be preprocessed, discretize the same geometry with
  //! the same patch layout, and \a fine must be a refinement of \a coarse.
  //! The restriction operator is the transpose of the returned matrix, which
  //! is what PETSc uses by default, so it is not built separately.
  //!
  //! Transfer::CHANGE_OF_BASIS exploits that adaptive refinement only inserts
  //! knot lines, so that the coarse space is contained in the fine one and the
  //! prolongation is the unique change-of-basis matrix between them. It is
  //! computed element by element and is exact. Transfer::L2_PROJECTION is the
  //! fallback for spaces which are not nested, where no such matrix exists.
  //! Every process builds the whole operator, since it holds every patch,
  //! and \a rowsOwned and \a colsOwned say how much of it belongs to this
  //! one. Those are the leading stretches of the two numberings on the first
  //! process and follow on from each other, the numbering being by global
  //! equation number, so they lay the operator out over the processes the
  //! way the matrices of the two levels are laid out.
  std::unique_ptr<SparseMatrix>
  prolongation(const SIMbase& coarse, const SIMbase& fine,
               const Operator& op,
               Transfer method = Transfer::CHANGE_OF_BASIS,
               int* rowsOwned = nullptr, int* colsOwned = nullptr);
}

#endif
