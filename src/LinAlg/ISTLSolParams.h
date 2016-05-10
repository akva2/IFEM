// $Id$
//==============================================================================
//!
//! \file ISTLSolParams.h
//!
//! \date Mar 10 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear solver parameters for ISTL matrices.
//! \details Includes linear solver method, preconditioner
//! and convergence criteria.
//!
//==============================================================================

#ifndef _ISTL_SOLPARAMS_H
#define _ISTL_SOLPARAMS_H

#include "ISTLSupport.h"

#include <iostream>
#include <set>
#include <string>
#include <vector>


class DomainDecomposition;
class LinSolParams;
class ProcessAdm;


/*! This implements a Schur-decomposition based preconditioner for the
 *  block system
 *  [A   B]
 *  [C   D]
 *
 *  The preconditioner is
 *  [Apre  ]
 *  [     P]
 *  Here Apre is some preconditioner for A and P some preconditioner for
 *  S = D - C*diag(A)^-1*B
!*/

namespace ISTL {

class BlockPreconditioner : public Preconditioner {
public:
  // define the category
  enum {
    //! \brief The category the preconditioner is part of.
    category=Dune::SolverCategory::sequential
  };

  //! \brief Constructor
  //! \param[in] A The system matrix
  //! \param[in] dd Domain decomposition
  BlockPreconditioner(const ISTL::Mat& A, const DomainDecomposition& dd_);

  //! \brief Destructor
  virtual ~BlockPreconditioner()
  {}

  //! \brief Preprocess preconditioner
  virtual void pre(ISTL::Vec& x, ISTL::Vec& b);

  //! \brief Applies the preconditioner
  //! \param[out] v The resulting vector
  //! \param[in] d The vector to apply the preconditioner to
  virtual void apply(ISTL::Vec& v, const ISTL::Vec& d);

  //! \brief Post-process function
  virtual void post(ISTL::Vec& x);

  //! \brief Obtain reference to a block preconditioner
  std::unique_ptr<ISTL::Preconditioner>& getBlockPre(size_t block)
  {
    return blockPre[block];
  }

  //! \brief Obtain reference to a block matrix
  ISTL::Mat& getBlock(size_t block)
  {
    return blocks[block];
  }

  //! \brief Obtain matrix adaptor for a block matrix
  ISTL::Operator& getBlockOp(size_t block)
  {
    if (!blockOp[block])
      blockOp[block].reset(new ISTL::Operator(blocks[block]));

    return *blockOp[block];
  }

protected:
  //! \brief Build block from block equations.
  static void extractBlock(ISTL::Mat& B, const ISTL::Mat& A,
                           const std::set<int>& eqs_row,
                           const std::set<int>& eqs_col);

  //! \brief Find A = B - C
  static void subtractMatrices(ISTL::Mat& A,
                               const ISTL::Mat& B,
                               const ISTL::Mat& C);

  std::vector<std::unique_ptr<ISTL::Preconditioner>> blockPre; //!< The preconditioners
  std::vector<std::unique_ptr<ISTL::Operator>> blockOp; //!< The matrix adaptors for the blocks
  std::vector<ISTL::Mat> blocks; //!< Matrix blocks
  const DomainDecomposition& dd; //!< Domain decomposition
};

}


/*!
  \brief Class for ISTL solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class ISTLSolParams
{
public:
  //! \brief Constructor.
  //! \param spar Linear solver parameters to use
  //! \param padm Process administrator to use
  ISTLSolParams(const LinSolParams& spar, const ProcessAdm& padm) :
    solParams(spar), adm(padm)
  {}


  //! \brief Setup solver and preconditioner.
  //! \param A Matrix to use
  //! \return tuple with (solver, preconditioner, op)
  std::tuple<std::unique_ptr<ISTL::InverseOperator>,
             std::unique_ptr<ISTL::Preconditioner>,
             std::unique_ptr<ISTL::Operator>>
    setupPC(ISTL::Mat& A);

  //! \brief Obtain linear solver parameters.
  const LinSolParams& get() const { return solParams; }

protected:
  //! \brief Internal helper function for setting up a preconditioner.
  //! \param A Matrix to construct preconditioner for
  //! \param op Matrix adaptor to use (used with AMG)
  //! \param block Block to read settings from
  ISTL::Preconditioner* setupPCInternal(ISTL::Mat& A,
                                        ISTL::Operator& op,
                                        size_t block);

  const LinSolParams& solParams; //!< Reference to linear solver parameters.
  const ProcessAdm& adm;      //!< Reference to process administrator.
};

#endif
