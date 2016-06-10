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
#include <dune/istl/operators.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/common/version.hh>

#include "LinSolParams.h"

#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "IFEM.h"


class DomainDecomposition;
class ProcessAdm;


/*! This implements a Schur-decomposition based preconditioner for the
 *  block system
 *  [A   B]
 *  [C   D]
 *
 *  The preconditioner is
 *  [Apre B]
 *  [     P]
 *  Here Apre is some preconditioner for A and P some preconditioner for
 *  S = D - C*diag(A)^-1*B. The B block may be dropped to make the preconditioner diagonal.
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
  BlockPreconditioner(const ISTL::Mat& A, const DomainDecomposition& dd_,
                      const std::string& schurType);

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
  Mat& getBlock(size_t block)
  {
    return blocks[block];
  }

  //! \brief Obtain matrix adaptor for a block matrix
  ISTL::MatrixAdapter& getBlockOp(size_t block)
  {
    if (!blockOp[block])
      blockOp[block].reset(new ISTL::MatrixAdapter(blocks[block]));

    return *blockOp[block];
  }

protected:
  //! \brief Build block from block equations.
  static void extractBlock(ISTL::Mat& B, const ISTL::Mat& A,
                           const std::set<int>& eqs_row,
                           const std::set<int>& eqs_col,
                           const std::map<int,int>& eqs_col_g2l);

  //! \brief Find A = B - C
  static void subtractMatrices(ISTL::Mat& A,
                               const ISTL::Mat& B,
                               const ISTL::Mat& C);

  std::vector<std::unique_ptr<ISTL::Preconditioner>> blockPre; //!< The preconditioners
  std::vector<std::unique_ptr<ISTL::MatrixAdapter>> blockOp; //!< The matrix adaptors for the blocks
  std::vector<ISTL::Mat> blocks; //!< Matrix blocks
  const DomainDecomposition& dd; //!< Domain decomposition
};

#ifdef HAVE_MPI
/**
 * \brief An overlapping schwarz operator.
 *
 * This operator represents a parallel matrix product using
 * sequential data structures together with a parallel index set
 * describing an overlapping domain decomposition and the communication.
 * \tparam M The type of the sequential matrix to use,
 * e.g. BCRSMatrix or another matrix type fulfilling the
 * matrix interface of ISTL.
 * \tparam X The type of the sequential vector to use for the left hand side,
 * e.g. BlockVector or another type fulfilling the ISTL
 * vector interface.
 * \tparam Y The type of the sequential vector to use for the right hand side,
 * e..g. BlockVector or another type fulfilling the ISTL
 * vector interface.
 * \tparam C The type of the communication object.
 * This must either be OwnerOverlapCopyCommunication or a type
 * implementing the same interface.
 *
 * This is a modified version which does summation of the vector after evaluation.
 */
template<class M, class X, class Y, class C>
class OverlappingSchwarzOperator : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
  //! \brief The type of the matrix we operate on.
  //!
  //! E.g. BCRSMatrix or another matrix type fulfilling the
  //! matrix interface of ISTL
  typedef M matrix_type;
  //! \brief The type of the domain.
  //!
  //! E.g. BlockVector or another type fulfilling the ISTL
  //! vector interface.
  typedef X domain_type;
  //! \brief The type of the range.
  //!
  //! E.g. BlockVector or another type fulfilling the ISTL
  //! vector interface.
  typedef Y range_type;
  //! \brief The field type of the range
  typedef typename X::field_type field_type;
  //! \brief The type of the communication object.
  //!
  //! This must either be OwnerOverlapCopyCommunication or a type
  //! implementing the same interface.
  typedef C communication_type;

  enum {
    //! \brief The solver category.
    category=Dune::SolverCategory::overlapping
  };

  /**
   * @brief constructor: just store a reference to a matrix.
   *
   * @param A The assembled matrix.
   * @param com The communication object for syncing overlap and copy
   * data points. (E.~g. OwnerOverlapCopyCommunication )
   */
  OverlappingSchwarzOperator (const matrix_type& A, const communication_type& com)
    : _A_(A), communication(com)
  {}

  //! apply operator to x:  \f$ y = A(x) \f$
  virtual void apply (const X& x, Y& y) const
  {
    _A_.umv(x,y);
    communication.addOwnerOverlapToAll(y,y);
  }

  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
  {
    Y y2(y.size());
    apply(x,y2);
    y.axpy(alpha,y2);
  }

  //! get the sequential assembled linear operator.
  virtual const matrix_type& getmat () const
  {
    return _A_;
  }

private:
  const matrix_type& _A_;
  const communication_type& communication;
};

typedef OverlappingSchwarzOperator<ISTL::Mat,ISTL::Vec,ISTL::Vec,Dune::OwnerOverlapCopyCommunication<int,int>> ParMatrixAdapter;


/**
 * @brief Block parallel preconditioner.
 *
 * This is essentially a wrapper that takes a sequential
 * preconditioner. In each step the sequential preconditioner
 * is applied and then all owner data points are updated on
 * all other processes.
 * \tparam M The type of the sequential matrix to use,
 * e.g. BCRSMatrix or another matrix type fulfilling the
 * matrix interface of ISTL.
 * \tparam X The type of the sequential vector to use for the left hand side,
 * e.g. BlockVector or another type fulfilling the ISTL
 * vector interface.
 * \tparam Y The type of the sequential vector to use for the right hand side,
 * e..g. BlockVector or another type fulfilling the ISTL
 * vector interface.
 * \tparam C The type of the communication object.
 * This must either be OwnerOverlapCopyCommunication or a type
 * implementing the same interface.
 * \tparam The type of the sequential preconditioner to use
 * for approximately solving the local matrix block consisting of unknowns
 * owned by the process. Has to implement the Preconditioner interface.
 */
template<class X, class Y, class C, class T=Dune::Preconditioner<X,Y> >
class pASM : public Dune::Preconditioner<X,Y> {
  public:
  //! \brief The domain type of the preconditioner.
  //!
  //! E.g. BlockVector or another type fulfilling the ISTL
  //! vector interface.
  typedef X domain_type;
  //! \brief The range type of the preconditioner.
  //!
  //! E.g. BlockVector or another type fulfilling the ISTL
  //! vector interface.
  typedef Y range_type;
  //! \brief The field type of the preconditioner.
  typedef typename X::field_type field_type;
  //! \brief The type of the communication object..
  //!
  //! This must either be OwnerOverlapCopyCommunication or a type
  //! implementing the same interface.
  typedef C communication_type;

  // define the category
  enum {
    //! \brief The category the preconditioner is part of.
    category=Dune::SolverCategory::overlapping
  };

  /*! \brief Constructor.

    constructor gets all parameters to operate the prec.
    \param p The sequential preconditioner.
    \param c The communication object for syncing overlap and copy
    data points. (E.~g. OwnerOverlapCopyCommunication )
    */
  pASM(T& p, const communication_type& c)
    : preconditioner(p), communication(c)
  {   }

  /*!
    \brief Prepare the preconditioner.

    \copydoc Preconditioner::pre(X&,Y&)
    */
  virtual void pre (X& x, Y& b)
  {
//    communication.copyOwnerToAll(x,x);     // make dirichlet values consistent
    preconditioner.pre(x,b);
  }

  /*!
    \brief Apply the preconditioner

    \copydoc Preconditioner::apply(X&,const Y&)
    */
  virtual void apply (X& v, const Y& d)
  {
    preconditioner.apply(v,d);
    communication.addOwnerOverlapToAll(v,v);
  }

  template<bool forward>
    void apply (X& v, const Y& d)
    {
      preconditioner.template apply<forward>(v,d);
      communication.addOwnerOverlapToAll(v,v);
    }

  /*!
    \brief Clean up.

    \copydoc Preconditioner::post(X&)
    */
  virtual void post (X& x)
  {
    preconditioner.post(x);
  }

  private:
  //! \brief a sequential preconditioner
  T& preconditioner;

  //! \brief the communication object
  const communication_type& communication;
};
#endif

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

  const LinSolParams& solParams; //!< Reference to linear solver parameters.
  const ProcessAdm& adm;      //!< Reference to process administrator.
};

#endif
