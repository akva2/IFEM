// $Id$
//==============================================================================
//!
//! \file ISTLSupport.h
//!
//! \date May 12 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief IFEM ISTL support
//!
//==============================================================================

#ifndef _ISTL_SUPPORT_H_
#define _ISTL_SUPPORT_H_

#include <vector>

#ifdef HAS_ISTL
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/schwarz.hh>
#ifdef HAVE_SUPERLU
#include <dune/istl/superlu.hh>
#endif
#ifdef HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif

namespace ISTL
{
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat; //!< A sparse system matrix
  typedef Dune::BlockVector<Dune::FieldVector<double,1>> Vec;  //!< A vector
//  typedef Dune::AssembledLinearOperator<Mat,Vec,Vec> Operator; //!< Linear operator abstraction
  typedef Dune::MatrixAdapter<Mat,Vec,Vec> Operator;      //!< A serial matrix operator
  typedef Dune::InverseOperator<Vec, Vec> InverseOperator;     //!< Linear system inversion abstraction
  typedef Dune::Preconditioner<Vec,Vec> Preconditioner;        //!< Preconditioner abstraction
  typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Dune::OwnerOverlapCopyCommunication<int,int>> ParMatrixAdapter; //!< A parallel matrix operator

  /*! \brief Wrapper template to avoid memory leaks */

  template<template<class M> class Pre>
  class IOp2Pre : public Dune::InverseOperator2Preconditioner<Pre<ISTL::Mat>, Dune::SolverCategory::sequential> {
    typedef Dune::InverseOperator2Preconditioner<Pre<ISTL::Mat>, Dune::SolverCategory::sequential> SolverType;
  public:
    IOp2Pre(Pre<ISTL::Mat>* iop) : SolverType(*iop)
    {
      m_op.reset(iop);
    }
  protected:
    std::unique_ptr<Pre<ISTL::Mat>> m_op;
  };

#if defined(HAVE_UMFPACK)
  typedef Dune::UMFPack<ISTL::Mat> LUType;
  typedef IOp2Pre<Dune::UMFPack> LU;
#elif defined(HAVE_SPUERLU)
  typedef Dune::SuperLUk<ISTL::Mat> LUType;
  typedef IOp2Pre<Dune::SuperLU> LU;
#endif

}
#endif

#endif
