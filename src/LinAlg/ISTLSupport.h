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
#include "dune/istl/bcrsmatrix.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/solvers.hh"

namespace ISTL
{
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat; //!< A sparse system matrix
  typedef Dune::BlockVector<Dune::FieldVector<double,1>> Vec;  //!< A vector
  typedef Dune::MatrixAdapter<Mat,Vec,Vec> Operator;           //!< Linear operator abstraction
  typedef Dune::InverseOperator<Vec, Vec> InverseOperator;     //!< Linear system inversion abstraction
  typedef Dune::Preconditioner<Vec,Vec> Preconditioner;        //!< Preconditioner abstraction
}
#endif

#endif
