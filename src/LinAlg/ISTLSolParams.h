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
#include "LinSolParams.h"

#include <iostream>
#include <set>
#include <string>
#include <vector>


class ProcessAdm;
class SettingMap;



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
  //! \param op Matrix operator to use
  //! \return tuple with (solver, preconditioner)
  std::tuple<std::unique_ptr<ISTL::InverseOperator>, std::unique_ptr<ISTL::Preconditioner>>
    setupPC(ISTL::Mat& A, ISTL::Operator& op);

  //! \brief Obtain linear solver parameters.
  const LinSolParams& get() const { return solParams; }

protected:
  const LinSolParams& solParams; //!< Reference to linear solver parameters.
  const ProcessAdm& adm;      //!< Reference to process administrator.
};

#endif
