// $Id$
//==============================================================================
//!
//! \file PiolaMapping.h
//!
//! \date Apr 13 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Utilities for Piola mapping transformations.
//!
//==============================================================================

#ifndef _PIOLA_MAPPING_H
#define _PIOLA_MAPPING_H

#include "matrixnd.h"

#include <vector>

class FiniteElement;
class Vec3;


namespace utl
{
  void piolaMapping(FiniteElement& fe,
                    const double detJ,
                    const matrix<Real>& J,
                    const matrix<Real>& Ji,
                    const std::vector<matrix<Real>>& dNdu,
                    const matrix3d<Real>& H);

  void piolaBasis(FiniteElement& fe,
                  const utl::matrix<Real>& J);

  void piolaGradient (FiniteElement& fe,
                      const Real detJ,
                      const matrix<Real>& J,
                      const matrix<Real>& Ji,
                      const std::vector<matrix<Real>>& dNdu,
                      const matrix3d<Real>& H);

  std::vector<matrix<Real>> jacobianGradient(const matrix<Real>& dudX,
                                             const matrix3d<Real>& d2Xdu2);

  vector<Real> determinantGradient(const matrix<Real>& J,
                                   const matrix<Real>& Ji,
                                   const matrix3d<Real>& H);
}

#endif
