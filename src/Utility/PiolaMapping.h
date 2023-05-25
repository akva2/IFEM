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
  bool piolaMapping(FiniteElement& fe,
                    const matrix<Real>& J,
                    const matrix<Real>& Ji,
                    const std::vector<matrix<Real>>& dNdu,
                    const matrix3d<Real>& H);

  void piolaBasis(FiniteElement& fe,
                  const utl::matrix<Real>& J);
}

#endif
