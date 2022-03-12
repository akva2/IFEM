// $Id$
//==============================================================================
//!
//! \file SplineFields2Dmx.h
//!
//! \date Mar 28 2011
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Class for mixed spline-based finite element vector fields in 2D.
//!
//==============================================================================

#ifndef _SPLINE_FIELDS_2D_MX_H
#define _SPLINE_FIELDS_2D_MX_H

#include "Fields.h"
#include <set>

class ASMs2Dmx;

namespace Go {
  class SplineSurface;
}


/*!
  \brief Class for mixed spline-based finite element vector fields in 2D.

  \details This class implements the methods required to evaluate a 2D mixed
  spline vector field at a given point in parametrical or physical coordinates.
*/

class SplineFields2Dmx : public Fields
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] basis Bases to use from patch
  //! \param[in] name Name of spline field
  SplineFields2Dmx(const ASMs2Dmx* patch, const RealArray& v,
                   char basis = 12, const char* name = nullptr);
  //! \brief Empty destructor.
  virtual ~SplineFields2Dmx() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  bool valueNode(size_t node, Vector& vals) const;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[out] vals Values in given physical coordinate
  virtual bool valueCoor(const Vec4& x, Vector& vals) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] vals Values in local point in given element
  bool valueFE(const ItgPoint& x, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const ItgPoint& x, Matrix& grad) const;

protected:
  const ASMs2Dmx* surf; //!< Patch description
  std::set<int>  bases; //!< Bases to use
};

#endif
