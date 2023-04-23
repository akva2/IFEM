// $Id$
//==============================================================================
//!
//! \file LRNURBSFields2D.h
//!
//! \date Apr 23 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR NURBS-based finite element vector fields in 2D.
//!
//==============================================================================

#ifndef _LRNURBS_FIELDS_2D_H
#define _LRNURBS_FIELDS_2D_H

#include "Fields.h"

class ASMu2D;

namespace LR {
  class Element;
  class LRSplineSurface;
}


/*!
  \brief Class for LR NURBS-based finite element vector fields in 2D.

  \details This class implements the methods required to evaluate a 2D LR NURBS
  vector field at a given point in parametrical or physical coordinates.
*/

class LRNURBSFields2D : public Fields
{
public:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] patch The spline patch on which the field is to be defined
  //! \param[in] v Array of control point field values
  //! \param[in] basis Basis to use from patch
  //! \param[in] nf Number of components to for field
  //! \param[in] name Name of spline field
  LRNURBSFields2D(const ASMu2D* patch, const RealArray& v,
                  char basis = 1, int nf = 0, const char* name = nullptr);

  //! \brief Construct directly from surface.
  //! \param[in] srf The spline surface to use
  //! \param[in] v Array of control point field values
  //! \param[in] ncmp Number of field components
  //! \param[in] name Name of spline field
  LRNURBSFields2D(const LR::LRSplineSurface* srf, const RealArray& v, int ncmp,
                  const char* name = nullptr);

  //! \brief Empty destructor.
  virtual ~LRNURBSFields2D() {}

  // Methods to evaluate the field
  //==============================

  //! \brief Computes the value in a given node/control point.
  //! \param[in] node Node number
  //! \param[out] vals Node values
  bool valueNode(size_t node, Vector& vals) const;

  //! \brief Computes the value at a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] vals Values in local point in given element
  bool valueFE(const ItgPoint& x, Vector& vals) const;

  //! \brief Computes the value at a given global coordinate.
  //! \param[in] x Global/physical coordinate for point
  //! \param[out] vals Values in given physical coordinate
  bool valueCoor(const Vec4& x, Vector& vals) const;

  //! \brief Computes the gradient for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] grad Gradient of solution in a given local coordinate
  bool gradFE(const ItgPoint& x, Matrix& grad) const;

  //! \brief Computes the hessian for a given local coordinate.
  //! \param[in] x Local coordinate of evaluation point
  //! \param[out] H Hessian of solution in a given local coordinate
  virtual bool hessianFE(const ItgPoint& x, Matrix3D& H) const;

protected:
  const LR::LRSplineSurface* basis; //!< Spline basis description
  const LR::LRSplineSurface* surf;  //!< Spline geometry description
};

#endif
