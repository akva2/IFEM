// $Id$
//==============================================================================
//!
//! \file ASMu2Dnurbs.h
//!
//! \date Nov 28 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D NURBS FE models.
//!
//==============================================================================

#ifndef _ASM_U2D_NURBS_H
#define _ASM_U2D_NURBS_H

#include "ASMu2D.h"


/*!
  \brief Driver for assembly of unstructured 2D NURBS FE models.
  \details This class contains methods common for 2D LR-NURBS patches.
*/

class ASMu2Dnurbs : public ASMu2D
{
public:
  //! \brief Default constructor.
  ASMu2Dnurbs(unsigned char n_s = 2, unsigned char n_f = 2);
  //! \brief Copy constructor.
  ASMu2Dnurbs(const ASMu2Dnurbs& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMu2Dnurbs() {}

  // Methods for model generation and refinement
  // ===========================================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream&);

  //! \brief Refines the parametrization by inserting tensor knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int dir, int nInsert);
  //! \brief Refines the parametrization by inserting extra tensor knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  //! \param[in] scale Scaling factor for the added knot values
  virtual bool refine(int dir, const RealArray& xi, double scale);
  //! \brief Raises the order of the tensor spline object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  virtual bool raiseOrder(int ru, int rv);

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel 1-based element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const;

protected:
  //! \brief Evaluates the basis functions and derivatives of order \a derivs
  //! of an element.
  virtual bool evaluateBasis(int iel, FiniteElement& el, int derivs = 0) const;

  //! \brief Evaluate basis functions in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisPtsSf& bas, int iel,
                            const LR::LRSplineSurface* spline = nullptr) const;

  //! \brief Evaluate basis functions and first derivatives in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisDerivsSf& bas, int iel,
                            const LR::LRSplineSurface* spline = nullptr) const;
  //! \brief Evaluate basis functions and two derivatives in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisDerivsSf2& bas, int iel) const;

  //! \brief Evaluate basis functions and two derivatives in a point.
  virtual void computeBasis(double u, double v,
                            Go::BasisDerivsSf3& bas, int iel) const;

  //! \brief Converts current tensor spline object to LR-spline.
  virtual LR::LRSplineSurface* createLRfromTensor();

private:
  bool noNurbs; //!< If true, we read a spline and thus forward to ASMu2D.
};

#endif
