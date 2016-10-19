// $Id$
//==============================================================================
//!
//! \file LagrangeInterpolator.h
//!
//! \date Oct 18 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Lagrangian interpolation utilities.
//!
//==============================================================================
#ifndef LAGRANGE_INTERPOLATOR_H_
#define LAGRANGE_INTERPOLATOR_H_

#include "MatVec.h"
#include <vector>

/*! \brief Lagrangian interpolation class.
*/

class LagrangeInterpolator {
public:
  //! \brief Constructor.
  //! \param grid The associated abscissa
  LagrangeInterpolator(const std::vector<double>& grid_) : grid(grid_) {}

  //! \brief Evaluate interpolation polynomial at coordinate.
  //! \param x Coordinate to evaluate interpolant at
  //! \param data Function values at interpolant abscissa
  double evaluate(double x, const std::vector<double>& data);

  //! \brief Get an interpolation matrix from the associated grid to a new grid.
  //! \param new_grid Coordinates on the new grid
  Matrix get(const std::vector<double>& new_grid);

  protected:
    std::vector<double> grid; //!< Grid points.
};


/*! \brief Lagrangian Hermite interpolation class.
*/

class LagrangeHermiteInterpolator {
public:
  //! \brief Constructor.
  //! \param grid The associated abscissa
  LagrangeHermiteInterpolator(const std::vector<double>& grid_) : grid(grid_) {}

  //! \brief Evaluate interpolation polynomial at coordinate.
  //! \param x Coordinate to evaluate interpolant at
  //! \param values Function values at interpolant abscissa
  //! \param tangent Tangent values at interpolant abscissa
  double evaluate(double x, const std::vector<double>& values,
                  const std::vector<double>& tangents);

  protected:
    std::vector<double> grid; //!< Grid points.
};

#endif
