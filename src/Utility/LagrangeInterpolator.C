// $Id$
//==============================================================================
//!
//! \file LagrangeInterpolator.C
//!
//! \date Oct 18 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Lagrangian interpolation utilities.
//!
//==============================================================================

#include "LagrangeInterpolator.h"
#include <cstddef>


static double Lagrange(double x, size_t nr,
                       const std::vector<double>& grid)
{
  double result = 1.0;
  for (size_t j = 0; j< grid.size(); ++j) {
    if (j != nr)
      result *= (x-grid[j])/(grid[nr]-grid[j]);
  }

  return result;
}


static double LagrangeDiv1(double x, size_t nr,
                           const std::vector<double>& grid)
{
  double sum = 0.0;
  for (size_t i = 0; i < grid.size(); ++i)
    if (i != nr)
      sum += 1.0 / (x - grid[i]);
  return sum * Lagrange(x, nr, grid);
//  double result = 0.0;
//  double div = 1.0;
//  for (size_t i = 0; i < grid.size(); ++i) {
//    double prod = 1.0;
//    for (size_t j = 0; j < grid.size(); ++j)
//      if (j != i)
//        prod *= (x - grid[j]);
//    result += prod;
//    if (i != nr)
//      div *= grid[nr] - grid[i];
//  }

//  return result / div;
}


double LagrangeInterpolator::evaluate(double x,
                                      const std::vector<double>& data)
{
  double result = 0.0;
  for (size_t i = 0; i < data.size(); ++i)
    result += data[i]*Lagrange(x, i, grid);

  return result;
}


Matrix LagrangeInterpolator::get(const std::vector<double>& new_grid)
{
  Matrix result(grid.size(), new_grid.size());
  for (size_t i = 0; i < grid.size(); ++i)
    for (size_t j = 0; j < new_grid.size(); ++j)
      result(i+1, j+1) =  Lagrange(new_grid[j], i, grid);

  return result;
}


double LagrangeHermiteInterpolator::evaluate(double x,
                                             const std::vector<double>& values,
                                             const std::vector<double>& tangents)
{
  double result = 0.0;
  for (size_t i = 0; i < grid.size(); ++i) {
    double l = Lagrange(x, i, grid);
    double l2 = l*l;
    result += (1.0 - 2.0*(x - grid[i])*LagrangeDiv1(x,i,grid))*l2*values[i] +
              (x - grid[i])*l2*tangents[i];
  }

  return result;
}
