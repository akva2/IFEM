//==============================================================================
//!
//! \file TestASMs1D.C
//!
//! \date Nov 9 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 1D spline FE models.
//!
//==============================================================================


#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "gtest/gtest.h"

#include <fstream>
#include <numeric>


TEST(TestASMs1D, ExtractBasis)
{
  std::stringstream str;
  str.str("100 1 0 0\n"
          "2 1\n"
          "4 4\n"
          "0 0 0 0 1 1 1 1\n"
          "0 0 1\n"
          "0 0 0\n"
          "0 0 0\n"
          "0.5 0 0.5\n");

  Go::ObjectHeader head;
  Go::SplineCurve crv;
  str >> head >> crv;

  auto&& value = [](double x)
    {
      return 0.5*pow(x,3.0)/(pow(1.0-x,3.0)+0.5*pow(x,3.0));
    };

  auto&& derivative = [](double x)
      { 
        return 6.0*pow(1-x,2.0)*pow(x,2.0)/pow(-pow(x,3.0)+6.0*pow(x,2)-6*x+2,2.0);
      };

  auto&& dderivative = [](double x)
      { 
        return -x*(-12.0*pow(x,5)+36*pow(x,4)-24*pow(x,3)-48*pow(x,2)+72*x-24)/pow(-pow(x,3.0)+6*pow(x,2)-6*x+2,3.0);
      };

  auto&& ddderivative=[](double x)
      {
        return 12.0*(3.0*pow(x,8.0) - 12.0*pow(x,7.0) + 10.0*pow(x,6.0) + 48.0*pow(x,5)-156.0*pow(x,4.0)+176.0*pow(x,3.0)-72.0*pow(x,2.0) + 4.0)/pow(pow(x,3.0)-6*pow(x,2.0)+6*x-2,4.0);
      };

  std::vector<double> N, dNdu, d2Ndu2, d3Ndu3;
  for (double x : {0.1, 0.9}) {
    crv.computeBasis(x, N, dNdu, d2Ndu2, d3Ndu3);

    std::cout << "Value is: " << value(x) << std::endl;
    std::cout << "Derivative is: " << derivative(x) << std::endl;
    std::cout << "DDerivative is: " << dderivative(x) << std::endl;
    std::cout << "DDDerivative is: " << ddderivative(x) << std::endl;
    std::vector<Go::Point> pt(4);
    crv.point(pt, x, 3);
    std::cout << "Basis functions are: ";
    for (auto& it : N)
      std::cout << it << " ";
    std::cout << std::endl;
    std::cout << "Basis derivatives are: ";
    for (auto& it : dNdu)
      std::cout << it << " ";
    std::cout << std::endl;
    std::cout << "Basis second derivatives are: ";
    for (auto& it : d2Ndu2)
      std::cout << it << " ";
    std::cout << std::endl;
    std::cout << "Basis third derivatives are: ";
    for (auto& it : d3Ndu3)
      std::cout << it << " ";
    std::cout << std::endl;
    std::cout << "Curve 0der value: " << pt[0] << std::endl;
    std::cout << "Curve 1der value: " << pt[1] << std::endl;
    std::cout << "Curve 2der value: " << pt[2] << std::endl;
    std::cout << "Curve 3der value: " << pt[3] << std::endl;
  }
}
