//==============================================================================
//!
//! \file TestASMs2D.C
//!
//! \date Feb 14 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 2D spline FE models.
//!
//==============================================================================

#include "ASMs2D.h"
#include "SIM2D.h"

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "gtest/gtest.h"
#include <fstream>
#include <numeric>


class TestASMs2D : public testing::Test,
                   public testing::WithParamInterface<int>
{
};


TEST_P(TestASMs2D, Connect)
{
  SIM2D sim(1);
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.createFEMmodel());
}

const std::vector<int> orientations2D = {0,1};
INSTANTIATE_TEST_CASE_P(TestASMs2D, TestASMs2D, testing::ValuesIn(orientations2D));


TEST(TestASMs2D, ExtractBasis)
{
  std::ifstream iff("srf.g2");
  Go::ObjectHeader head;
  Go::SplineSurface srf;
  iff >> head >> srf;

  auto&& value = [](double x, double y)
    {
      return pow(x,3.0)*pow(y-1,3.0)/(2*((pow(x,3.0)*pow(y-1.0,3.0)/2 - pow(x-1,3.0)*pow(y-1,3.0))));
    };

  auto&& derx = [](double x, double y)
    {
      return 6*pow(x,2.0)*pow(x-1,2.0)/pow(pow(x,3.0)-6*pow(x,2.0)+6*x-2,2.0);
    };

  auto&& derxx = [](double x, double y)
    {
      return -(12*x*(pow(x,5.0)-3*pow(x,4.0)+2*pow(x,3.0)+4*pow(x,2)-6*x+2))/pow(pow(x,3.0)-6*pow(x,2.0)+6*x-2,3.0);
    };

  auto&& derxxx = [](double x, double y)
    {
      return 12.0*(3.0*pow(x,8.0) - 12.0*pow(x,7.0) + 10.0*pow(x,6.0) + 48.0*pow(x,5)-156.0*pow(x,4.0)+176.0*pow(x,3.0)-72.0*pow(x,2.0) + 4.0)/pow(pow(x,3.0)-6*pow(x,2.0)+6*x-2,4.0);
    };

  auto&& trunc = [](double x)
  {
    return fabs(x) < 1e-11 ? 0.0 : x;
  };

  for (double x : {0.1, 0.9}) {
    Go::BasisDerivsSf3 derivs;
    srf.computeBasis(x, x, derivs);

    std::cout << "Value is: " << value(x,x) << std::endl;
    std::cout << "Derivative is: " << derx(x,x) << std::endl;
    std::cout << "DDerivative is: " << derxx(x,x) << std::endl;
    std::cout << "DDDerivative is: " << derxxx(x,x) << std::endl;
    std::vector<Go::Point> pt(10);
    srf.point(pt, x, x, 3);
    std::cout << "N: ";
    for (auto& it : derivs.basisValues)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "dN/du: ";
    for (auto& it : derivs.basisDerivs_u)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "dN/dv: ";
    for (auto& it : derivs.basisDerivs_v)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "d2N/du2: ";
    for (auto& it : derivs.basisDerivs_uu)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "d2N/dv2: ";
    for (auto& it : derivs.basisDerivs_vv)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "d2N/dudv: ";
    for (auto& it : derivs.basisDerivs_uv)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "d3N/du3: ";
    for (auto& it : derivs.basisDerivs_uuu)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "d3N/dv3: ";
    for (auto& it : derivs.basisDerivs_vvv)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "d3N/du2dv: ";
    for (auto& it : derivs.basisDerivs_uuv)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    std::cout << "d3N/dudv2: ";
    for (auto& it : derivs.basisDerivs_uvv)
      std::cout << trunc(it) << " ";
    std::cout << std::endl;
    const std::array<std::string,10> lab = {"value", "d/du", "d/dv", "d2/du2", "d2/dudv",
                                            "d2/dv2", "d3/du3", "d3/du2dv", 
                                             "d3/dudv2", "d3/dv3"};
    auto it = lab.begin();
    for (auto& p : pt) {
      std::cout << "Surface " << *it << ": " 
                << trunc(p[0]) << " " << trunc(p[1]) << std::endl;
      ++it;
    }
  }
}
