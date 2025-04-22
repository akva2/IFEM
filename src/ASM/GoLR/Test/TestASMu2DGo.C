//==============================================================================
//!
//! \file TestASMu2DGo.C
//!
//! \date Apr 8 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for assembly of unstructured 2D Go spline FE models.
//!
//==============================================================================

#include "ASMuSquareGo.h"
#include "GaussQuadrature.h"
#include "SIM2D.h"

#include "gtest/gtest.h"
#include <fstream>
#include <numeric>


struct EdgeTest
{
  int edge;
  int edgeIdx;
  std::array<int,2> c1;
  std::array<int,2> c2;
};


class TestASMu2D : public testing::Test,
                   public testing::WithParamInterface<EdgeTest>
{
};


class TestuSIM2D : public SIM2D
{
public:
  TestuSIM2D() : SIM2D(1)
  {
    opt.discretization = ASM::LRSplineGo;
    EXPECT_TRUE(this->read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
    EXPECT_TRUE(this->createFEMmodel());
  }
  virtual ~TestuSIM2D() {}
};


TEST_P(TestASMu2DGo, ConstrainEdge)
{
  TestuSIM2D sim;
  ASMu2DGo* pch = static_cast<ASMu2DGo*>(sim.getPatch(1));
  ASSERT_TRUE(pch != nullptr);

  pch->constrainEdge(GetParam().edgeIdx, false, 1, 1, 1);

  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1);

  for (int node : glbNodes)
    EXPECT_TRUE(pch->findMPC(node,1) != nullptr);
}
