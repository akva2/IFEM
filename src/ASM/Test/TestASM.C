//==============================================================================
//!
//! \file TestASM.C
//!
//! \date Mar 29 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for ASM modules
//!
//==============================================================================

#include "SAM.h"
#include "SIM2D.h"
#include "ASMbase.h"

#include "gtest/gtest.h"

#include <fstream>


static void check_intmatrices_equal (const IntVec& first, const IntVec& second)
{
  ASSERT_TRUE(first.size() == second.size());
  auto it = second.begin();
  for (const auto& it2 : first) {
    ASSERT_EQ(it2, *it);
    ++it;
  }
};


TEST(TestASM, ASMs2DBoundaryNodes)
{
  SIM2D sim(1);
  sim.read("src/ASM/Test/refdata/asm_2D_2P.xinp");
  sim.preprocess();

  IntVec genNodes;
  sim.getPatch(1)->getBoundaryNodes(2, genNodes);
  IntVec expNodes(5);
  int n = 0;
  std::generate(expNodes.begin(), expNodes.end(), [&n]{return n += 5; });
  check_intmatrices_equal(genNodes, expNodes);
  genNodes.clear();
  sim.getPatch(2)->getBoundaryNodes(2, genNodes);
  n = 25;
  std::generate(expNodes.begin(), expNodes.end(), [&n]{return n += 4; });
  check_intmatrices_equal(genNodes, expNodes);
}


TEST(TestASM, ASMs2DmxBoundaryNodes)
{
  SIM2D sim({1,1});
  sim.read("src/ASM/Test/refdata/asm_2D_2P.xinp");
  sim.preprocess();

  IntVec genNodes;
  sim.getPatch(1)->getBoundaryNodes(2, genNodes, 1);
  IntVec expNodes(6);
  int n = 0;
  std::generate(expNodes.begin(), expNodes.end(), [&n]{return n += 6; });
  check_intmatrices_equal(genNodes, expNodes);

  genNodes.clear();
  sim.getPatch(1)->getBoundaryNodes(2, genNodes, 2);
  n = 41-5;
  IntVec expNodes2(5);
  std::generate(expNodes2.begin(), expNodes2.end(), [&n]{return n += 5; });
  check_intmatrices_equal(genNodes, expNodes2);
  expNodes.insert(expNodes.end(), expNodes2.begin(), expNodes2.end());
  genNodes.clear();
  sim.getPatch(1)->getBoundaryNodes(2, genNodes, 0);
  check_intmatrices_equal(genNodes, expNodes);

  genNodes.clear();
  expNodes.resize(6);
  sim.getPatch(2)->getBoundaryNodes(2, genNodes, 1);
  n = 66-5;
  std::generate(expNodes.begin(), expNodes.end(), [&n]{return n += 5; });
  check_intmatrices_equal(genNodes, expNodes);

  genNodes.clear();
  sim.getPatch(2)->getBoundaryNodes(2, genNodes, 2);
  n = 95-4;
  std::generate(expNodes2.begin(), expNodes2.end(), [&n]{return n += 4; });
  check_intmatrices_equal(genNodes, expNodes2);
  expNodes.insert(expNodes.end(), expNodes2.begin(), expNodes2.end());

  genNodes.clear();
  sim.getPatch(2)->getBoundaryNodes(2, genNodes, 0);
  check_intmatrices_equal(genNodes, expNodes);
}
