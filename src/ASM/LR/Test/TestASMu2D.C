//==============================================================================
//!
//! \file TestASMu2D.C
//!
//! \date Jul 14 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for assembly of unstructured 2D spline FE models.
//!
//==============================================================================

#include "ASMu2D.h"
#include "IntegrandBase.h"
#include "SIM2D.h"

#include "gtest/gtest.h"
#include "LRSpline/LRSplineSurface.h"

#include <fstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif


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


class DummyIntegrand : public IntegrandBase {};


static ASMu2D* getPatch (SIMinput& sim, const char* file)
{
  sim.opt.discretization = ASM::LRSpline;
  EXPECT_TRUE(sim.read(file));
  EXPECT_TRUE(sim.createFEMmodel());
  return static_cast<ASMu2D*>(sim.getPatch(1));
}


template<class Dim>
class AdaptiveTestSIM : public Dim
{
public:
  AdaptiveTestSIM(const std::string& input) :  Dim(1), file(input)
  { this->myProblem = new DummyIntegrand; }

  bool adaptMesh(const std::vector<int>& elems)
  {
    LR::RefineData prm;
    prm.options.resize(3);
    prm.options[2] = 2;
    prm.elements = elems;

    if (!this->refine(prm))
      return false;
    this->clearProperties();
    if (!this->read(file.c_str()))
      return false;
    return this->preprocess();
  }
protected:
  std::string file;
};


TEST_P(TestASMu2D, ConstrainEdge)
{
  SIM2D sim(1);
  ASMu2D* pch = getPatch(sim, "src/ASM/LR/Test/refdata/boundary_nodes.xinp");
  ASSERT_TRUE(pch != nullptr);
  pch->constrainEdge(GetParam().edgeIdx, false, 1, 1, 1);
  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1, 1);
  for (int& it : glbNodes)
    ASSERT_TRUE(pch->findMPC(it, 1) != nullptr);
}


TEST_P(TestASMu2D, ConstrainEdgeOpen)
{
  SIM2D sim(1);
  ASMu2D* pch = getPatch(sim, "src/ASM/LR/Test/refdata/boundary_nodes.xinp");
  ASSERT_TRUE(pch != nullptr);
  pch->constrainEdge(GetParam().edgeIdx, true, 1, 1, 1);
  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1, 1);
  int crn = pch->getCorner(GetParam().c1[0], GetParam().c1[1], 1);
  ASSERT_TRUE(pch->findMPC(crn, 1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  crn = pch->getCorner(GetParam().c2[0], GetParam().c2[1], 1);
  ASSERT_TRUE(pch->findMPC(crn, 1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  for (int& it : glbNodes)
    ASSERT_TRUE(pch->findMPC(it, 1) != nullptr);
}


static const std::vector<EdgeTest> edgeTestData =
        {{1, -1, {-1, -1}, {-1 , 1}},
         {2,  1, { 1, -1}, { 1 , 1}},
         {3, -2, {-1, -1}, { 1, -1}},
         {4,  2, {-1,  1}, { 1,  1}}};


INSTANTIATE_TEST_CASE_P(TestASMu2D, TestASMu2D,
                        testing::ValuesIn(edgeTestData));

TEST(TestASMu2D, ThreadGroups)
{
  AdaptiveTestSIM<SIM2D> sim("src/ASM/LR/Test/refdata/threadgroups.xinp");
#ifdef USE_OPENMP
  omp_set_num_threads(2);
#endif
  ASMu2D* pch = getPatch(sim, "src/ASM/LR/Test/refdata/threadgroups.xinp");
  const ThreadGroups& g = pch->getThreadGroups();
  sim.preprocess();
#ifdef USE_OPENMP
  static const std::vector<std::vector<std::vector<int>>> eq_groups =
   {{{0,1,2,3,4}, {10,11,12,13,14}},
    {{5,6,7,8,9}, {15,16,17,18,19,20,21,22,23,24}}};

  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j< 2; ++j) {
      for (size_t e = 0; e < g[i][j].size(); ++e)
        ASSERT_EQ(g[i][j][e], eq_groups[i][j][e]);
    }
  ASSERT_TRUE(sim.adaptMesh({0}));
  static const std::vector<std::vector<std::vector<int>>> ref1_groups =
    {{{0, 25, 26, 27, 1, 2, 3, 4}, {10, 11, 12, 13, 14}},
     {{5,6,7,8,9}, {15,16,17,18,19,20,21,22,23,24}}};

  const ThreadGroups& g2 = pch->getThreadGroups();
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j< 2; ++j) {
      for (size_t e = 0; e < g[i][j].size(); ++e)
        ASSERT_EQ(g2[i][j][e], ref1_groups[i][j][e]);
    }
#else
  ASSERT_EQ(g[0][0].size(), 16);
  ASSERT_EQ(g[0][1].size(), 0);
  ASSERT_EQ(g[1].size(), 0);
#endif
}
