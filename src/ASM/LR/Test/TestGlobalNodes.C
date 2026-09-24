//==============================================================================
//!
//! \file TestGlobalNodes.C
//!
//! \date Jun 21 2024
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for simple global node establishment for unstructured FE models.
//!
//==============================================================================

#include "LR/GlobalNodes.h"
#include "ASMuCube.h"
#include "ASMuSquare.h"

#include "LRSpline/LRSplineVolume.h"

#include "Catch2Support.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <utility>
#include <vector>


TEST_CASE("TestGlobalNodes.2D")
{
  ASMuSquare pch1;
  pch1.generateFEMTopology();
  ASMuSquare pch2(2, 1.0, 0.0);
  pch2.generateFEMTopology();
  ASMuSquare pch3(2, 1.0, 1.0);
  pch3.generateFEMTopology();

  std::vector<const LR::LRSpline*> splines{pch1.getBasis(),
                                           pch2.getBasis(),
                                           pch3.getBasis()};
  std::vector<ASM::Interface> ifs;
  ifs.push_back(ASM::Interface{1, 2, 2, 1, 0, 1, 1, 0});
  ifs.push_back(ASM::Interface{2, 3, 4, 3, 0, 1, 1, 0});

  auto nodes = GlobalNodes::calcGlobalNodes(splines, ifs);
  const auto ref = std::vector{
    GlobalNodes::IntVec{0, 1, 2, 3},
    GlobalNodes::IntVec{1, 4, 3, 5},
    GlobalNodes::IntVec{3, 5, 6, 7},
  };

  REQUIRE(nodes == ref);
}


TEST_CASE("TestGlobalNodes.3D")
{
  ASMuCube pch1;
  pch1.generateFEMTopology();
  ASMuCube pch2(2, 1.0, 0.0, 0.0);
  pch2.generateFEMTopology();
  ASMuCube pch3(2, 1.0, 1.0, 1.0);
  pch3.generateFEMTopology();

  std::vector<const LR::LRSpline*> splines{pch1.getBasis(),
                                           pch2.getBasis(),
                                           pch3.getBasis()};
  std::vector<ASM::Interface> ifs;
  ifs.push_back(ASM::Interface{1, 2, 2, 1, 0, 2, 1, 0});
  ifs.push_back(ASM::Interface{2, 3, 6, 5, 0, 2, 1, 0});

  auto nodes = GlobalNodes::calcGlobalNodes(splines, ifs);

  const auto ref = std::vector{
    GlobalNodes::IntVec{0,  1, 2,  3,  4,  5,  6,  7},
    GlobalNodes::IntVec{1,  8, 3,  9,  5, 10,  7, 11},
    GlobalNodes::IntVec{5, 10, 7, 11, 12, 13, 14, 15},
  };

  REQUIRE(nodes == ref);
}


TEST_CASE("TestGlobalNodes.Edges3D")
{
  ASMuCube pch;
  pch.generateFEMTopology();
  const LR::LRSpline& lr = *pch.getBasis();

  // The local edges of a tri-variate patch are numbered such that edges 1-4
  // run along u, 5-8 along v and 9-12 along w, see ASMs3D::getBoundary1Nodes.
  // Each of them is therefore the intersection of the two faces that meet
  // there. This pins the edge numbering of GlobalNodes::getBoundaryNodes,
  // which has no other test coverage and which multi-patch refinement relies
  // on to identify the functions on a shared boundary.
  const std::array<std::pair<int,int>,12> facesOfEdge = {{
    {3,5}, {4,5}, {3,6}, {4,6},   // along u: south/north x bottom/top
    {1,5}, {2,5}, {1,6}, {2,6},   // along v: west/east x bottom/top
    {1,3}, {2,3}, {1,4}, {2,4}    // along w: west/east x south/north
  }};

  std::vector<GlobalNodes::IntVec> edges(12);
  for (int lidx = 1; lidx <= 12; lidx++)
  {
    GlobalNodes::IntVec a = GlobalNodes::getBoundaryNodes(lr, 2, facesOfEdge[lidx-1].first, 0);
    GlobalNodes::IntVec b = GlobalNodes::getBoundaryNodes(lr, 2, facesOfEdge[lidx-1].second, 0);
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    GlobalNodes::IntVec expected;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::back_inserter(expected));
    REQUIRE(!expected.empty());

    edges[lidx-1] = GlobalNodes::getBoundaryNodes(lr, 1, lidx, 0);
    GlobalNodes::IntVec sorted = edges[lidx-1];
    std::sort(sorted.begin(), sorted.end());
    REQUIRE(sorted == expected);
  }

  // No two edges may consist of the same functions
  for (int i = 0; i < 12; i++)
    for (int j = i+1; j < 12; j++)
      REQUIRE(edges[i] != edges[j]);
}
