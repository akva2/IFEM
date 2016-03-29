//==============================================================================
//!
//! \file TestSAM.C
//!
//! \date Mar 29 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for SAM 
//!
//==============================================================================

#include "ASMbase.h"
#include "SAM.h"
#include "SIM2D.h"
#include "tinyxml.h"

#include "gtest/gtest.h"

#include <fstream>


typedef std::vector< std::vector<int> > IntMat;


static IntMat readIntMatrix(size_t r, const std::string& file)
{
  std::vector< std::vector<int> > result;
  result.resize(r);
  std::ifstream f(file);
  for (size_t i=0;i<r;++i) {
    size_t size;
    f >> size;
    result[i].resize(size);
    for (size_t j=0;j<size;++j)
      f >> result[i][j];
  }

  return result;
}


auto&& check_intmatrices_equal = [](const std::vector<std::set<int>>& A,
                                    const std::string& path)
{
  IntMat B = readIntMatrix(A.size(), path);
  size_t i = 0;
  for (const auto& it : A) {
    size_t j = 0;
    for (const auto& it2 : it)
      ASSERT_EQ(it2, B[i][j++]);
    ++i;
  }
};


TEST(TestSAM, SingleBasis1P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_1P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);

  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_singlebasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, SingleBasisDirichlet1P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);

  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_1P.ref");

  int eq = 1;
  for (int i = 1; i <= 3; ++i) {
    ASSERT_EQ(sam->getEquation(3*(i-1)+1, 1), eq++);
    ASSERT_EQ(sam->getEquation(3*(i-1)+2, 1), eq++);
    ASSERT_EQ(sam->getEquation(3*(i-1)+3, 1), 0);
  }
}


TEST(TestSAM, MixedBasis1P)
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_1P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);
  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_1P.ref");

  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, MixedBasisDirichlet1P)
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_1P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);
  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_1P.ref");

  int eq = 1;
  for (int i = 1; i <= 4; ++i) {
    ASSERT_EQ(sam->getEquation(4*(i-1)+1, 1), eq++);
    ASSERT_EQ(sam->getEquation(4*(i-1)+2, 1), eq++);
    ASSERT_EQ(sam->getEquation(4*(i-1)+3, 1), eq++);
    ASSERT_EQ(sam->getEquation(4*(i-1)+4, 1), 0);
  }

  for (int i = 17; i <= 25; ++i)
    ASSERT_EQ(sam->getEquation(i, 1), eq++);
}


TEST(TestSAM, SingleBasis2P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_2P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);
  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_singlebasis_2P.ref");
  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, MixedBasis2P)
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_2P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);
  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_2P.ref");
  for (int i = 1; i <= sam->getNoNodes(); ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);
}


TEST(TestSAM, SingleBasisDirichlet2P)
{
  SIM2D sim(1);
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);
  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_singlebasis_dir_2P.ref");
  ASSERT_EQ(sam->getEquation(1, 1), 1);
  ASSERT_EQ(sam->getEquation(2, 1), 2);
  ASSERT_EQ(sam->getEquation(3, 1), 3);
  ASSERT_EQ(sam->getEquation(4, 1), 4);
  ASSERT_EQ(sam->getEquation(5, 1), 0);
  ASSERT_EQ(sam->getEquation(6, 1), 0);
}


TEST(TestSAM, MixedBasisDirichlet2P)
{
  SIM2D sim({1,1});
  sim.read("src/LinAlg/Test/refdata/sam_2D_dir_2P.xinp");
  sim.preprocess();

  std::vector<std::set<int>> dofc;
  const SAM* sam = sim.getSAM();
  sam->getDofCouplings(dofc);

  check_intmatrices_equal(dofc, "src/LinAlg/Test/refdata/sam_2D_mixedbasis_dir_2P.ref");
  for (int i = 1; i <= 13; ++i)
    ASSERT_EQ(sam->getEquation(i, 1), i);

  const ASMbase* pch = sim.getPatch(2);
  int eq = 14;
  for (int i = 1; i <= 3; ++i) {
    ASSERT_EQ(sam->getEquation(pch->getNodeID((i-1)*3+2), 1), eq++); 
    ASSERT_EQ(sam->getEquation(pch->getNodeID((i-1)*3+3), 1), 0); 
  }
  ASSERT_EQ(sam->getEquation(20, 1), eq++);
  ASSERT_EQ(sam->getEquation(21, 1), eq++);
}
