//==============================================================================
//!
//! \file TestMultiPatchLRRefine.C
//!
//! \date Mar 31 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for multi-patch LR refinement.
//!
//==============================================================================

#include "ASMbase.h"
#include "ASMunstruct.h"
#include "IntegrandBase.h"
#include "SIM2D.h"
#include "SIM3D.h"

#include "gtest/gtest.h"

#include <fstream>


// Dummy SIM class.
template<class Dim>
class RefineSim : public Dim
{
public:
  RefineSim() : Dim(Dim::dimension)
  { 
    Dim::opt.discretization = ASM::LRSpline;
  }
  virtual ~RefineSim() {}
};


class TestMultiPatchLRRefine2D: public testing::Test,
                                public testing::WithParamInterface<int>
{
};


TEST_P(TestMultiPatchLRRefine2D, Refine)
{
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam() << ".xinp";

  RefineSim<SIM2D> sim;
  EXPECT_TRUE(sim.read(str.str().c_str()));
  EXPECT_TRUE(sim.preprocess());

  for (size_t i = 0; i < 3; ++i) {
    LR::RefineData prm;
    sim.getPatch(2)->getBoundaryNodes(1, prm.elements);

    prm.options.resize(3);
    prm.options[0] = 1;
    prm.options[1] = 1;
    prm.options[2] = 2;
    sim.refine(prm);
    sim.clearProperties();

    bool result = sim.read(str.str().c_str());
    if (!result) {
      std::ofstream os("failure.lr");
      sim.dumpGeometry(os);
    }

    EXPECT_TRUE(result);
  }
}


class TestMultiPatchLRRefine3D: public testing::Test,
                                public testing::WithParamInterface<int>
{
};


TEST_P(TestMultiPatchLRRefine3D, Refine)
{
  std::stringstream str;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam() << ".xinp";

  RefineSim<SIM3D> sim;
  EXPECT_TRUE(sim.read(str.str().c_str()));
  EXPECT_TRUE(sim.preprocess());

  for (size_t i = 0; i < 3; ++i) {
    LR::RefineData prm;
    sim.getPatch(2)->getBoundaryNodes(1, prm.elements);

    prm.options.resize(3);
    prm.options[0] = 1;
    prm.options[1] = 1;
    prm.options[2] = 2;
    sim.refine(prm);
    sim.clearProperties();

    bool result = sim.read(str.str().c_str());
    if (!result) {
      std::ofstream os("failure.lr");
      sim.dumpGeometry(os);
    }
    EXPECT_TRUE(result);
  }
}


const std::vector<int> orientations2D = {0,1};
INSTANTIATE_TEST_CASE_P(TestMultiPatchLRRefine2D, TestMultiPatchLRRefine2D, testing::ValuesIn(orientations2D));

const std::vector<int> orientations3D = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
INSTANTIATE_TEST_CASE_P(TestMultiPatchLRRefine3D, TestMultiPatchLRRefine3D, testing::ValuesIn(orientations3D));
