//==============================================================================
//!
//! \file TestASMs2Dmx.C
//!
//! \date Jun 24 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for mixed 2D assembly handler
//!
//==============================================================================

#include "ASMs2Dmx.h"
#include "SIM2D.h"
#include "IntegrandBase.h"
#include "AlgEqSystem.h"
#include "ElmMats.h"
#include "FiniteElement.h"
#include "gtest/gtest.h"

class TestNitscheIntegrand : public IntegrandBase {
public:
  TestNitscheIntegrand(unsigned short int n = 0) : IntegrandBase(n) {}

  bool hasInteriorTerms() const override { return false; }
  bool hasBoundaryTerms() const override { return false; }
  bool hasNitscheTerms() const override { return true; }

  LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                  size_t, bool) const override
  {
    ElmMats* result = new ElmMats;
    result->resize(1,1);
    size_t siz = 2*(nen[0]+nen[1]);
    result->A[0].resize(siz, siz);
    result->b[0].resize(siz);

    return result;
  }

  bool evalBouMx(LocalIntegral& elmInt,
                 const MxFiniteElement& fe,
                 const TimeDomain& time,
                  const Vec3& X, const Vec3& normal) const override
  {
    ElmMats& A = static_cast<ElmMats&>(elmInt);

    for (size_t i = 1; i <= fe.basis(1).size(); ++i) {
      A.b[0](i+fe.basis(1).size() + fe.basis(2).size()) = -1.0;
      A.b[0](i) = 1.0;
    }

    return true;
  }
};


class TestNitscheSIM : public SIM2D
{
public:
  TestNitscheSIM() : SIM2D({1,1})
  {
    myProblem = new TestNitscheIntegrand;
  }
};


TEST(TestASMs2Dmx, Nitsche)
{
  TestNitscheSIM sim;
  if (!sim.read("src/ASM/Test/refdata/ASMs2Dmx_Nitsche.xinp"))
    ASSERT_EQ(0, 1);

  if (!sim.preprocess())
    ASSERT_EQ(0, 2);

  sim.initSystem(sim.opt.solver, 1, 1);
  sim.setMode(SIM::STATIC);

  if (!sim.assembleSystem())
    ASSERT_EQ(0, 3);

  Vector load;
  sim.extractLoadVec(load);
  const std::array<int,25> data1 {0, 1, 1, 1, 0, 2, 2, 2, 0, 2, 2, 2, 0, 1, 1, 1};
  const std::array<int,25> data2 {1, 1, 1, 0, 2, 2, 2, 0, 2, 2, 2, 0, 1, 1, 1, 0};
  for (size_t i = 0; i < data1.size(); ++i) {
    ASSERT_EQ(load(i+1), data1[i]);
    ASSERT_EQ(load(i+26), -data2[i]);
  }
}
