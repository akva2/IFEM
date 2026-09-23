//==============================================================================
//!
//! \file TestMultigridTransfer.C
//!
//! \date Sep 23 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for the geometric multigrid transfer operators.
//!
//==============================================================================

#include "MultigridTransfer.h"

#include "ASMbase.h"
#include "ASMunstruct.h"
#include "SAM.h"
#include "SIM2D.h"
#include "SparseMatrix.h"
#include "Vec3.h"

#include "Catch2Support.h"

#include <memory>
#include <numeric>
#include <vector>

#ifdef HAS_LRSPLINE

namespace {

const char* infile = "src/SIM/Test/refdata/mgtransfer.xinp";


//! \brief Reads the test model into a simulator and preprocesses it.
//! \param sim The simulator to set up
//! \param[in] mesh If given, adopt the mesh of this simulator
//!
//! \details Passing \a mesh mimics what the multigrid driver does when it
//! snapshots a level before the adaptive loop refines it further: a fresh
//! simulator reads the same input file, is given a private copy of the mesh
//! as it is at that point, and is then preprocessed on that mesh.
bool setup (SIM2D& sim, const SIM2D* mesh = nullptr)
{
  sim.opt.discretization = ASM::LRSpline;
  if (!sim.read(infile))
    return false;

  if (mesh)
    for (size_t i = 0; i < sim.getFEModel().size(); i++)
      if (!sim.getFEModel()[i]->copyMeshFrom(*mesh->getFEModel()[i]))
        return false;

  return sim.preprocess();
}


//! \brief Re-establishes the FE data of a simulator after a refinement.
bool reprocess (SIM2D& sim)
{
  sim.clearProperties();
  return sim.read(infile) && sim.preprocess();
}


//! \brief Returns one coordinate of the control points, indexed by equation.
//! \param[in] sim The simulator to extract the control points of
//! \param[in] dir Coordinate direction to extract
//!
//! \details Refinement leaves the geometry map unchanged, so the control
//! points of the coarse mesh and those of the fine mesh represent the same
//! function. They are therefore an exact test case for the prolongation.
Vector controlPoints (const SIM2D& sim, int dir)
{
  const ASMbase* pch = sim.getFEModel().front();
  Vector X(sim.getSAM()->getNoEquations());
  for (size_t i = 1; i <= pch->getNoNodes(); i++)
    if (int eq = sim.getSAM()->getEquation(pch->getNodeID(i),1); eq > 0)
      X(eq) = pch->getCoord(i)[dir];

  return X;
}


//! \brief Multiplies a prolongation operator with a vector.
Vector prolongate (const SparseMatrix& P, const Vector& u)
{
  Vector result(P.rows());
  for (const auto& [ij,v] : P.getValues())
    result(ij.first) += v*u(ij.second);

  return result;
}

}


TEST_CASE("TestMultigridTransfer.Identity")
{
  // Two simulators on the same mesh must give the identity as prolongation.
  // This also covers ASMbase::copyMeshFrom, since the second simulator gets
  // its mesh from the first one rather than from the input file.
  SIM2D coarse(1), fine(1);
  REQUIRE(setup(coarse));
  REQUIRE(setup(fine,&coarse));

  const size_t neq = coarse.getSAM()->getNoEquations();
  REQUIRE(fine.getSAM()->getNoEquations() == neq);

  std::unique_ptr<SparseMatrix> P = MG::prolongation(coarse,fine,{"A",1,0});
  REQUIRE(P != nullptr);
  REQUIRE(P->rows() == neq);
  REQUIRE(P->cols() == neq);

  for (size_t i = 1; i <= neq; i++)
    for (size_t j = 1; j <= neq; j++)
      REQUIRE_THAT((*P)(i,j), WithinAbs(i == j ? 1.0 : 0.0, 1e-12));
}


TEST_CASE("TestMultigridTransfer.Exactness")
{
  // The coarse space is contained in the fine one, so the prolongation must
  // reproduce a coarse function exactly in the fine basis. The geometry map
  // is such a function, and refinement leaves it unchanged, which makes the
  // control points a reference the operator can be checked against.
  //
  // Local refinement is the interesting case: it leaves overloaded elements
  // behind, which a plain element-by-element change of basis cannot resolve,
  // whereas uniform refinement does not. Both are covered here.
  const bool uniform = GENERATE(true,false);

  SECTION(uniform ? "Uniform refinement" : "Local refinement") {
  SIM2D coarse(1), fine(1);
  REQUIRE(setup(coarse));
  REQUIRE(setup(fine,&coarse));

  // With the structured mesh strategy, options[2] == 2, the indices in
  // RefineData::elements name basis functions rather than elements.
  for (int step = 0; step < (uniform ? 1 : 3); step++)
  {
    LR::RefineData prm;
    prm.options = {10,1,2};
    size_t nfunc = fine.getFEModel().front()->getNoNodes();
    for (size_t i = 0; i < nfunc; i++)
      if (uniform || i%7 == 0)
        prm.elements.push_back(i);

    REQUIRE(fine.refine(prm));
    REQUIRE(reprocess(fine));
  }

  REQUIRE(fine.getSAM()->getNoEquations() >
          coarse.getSAM()->getNoEquations());

  std::unique_ptr<SparseMatrix> P = MG::prolongation(coarse,fine,{"A",1,0});
  REQUIRE(P != nullptr);
  REQUIRE(P->rows() == (size_t)fine.getSAM()->getNoEquations());
  REQUIRE(P->cols() == (size_t)coarse.getSAM()->getNoEquations());

  for (int dir = 0; dir < 2; dir++)
  {
    Vector Xf = prolongate(*P,controlPoints(coarse,dir));
    Vector Xref = controlPoints(fine,dir);
    REQUIRE(Xf.size() == Xref.size());
    for (size_t i = 1; i <= Xf.size(); i++)
      REQUIRE_THAT(Xf(i), WithinAbs(Xref(i), 1e-10));
  }

  // The rows of a prolongation between nested spline spaces sum to one,
  // since both bases are partitions of unity.
  Vector rowSum(P->rows());
  for (const auto& [ij,v] : P->getValues())
    rowSum(ij.first) += v;
  for (size_t i = 1; i <= rowSum.size(); i++)
    REQUIRE_THAT(rowSum(i), WithinAbs(1.0, 1e-10));
  }
}

#endif
