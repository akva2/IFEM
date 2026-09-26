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
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "ASMunstruct.h"
#include "SAM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SparseMatrix.h"
#include "Vec3.h"

#include "Catch2Support.h"

#include <memory>
#include <numeric>
#include <vector>

namespace {

//! \brief Returns one coordinate of the control points, indexed by equation.
//! \param[in] sim The simulator to extract the control points of
//! \param[in] dir Coordinate direction to extract
//!
//! \details Refinement leaves the geometry map unchanged, so the control
//! points of the coarse mesh and those of the fine mesh represent the same
//! function. They are therefore an exact test case for the prolongation.
Vector controlPoints (const SIMbase& sim, int dir)
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


//! \brief Raises the order of and refines every patch of a simulator.
//! \param sim The simulator whose patches to work on
//! \param[in] nInsert Number of knots to insert into each knot span
//! \param[in] raise Orders to raise the basis by first
//!
//! \details This is what the \a raiseorder and \a refine tags of the input
//! file do, carried out here instead so that the coarse and the fine level
//! can be read from one and the same file and differ only in how far it was
//! taken. Both are uniform over the patch and over its directions.
bool refineAll (SIMbase& sim, int nInsert, int raise = 0)
{
  for (ASMbase* pch : sim.getFEModel())
    if (ASMs2D* p2 = dynamic_cast<ASMs2D*>(pch); p2)
    {
      if (raise > 0 && !p2->raiseOrder(raise,raise))
        return false;
      for (int d = 0; d < 2 && nInsert > 0; d++)
        if (!p2->uniformRefine(d,nInsert))
          return false;
    }
    else if (ASMs3D* p3 = dynamic_cast<ASMs3D*>(pch); p3)
    {
      if (raise > 0 && !p3->raiseOrder(raise,raise,raise))
        return false;
      for (int d = 0; d < 3 && nInsert > 0; d++)
        if (!p3->uniformRefine(d,nInsert))
          return false;
    }
    else
      return false;

  return true;
}

}


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

  std::unique_ptr<MG::Prolongation> res = MG::prolongation(coarse,fine,{"A",1,0});
  REQUIRE(res != nullptr);
  const SparseMatrix* P = res->P.get();
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

  std::unique_ptr<MG::Prolongation> res = MG::prolongation(coarse,fine,{"A",1,0});
  REQUIRE(res != nullptr);
  const SparseMatrix* P = res->P.get();
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


namespace {

//! \brief Reads a tensor product test model and preprocesses it.
//! \param sim The simulator to set up
//! \param[in] file The input file to read
//! \param[in] nRef Knots inserted into each knot span before preprocessing
//! \param[in] raise Orders the basis is raised by before preprocessing
//!
//! \details The levels of a tensor product hierarchy are the same geometry
//! read as many times, refined to a different depth each time. Doing the
//! refinement here rather than in the input file is what lets them share one.
template<class Dim>
bool setupTensor (Dim& sim, const char* file, int nRef, int raise = 0)
{
  sim.opt.discretization = ASM::Spline;
  if (!sim.read(file))
    return false;

  return refineAll(sim,nRef,raise) && sim.preprocess();
}

}


TEST_CASE("TestMultigridTransfer.TensorIdentity")
{
  // Two simulators on the same tensor mesh must give the identity.
  SIM2D coarse(1), fine(1);
  REQUIRE(setupTensor(coarse,"src/SIM/Test/refdata/mgtensor.xinp",0));
  REQUIRE(setupTensor(fine,"src/SIM/Test/refdata/mgtensor.xinp",0));

  const size_t neq = coarse.getSAM()->getNoEquations();
  REQUIRE(fine.getSAM()->getNoEquations() == neq);

  std::unique_ptr<MG::Prolongation> res = MG::prolongation(coarse,fine,{"A",1,0});
  REQUIRE(res != nullptr);
  const SparseMatrix* P = res->P.get();
  REQUIRE(P != nullptr);
  REQUIRE(P->rows() == neq);
  REQUIRE(P->cols() == neq);

  for (size_t i = 1; i <= neq; i++)
    for (size_t j = 1; j <= neq; j++)
      REQUIRE_THAT((*P)(i,j), WithinAbs(i == j ? 1.0 : 0.0, 1e-12));
}


TEST_CASE("TestMultigridTransfer.TensorExactness")
{
  // The coarse space is contained in the fine one, so the prolongation must
  // reproduce a coarse function exactly in the fine basis. The geometry map
  // is such a function, and refinement leaves it unchanged, which makes the
  // control points a reference the operator can be checked against.
  const int nsd = GENERATE(2,3);

  SECTION(nsd == 2 ? "Surface" : "Volume") {
  const char* file = nsd == 2 ? "src/SIM/Test/refdata/mgtensor.xinp"
                              : "src/SIM/Test/refdata/mgtensor3D.xinp";

  std::unique_ptr<SIMbase> coarse, fine;
  if (nsd == 2) {
    coarse = std::make_unique<SIM2D>(1);
    fine   = std::make_unique<SIM2D>(1);
    REQUIRE(setupTensor(static_cast<SIM2D&>(*coarse),file,0));
    REQUIRE(setupTensor(static_cast<SIM2D&>(*fine),file,2));
  }
  else {
    coarse = std::make_unique<SIM3D>(1);
    fine   = std::make_unique<SIM3D>(1);
    REQUIRE(setupTensor(static_cast<SIM3D&>(*coarse),file,0));
    REQUIRE(setupTensor(static_cast<SIM3D&>(*fine),file,2));
  }

  REQUIRE(fine->getSAM()->getNoEquations() >
          coarse->getSAM()->getNoEquations());

  std::unique_ptr<MG::Prolongation> res =
    MG::prolongation(*coarse,*fine,{"A",1,0});
  REQUIRE(res != nullptr);
  const SparseMatrix* P = res->P.get();
  REQUIRE(P != nullptr);
  REQUIRE(P->rows() == (size_t)fine->getSAM()->getNoEquations());
  REQUIRE(P->cols() == (size_t)coarse->getSAM()->getNoEquations());

  for (int dir = 0; dir < nsd; dir++)
  {
    Vector Xf = prolongate(*P,controlPoints(*coarse,dir));
    Vector Xref = controlPoints(*fine,dir);
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


TEST_CASE("TestMultigridTransfer.TensorOrderReduction")
{
  // Spaces of a different order are not nested, so what comes out are the
  // two factors of the projection standing in for a change of basis. Both
  // bases contain the constants, which the projection must therefore carry
  // over unchanged: the constant one solves M x = B 1 with x the constant
  // one as well, which is to say the two have the same row sums.
  SIM2D coarse(1), fine(1);
  REQUIRE(setupTensor(coarse,"src/SIM/Test/refdata/mgtensor.xinp",0));
  REQUIRE(setupTensor(fine,"src/SIM/Test/refdata/mgtensor.xinp",0,1));

  std::unique_ptr<MG::Prolongation> res = MG::prolongation(coarse,fine,{"A",1,0});
  REQUIRE(res != nullptr);
  REQUIRE(res->isProjection());
  REQUIRE(res->P == nullptr);
  REQUIRE(res->B->rows() == (size_t)fine.getSAM()->getNoEquations());
  REQUIRE(res->B->cols() == (size_t)coarse.getSAM()->getNoEquations());
  REQUIRE(res->mass->rows() == res->B->rows());

  Vector rowB(res->B->rows()), rowM(res->mass->rows());
  for (const auto& [ij,v] : res->B->getValues())
    rowB(ij.first) += v;
  for (const auto& [ij,v] : res->mass->getValues())
    rowM(ij.first) += v;

  for (size_t i = 1; i <= rowB.size(); i++)
    REQUIRE_THAT(rowB(i), WithinAbs(rowM(i), 1e-12));
}
