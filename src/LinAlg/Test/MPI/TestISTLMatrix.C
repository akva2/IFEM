//==============================================================================
//!
//! \file TestISTLMatrix.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for parallel ISTL matrices
//!
//==============================================================================

#include "AlgEqSystem.h"
#include "ISTLMatrix.h"
#include "ISTLSolParams.h"
#include "SIM2D.h"
#include "ASMs2D.h"
#include "IntegrandBase.h"
#include "SAM.h"
#include "IFEM.h"

#include "gtest/gtest.h"

#include <fstream>


typedef std::vector<int> IntVec;

static IntVec readIntVector(const std::string& file)
{
  std::vector<int> result;
  std::ifstream f(file);
  size_t size;
  f >> size;
  result.resize(size);
  for (size_t j=0;j<size;++j)
    f >> result[j];

  return result;
}


class DummyIntegrandForDummies : public IntegrandBase {
public:
  DummyIntegrandForDummies(unsigned short int n = 0) : IntegrandBase(n) {}
};


class InspectMatrixSIM : public SIM2D {
public:
  InspectMatrixSIM(unsigned char n1 = 2, bool check = false) :
    SIM2D(n1, check) { myProblem = new DummyIntegrandForDummies; }

  InspectMatrixSIM(const std::vector<unsigned char>& n, bool check = false) :
    SIM2D(n, check) { myProblem = new DummyIntegrandForDummies; }

  SystemMatrix* getMatrix() { return myEqSys->getMatrix(0); }
  SystemVector* getVector() { return myEqSys->getVector(0); }
};


TEST(TestISTLMatrix, AssembleMPI)
{
  InspectMatrixSIM sim(1);
  sim.read("src/LinAlg/Test/refdata/petsc_test.xinp");
  sim.opt.solver = SystemMatrix::ISTL;
  sim.preprocess();
  sim.initSystem(SystemMatrix::ISTL);

  Matrix stencil(4,4);
  stencil(1,1) = stencil(2,2) = stencil(3,3) = stencil(4,4) = 1.0;

  for (int iel = 1; iel <= sim.getSAM()->getNoElms(); ++iel)
    sim.getMatrix()->assemble(stencil, *sim.getSAM(), iel);

  sim.getMatrix()->beginAssembly();
  sim.getMatrix()->endAssembly();

  // now inspect the matrix
  const ProcessAdm& adm = sim.getProcessAdm();
  ISTL::Mat& mat = static_cast<ISTLMatrix*>(sim.getMatrix())->getMatrix();
  ISTL::Vec b(mat.N()), b2(mat.N()), b3(mat.N()), b4(mat.N());
  double dot = 0.0;

  try {
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    Comm comm(*adm.getCommunicator());
    comm.indexSet().beginResize();
    typedef Dune::ParallelLocalIndex<Dune::OwnerOverlapCopyAttributeSet::AttributeSet> LI;
    for (size_t i = 0; i < adm.dd.getMLGEQ().size(); ++i) {
      int gid = adm.dd.getGlobalEq(i+1);
      comm.indexSet().add(gid-1, LI(i, gid >= adm.dd.getMinEq() ?
                                       Dune::OwnerOverlapCopyAttributeSet::owner :
                                       Dune::OwnerOverlapCopyAttributeSet::overlap));
    }
    comm.indexSet().endResize();
    comm.remoteIndices().setIncludeSelf(true);
    comm.remoteIndices().template rebuild<false>();

    ISTL::ParMatrixAdapter op(mat, comm);
    Dune::OverlappingSchwarzScalarProduct<ISTL::Vec,Comm> sp(comm);

    b = 1.0;
    op.apply(b, b2);
    b3 = 1.0;
    op.applyscaleadd(0.5, b, b3);
    dot = sp.dot(b2, b2);
//    ISTL::LU lpc(new ISTL::LUType(mat));
//    ISTL::LU lpc(new ISTL::LUType(mat));
    ISTL::GJ lpc(mat, 1, 1.0);
    ISTL::pASM<ISTL::Vec,ISTL::Vec,Comm,ISTL::GJ> pc(lpc,comm);
    Dune::RestartedGMResSolver<ISTL::Vec> cg(op, sp, pc, 1e-12, 1000, 1000, adm.getProcId()==0?2:0);
//    Dune::CGSolver<ISTL::Vec> cg(op, sp, pc, 1e-12, 1000, adm.getProcId()==0?2:0);
//    Dune::MINRESSolver<ISTL::Vec> cg(op, sp, pc, 1e-12, 1000, adm.getProcId()==0?2:0);
//    Dune::BiCGSTABSolver<ISTL::Vec> cg(op, sp, pc, 1e-12, 1000, adm.getProcId()==0?2:0);
    Dune::InverseOperatorResult r;
    b4 = 0;
    b = b2;
    cg.apply(b4, b, 1e-14, r);
    IFEM::cout << b4 << std::endl;
  } catch (Dune::ISTLError e) {
    std::cerr << e << std::endl;
    ASSERT_TRUE(false);
  }

  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");
  for (size_t i = 1; i <= adm.dd.getMLGEQ().size(); ++i) {
    ASSERT_FLOAT_EQ(v[adm.dd.getGlobalEq(i)-1], b2[i-1]);
    ASSERT_FLOAT_EQ(1.0+0.5*v[adm.dd.getGlobalEq(i)-1], b3[i-1]);
  }
  ASSERT_FLOAT_EQ(dot, 196);
}
