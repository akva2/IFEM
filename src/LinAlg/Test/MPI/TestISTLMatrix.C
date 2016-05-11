//==============================================================================
//!
//! \file TestPETScMatrix.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for parallel PETSc matrices
//!
//==============================================================================

#include "AlgEqSystem.h"
#include "ISTLMatrix.h"
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


TEST(TestISTLMatrix, Assemble)
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

  IFEM::cout << "past haere" << std::endl;

  // now inspect the matrix
  const ProcessAdm& adm = sim.getProcessAdm();
  ISTL::Mat& mat = static_cast<ISTLMatrix*>(sim.getMatrix())->getMatrix();
  ISTL::Vec b(mat.N()), b2(mat.N());

  try {
  Dune::IndexInfoFromGrid<int, int> index;
  for (size_t i = 0; i < adm.dd.getMLGEQ().size(); ++i) {
    int gid = adm.dd.getGlobalEq(i+1);

    index.addLocalIndex(std::make_tuple(gid-1, i, gid >= adm.dd.getMinEq() ? 
						  Dune::OwnerOverlapCopyAttributeSet::owner :
                                                  Dune::OwnerOverlapCopyAttributeSet::overlap));
  }
  Dune::OwnerOverlapCopyCommunication<int,int> comm(index, *adm.getCommunicator());
  IFEM::cout << "past hare" << std::endl;
  ISTL::ParMatrixAdapter op(mat, comm);

  IFEM::cout << "past ehre" << std::endl;

  b = 1.0;
  op.apply(b, b2);
  for(auto& it : b2)
    IFEM::cout << it << " ";
  IFEM::cout << std::endl;
  comm.addOwnerOverlapToAll(b2,b2);
  for(auto& it : b2)
    IFEM::cout << it << " ";
  IFEM::cout << std::endl;
  } catch (Dune::ISTLError e) {
    IFEM::cout << e << std::endl;
  }

  if (adm.getProcId() < 1) {
  IntVec v = readIntVector("src/LinAlg/Test/refdata/petsc_matrix_diagonal.ref");
  for (int i = adm.dd.getMinEq(); i <= adm.dd.getMaxEq(); ++i) {
    IFEM::cout << "check " << i-1 << " " << v[i-1] << " " << i-adm.dd.getMinEq() << " " << b2[i-adm.dd.getMinEq()] << std::endl;
    ASSERT_FLOAT_EQ(v[i-1], b2[i-adm.dd.getMinEq()]);
  }
  }
}
