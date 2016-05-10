// $Id$
//==============================================================================
//!
//! \file ISTLSolParams.C
//!
//! \date Mar 10 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear solver parameters for ISTL.
//!
//==============================================================================

#include "ISTLSolParams.h"
#include "ASMstruct.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#include "SAMpatch.h"
#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/matrixmatrix.hh>


/*! \brief Helper template for setting up a solver with the appropriate
           preconditioner type.
    \details We cannot instance using dynamic polymorphism, the solver need
             access to the real type for the preconditioner. We can however
             call the solver in the interface class scope afterwards.
 */
template<class Prec>
static Dune::InverseOperator<ISTL::Vec,ISTL::Vec>*
  setupWithPreType(const LinSolParams& solParams,
                   ISTL::Operator& op,
                   ISTL::Preconditioner& prec)
{
  Prec& pre = static_cast<Prec&>(prec);

  std::string type = solParams.getStringValue("type");
  double rtol = solParams.getDoubleValue("rtol");
  int maxits = solParams.getIntValue("maxits");
  int verbosity = solParams.getIntValue("verbosity");
  int restart = solParams.getIntValue("gmres_restart_iterations");
  if (type == "bcgs")
    return new Dune::BiCGSTABSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "cg")
    return new Dune::CGSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "minres")
    return new Dune::MINRESSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "gmres")
    return new Dune::RestartedGMResSolver<ISTL::Vec>(op, pre, rtol, restart,
                                                     maxits, verbosity);

  return nullptr;
}


#include "Profiler.h"


namespace ISTL {

BlockPreconditioner::BlockPreconditioner(const ISTL::Mat& A,
                                         const DomainDecomposition& dd_) :
  dd(dd_)
{
  blocks.resize(dd.getNoBlocks()*dd.getNoBlocks());
  size_t k = 0;
  for (size_t i = 0; i < dd.getNoBlocks(); ++i)
    for (size_t j = 0; j < dd.getNoBlocks(); ++j)
      extractBlock(blocks[k++], A, dd.getBlockEqs(i), dd.getBlockEqs(j));

  // Scale blocks[2] with inverse diagonal of blocks[0]
  for (auto row = blocks[2].begin(); row != blocks[2].end(); ++row)
    for (auto it = row->begin(); it != row->end(); ++it)
      *it /= blocks[0][it.index()][it.index()];

  // blocks[1] = D - C*diag(A)^-1*B
  ISTL::Mat SP;
  Dune::matMultMat(SP, blocks[2], blocks[1]);
  subtractMatrices(blocks[1], blocks[3], SP);

  // Clear unnecessary block matrices
  blocks.resize(2);

  blockPre.resize(dd.getNoBlocks());
  blockOp.resize(dd.getNoBlocks());
}


void BlockPreconditioner::pre(ISTL::Vec& x, ISTL::Vec& b)
{
  ISTL::Vec tempx, tempb;
  for (size_t block = 0; block < dd.getNoBlocks(); ++block) {
    ISTL::Vec tempx, tempb;
    tempx.resize(dd.getBlockEqs(block).size());
    tempb.resize(dd.getBlockEqs(block).size());
    size_t i = 0;
    for (auto& it : dd.getBlockEqs(block))
      tempx[i] = x[it-1], tempb[i++] = b[it-1];

    blockPre[block]->pre(tempx, tempb);

    i = 0;
    for (auto& it : dd.getBlockEqs(block))
      x[it-1] = tempx[i], b[it-1] = tempb[i++];
  }
}


void BlockPreconditioner::apply(ISTL::Vec& v, const ISTL::Vec& d)
{
  // backwards to allow for a non-diagonal preconditioner later
  for (int block = dd.getNoBlocks()-1; block >= 0; --block) {
    ISTL::Vec tempx, tempb;
    tempx.resize(dd.getBlockEqs(block).size());
    tempb.resize(dd.getBlockEqs(block).size());
    size_t i = 0;
    for (auto& it : dd.getBlockEqs(block))
      tempb[i++] = d[it-1];

    tempx = 0;
    blockPre[block]->apply(tempx, tempb);

    i = 0;
    for (auto& it : dd.getBlockEqs(block))
      v[it-1] = tempx[i++];
  }
}


void BlockPreconditioner::post(ISTL::Vec& x)
{
  // Not necessary?
  return;
}


void BlockPreconditioner::extractBlock(ISTL::Mat& B, const ISTL::Mat& A,
                                       const std::set<int>& eqs_row1,
                                       const std::set<int>& eqs_col1)
{
  std::vector<int> eqs_row(eqs_row1.begin(), eqs_row1.end());
  std::vector<int> eqs_col(eqs_col1.begin(), eqs_col1.end());
  size_t sum=0;
  std::vector<std::set<int>> adj;
  adj.resize(eqs_row.size());
  size_t i = 0;

  std::vector<int> eq2bc(A.M(), -1);
  for (auto& it3 : eqs_row) {
    auto it = A.begin()+it3-1;
    for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
      if (eq2bc[it2.index()] == -1) {
        auto pos = eqs_col1.find(it2.index()+1);
        if (pos != eqs_col1.end())
          eq2bc[it2.index()] = std::distance(eqs_col1.begin(), pos);
        else
          eq2bc[it2.index()] = -2;
      }
      if (eq2bc[it2.index()] > -1) {
        adj[i].insert(eq2bc[it2.index()]);
        ++sum;
      }
    }
    ++i;
  }

  B.setSize(eqs_row.size(), eqs_col.size(), sum);
  B.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < adj.size(); i++)
    B.setrowsize(i,adj[i].size());
  B.endrowsizes();

  for (size_t i = 0; i < adj.size(); i++)
    for (auto& it : adj[i])
      B.addindex(i, it);

  B.endindices();
  B = 0;
  for (auto it = B.begin(); it != B.end(); ++it)
    for (auto it2 = it->begin(); it2 != it->end(); ++it2)
      *it2 = A[eqs_row[it.index()]-1][eqs_col[it2.index()]-1];
}


void BlockPreconditioner::subtractMatrices(ISTL::Mat& A, const ISTL::Mat& B,
                                           const ISTL::Mat& C)
{
  size_t sum=0;
  std::vector<std::set<int>> adj;
  adj.resize(B.N());
  for (auto row = B.begin(); row != B.end(); ++row) {
    for (auto it = row->begin(); it != row->end(); ++it)
      adj[row.index()].insert(it.index());
  }
  for (auto row = C.begin(); row != C.end(); ++row) {
    for (auto it = row->begin(); it != row->end(); ++it)
      adj[row.index()].insert(it.index());
    sum += adj[row.index()].size();
  }

  A.setSize(B.N(), B.M(), sum);
  A.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < adj.size(); i++)
    A.setrowsize(i,adj[i].size());
  A.endrowsizes();

  for (size_t i = 0; i < adj.size(); i++)
    for (auto& it : adj[i])
      A.addindex(i, it);

  A.endindices();
  A = 0;
  for (auto row = B.begin(); row != B.end(); ++row)
    for (auto it = row->begin(); it != row->end(); ++it)
      A[row.index()][it.index()] = *it;

  for (auto row = C.begin(); row != C.end(); ++row)
    for (auto it = row->begin(); it != row->end(); ++it)
      A[row.index()][it.index()] -= *it;
}

} // namespace ISTL


ISTL::Preconditioner* ISTLSolParams::setupPCInternal(ISTL::Mat& A,
                                                     ISTL::Operator& op,
                                                     size_t block)
{
  std::string prec = solParams.getBlock(block).getStringValue("pc");
  if (prec == "ilu") {
    int fill_level = solParams.getBlock(block).getIntValue("ilu_fill_level");
    if (fill_level == 0)
      return new Dune::SeqILU0<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1.0);
    else
      return new Dune::SeqILUn<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, fill_level, 1.0);
  } else if (prec == "sor")
    return new Dune::SeqSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0);
  else if (prec == "ssor")
    return new Dune::SeqSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0);
  else if (prec == "jacobi")
    return new Dune::SeqJac<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0);
  else if (prec == "gs")
    return new Dune::SeqGS<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0);
  else if (prec == "asm" || prec == "asmlu") {
    size_t nx = solParams.getBlock(block).getIntValue("nx");
    nx = std::max(1ul, nx);
    size_t ny = solParams.getBlock(block).getIntValue("ny");
    ny = std::max(1ul, ny);
    size_t nz = solParams.getBlock(block).getIntValue("nz");
    nz = std::max(1ul, nz);
    int overlap = solParams.getBlock(block).getIntValue("overlap");

    const SAMpatch* samp = adm.dd.getSAM();
    std::vector<std::set<int>>  locSubdDofs(nx*ny*nz*samp->getNoPatches());
    size_t d = 0;
    for (const auto& it : *samp) {
      const ASMstruct* pch = dynamic_cast<const ASMstruct*>(it);
      if (!pch)
        break;
      int n1, n2, n3;
      pch->getNoStructElms(n1,n2,n3);
      const_cast<DomainDecomposition&>(adm.dd).calcAppropriateGroups(n1, n2, n3, nx, ny, nz, overlap);
      for (size_t g = 0; g < adm.dd.getNoSubdomains(); ++g, ++d) {
        for (const auto& iEl : adm.dd[g]) {
          IntVec eqns;
          samp->getElmEqns(eqns, it->getElmID(iEl+1));
          for (auto& it : eqns) {
            if (it > 0)
              locSubdDofs[d].insert(it-1);
          }
        }
      }
    }
    if (prec == "asmlu") {
      Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec,
                                 Dune::AdditiveSchwarzMode,
                                 Dune::SuperLU<ISTL::Mat>>::subdomain_vector ddofs(locSubdDofs.size());
      for (size_t i = 0; i < locSubdDofs.size(); ++i)
        for (const auto& it : locSubdDofs[i])
          ddofs[i].insert(it);

      return new Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec,
                                             Dune::AdditiveSchwarzMode,
                                             Dune::SuperLU<ISTL::Mat>>(A, ddofs);
    } else {
      Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec>::subdomain_vector ddofs(locSubdDofs.size());
      for (size_t i = 0; i < locSubdDofs.size(); ++i)
        ddofs[i].insert(locSubdDofs[i].begin(), locSubdDofs[i].end());

      return new Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec>(A, ddofs);
    }
  } else if (prec == "amg") {
    // The coupling metric used in the AMG
    typedef Dune::Amg::FirstDiagonal CouplingMetric;
    // The coupling criterion used in the AMG
    typedef Dune::Amg::SymmetricCriterion<ISTL::Mat, CouplingMetric> CritBase;
    // The coarsening criterion used in the AMG
    typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;
    Criterion crit;
    typedef Dune::Amg::AMG<ISTL::Operator, ISTL::Vec, Dune::SeqSSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>> AMG;
    AMG::SmootherArgs args;
    args.relaxationFactor = 1.0;
    args.iterations = std::max(1, solParams.getBlock(0).getIntValue("multigrid_no_smooth"));

    return new AMG(op, crit, args);
  }

  return nullptr;
}


std::tuple<std::unique_ptr<ISTL::InverseOperator>,
           std::unique_ptr<ISTL::Preconditioner>,
           std::unique_ptr<ISTL::Operator>> ISTLSolParams::setupPC(ISTL::Mat& A)
{
  std::unique_ptr<ISTL::InverseOperator> solver;
  std::unique_ptr<ISTL::Preconditioner> pre;
  std::unique_ptr<ISTL::Operator> op;

  if (solParams.getNoBlocks() > 2) {
    std::cerr << "*** ISTL ** More than two blocks are not implemented." << std::endl;
    return std::make_tuple(nullptr, nullptr, nullptr);
  }

  op.reset(new ISTL::Operator(A));

  if (solParams.getNoBlocks() > 1) {
    ISTL::BlockPreconditioner* bpre = new ISTL::BlockPreconditioner(A, adm.dd);
    pre.reset(bpre);
    for (size_t i = 0; i < adm.dd.getNoBlocks(); ++i)
      bpre->getBlockPre(i).reset(setupPCInternal(bpre->getBlock(i),
                                                 bpre->getBlockOp(i), i));
    solver.reset(setupWithPreType<ISTL::BlockPreconditioner>(solParams, *op, *pre));
  } else {
    pre.reset(setupPCInternal(A, *op, 0));
    std::string prec = solParams.getBlock(0).getStringValue("pc");
    if (prec == "ilu") {
      if (solParams.getBlock(0).getIntValue("ilu_fill_level") == 0)
        solver.reset(setupWithPreType<Dune::SeqILU0<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, *op, *pre));
      else
        solver.reset(setupWithPreType<Dune::SeqILUn<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, *op, *pre));
    } else if (prec == "sor")
      solver.reset(setupWithPreType<Dune::SeqSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, *op, *pre));
    else if (prec == "ssor")
      solver.reset(setupWithPreType<Dune::SeqSSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, *op, *pre));
    else if (prec == "jacobi")
      solver.reset(setupWithPreType<Dune::SeqJac<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, *op, *pre));
    else if (prec == "gs")
      solver.reset(setupWithPreType<Dune::SeqGS<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, *op, *pre));
    else if (prec == "asm")
      solver.reset(setupWithPreType<Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec>>(solParams, *op, *pre));
    else if (prec == "asmlu")
        solver.reset(setupWithPreType<Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec,
                                                                  Dune::AdditiveSchwarzMode, Dune::SuperLU<ISTL::Mat>>>(solParams, *op, *pre));
    else if (prec == "amg") {
      // The coupling metric used in the AMG
      typedef Dune::Amg::FirstDiagonal CouplingMetric;
      // The coupling criterion used in the AMG
      typedef Dune::Amg::SymmetricCriterion<ISTL::Mat, CouplingMetric> CritBase;
      // The coarsening criterion used in the AMG
      typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;
      Criterion crit;
      typedef Dune::Amg::AMG<ISTL::Operator, ISTL::Vec, Dune::SeqSSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>> AMG;
      solver.reset(setupWithPreType<AMG>(solParams, *op, *pre));
    }
  }

  return std::make_tuple(std::move(solver), std::move(pre), std::move(op));
}
