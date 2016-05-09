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


/*! \brief Helper template for setting up a solver with the appropriate
           preconditioner type.
    \details We cannot instance using dynamic polymorphism, the solver need
             access to the real type for the preconditioner. We can however
             call the solver in the interface class scope afterwards.
 */
template<class Prec>
static Dune::InverseOperator<ISTL::Vec,ISTL::Vec>*
  setupWithPreType(const LinSolParams& solParams,
                   Dune::MatrixAdapter<ISTL::Mat,ISTL::Vec,ISTL::Vec>& op,
                   Dune::Preconditioner<ISTL::Vec,ISTL::Vec>& prec)
{
  Prec& pre = static_cast<Prec&>(prec);

  std::string type = solParams.getStringValue("type");
  double rtol = solParams.getDoubleValue("rtol");
  int maxits = solParams.getIntValue("maxits");
  const int verbosity = 2;
  if (type == "bcgs")
    return new Dune::BiCGSTABSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "cg")
    return new Dune::CGSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "minres")
    return new Dune::MINRESSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "gmres")
    return new Dune::RestartedGMResSolver<ISTL::Vec>(op, pre, rtol, solParams.getIntValue("gmres_restart_iterations"),
                                                     maxits, verbosity);

  return nullptr;
}


std::tuple<std::unique_ptr<ISTL::InverseOperator>,
           std::unique_ptr<ISTL::Preconditioner>> ISTLSolParams::setupPC(ISTL::Mat& A, ISTL::Operator& op)
{
  std::unique_ptr<ISTL::InverseOperator> solver;
  std::unique_ptr<ISTL::Preconditioner> pre;
  std::string prec = solParams.getBlock(0).getStringValue("pc");
  if (prec == "ilu") {
    if (solParams.getBlock(0).getIntValue("ilu_fill_level") == 0) {
      pre.reset(new Dune::SeqILU0<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1.0));
      solver.reset(setupWithPreType<Dune::SeqILU0<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, op, *pre));
    } else {
      pre.reset(new Dune::SeqILUn<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, solParams.getBlock(0).getIntValue("ilu_fill_level"), 1.0));
      solver.reset(setupWithPreType<Dune::SeqILUn<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, op, *pre));
    }
  } else if (prec == "sor") {
    pre.reset(new Dune::SeqSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0));
    solver.reset(setupWithPreType<Dune::SeqSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, op, *pre));
  } else if (prec == "ssor") {
    pre.reset(new Dune::SeqSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0));
    solver.reset(setupWithPreType<Dune::SeqSSOR<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, op, *pre));
  } else if (prec == "jacobi") {
    pre.reset(new Dune::SeqJac<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0));
    solver.reset(setupWithPreType<Dune::SeqJac<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, op, *pre));
  } else if (prec == "gs") {
    pre.reset(new Dune::SeqGS<ISTL::Mat,ISTL::Vec,ISTL::Vec>(A, 1, 1.0));
    solver.reset(setupWithPreType<Dune::SeqGS<ISTL::Mat,ISTL::Vec,ISTL::Vec>>(solParams, op, *pre));
  }
  else if (prec == "asm" || prec == "asmlu") {
    size_t nx = solParams.getBlock(0).getIntValue("nx");
    nx = std::max(1ul, nx);
    size_t ny = solParams.getBlock(0).getIntValue("ny");
    ny = std::max(1ul, ny);
    size_t nz = solParams.getBlock(0).getIntValue("nz");
    nz = std::max(1ul, nz);
    int overlap = solParams.getBlock(0).getIntValue("overlap");

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

      pre.reset(new Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec,
                                                Dune::AdditiveSchwarzMode, Dune::SuperLU<ISTL::Mat>>(A, ddofs));
      solver.reset(setupWithPreType<Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec,
                                                                Dune::AdditiveSchwarzMode, Dune::SuperLU<ISTL::Mat>>>(solParams, op, *pre));
    } else {
      Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec>::subdomain_vector ddofs(locSubdDofs.size());
      for (size_t i = 0; i < locSubdDofs.size(); ++i)
        ddofs[i].insert(locSubdDofs[i].begin(), locSubdDofs[i].end());

      pre.reset(new Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec>(A, ddofs));
      solver.reset(setupWithPreType<Dune::SeqOverlappingSchwarz<ISTL::Mat, ISTL::Vec>>(solParams, op, *pre));
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

    pre.reset(new AMG(op, crit, args));
    solver.reset(setupWithPreType<AMG>(solParams, op, *pre));
  }

  return std::make_tuple(std::move(solver), std::move(pre));
}
