// $Id$
//==============================================================================
//!
//! \file ISTLMatrix.C
//!
//! \date Mar 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Representation of the system matrix in ISTL format.
//!
//==============================================================================

#include "ISTLMatrix.h"
#include "LinSolParams.h"
#include "LinAlgInit.h"
#include "ProcessAdm.h"
#include "SIMenums.h"
#include "ASMstruct.h"
#include "Profiler.h"
#include "SAMpatch.h"
#include "DomainDecomposition.h"
#include <dune/common/version.hh>
#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/paamg/amg.hh>
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/paamg/twolevelmethod.hh>
#endif


ISTLVector::ISTLVector(const ProcessAdm& padm) : adm(padm)
{
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, size_t n) : adm(padm)
{
  x.resize(n);
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ProcessAdm& padm, const Real* values, size_t n) : adm(padm)
{
  x.resize(n);
  this->restore(values);
  LinAlgInit::increfs();
}


ISTLVector::ISTLVector(const ISTLVector& vec) : adm(vec.adm)
{
  x = vec.x;
  LinAlgInit::increfs();
}


ISTLVector::~ISTLVector()
{
  LinAlgInit::decrefs();
}


void ISTLVector::init(Real value)
{
  StdVector::init(value);
  x = value;
}

size_t ISTLVector::dim() const
{
  return x.size();
}


void ISTLVector::redim(size_t n)
{
  x.resize(n);
  StdVector::redim(n);
}


bool ISTLVector::beginAssembly()
{
  for (size_t i = 0; i < size(); ++i)
    x[i] = (*this)[i];

  return true;
}


bool ISTLVector::endAssembly()
{
  return true;
}


Real ISTLVector::L1norm() const
{
  return x.one_norm();
}


Real ISTLVector::L2norm() const
{
  return x.two_norm();
}


Real ISTLVector::Linfnorm() const
{
  return x.infinity_norm();
}


ISTLMatrix::ISTLMatrix (const ProcessAdm& padm, const LinSolParams& spar,
                        LinAlg::LinearSystemType ltype) :
 SparseMatrix(SUPERLU, 1), op(A), adm(padm),
 solParams(spar), linsysType(ltype)
{
  LinAlgInit::increfs();

  setParams = true;
  nLinSolves = 0;
}


ISTLMatrix::ISTLMatrix (const ISTLMatrix& B) :
  op(A), adm(B.adm), solParams(B.solParams)
{
  A = B.A;

  LinAlgInit::increfs();

  setParams = true;
}


ISTLMatrix::~ISTLMatrix ()
{
  LinAlgInit::decrefs();
}


void ISTLMatrix::initAssembly (const SAM& sam, bool b)
{
  SparseMatrix::initAssembly(sam, b);
  SparseMatrix::preAssemble(sam, b);

  std::vector<std::set<int>> dofc;
  sam.getDofCouplings(dofc);
  samp = dynamic_cast<const SAMpatch*>(&sam);

  // Set correct number of rows and columns for matrix.
  size_t sum = 0;
  for (const auto& it : dofc)
    sum += it.size();

  A.setSize(rows(), cols(), sum);
  A.setBuildMode(Mat::row_wise);

  for (auto row = A.createbegin(); row != A.createend(); ++row)
    for (const auto& it : dofc[row.index()])
      row.insert(it-1);
}

bool ISTLMatrix::beginAssembly()
{
  this->optimiseSLU();
  for (size_t j = 0; j < cols(); ++j)
    for (int i = IA[j]; i < IA[j+1]; ++i)
      A[JA[i]][j] = SparseMatrix::A[i];

  return true;
}


bool ISTLMatrix::endAssembly()
{
  return true;
}


void ISTLMatrix::init ()
{
  SparseMatrix::init();

  // Set all matrix elements to zero
  A = 0;
}


/*! \brief Helper template for setting up a solver with the appropriate
           preconditioner type.
    \details We cannot instance using dynamic polymorphism, the solver need
             access to the real type for the preconditioner. We can however
             call the solver in the interface class scope afterwards.
 */
template<class Prec>
static Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*
  setupWithPreType(const LinSolParams& solParams,
                   Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>& op,
                   Dune::Preconditioner<ISTLMatrix::Vec,ISTLMatrix::Vec>& prec)
{
  Prec& pre = static_cast<Prec&>(prec);

  std::string type = solParams.getStringValue("type");
  double rtol = solParams.getDoubleValue("rtol");
  int maxits = solParams.getIntValue("maxits");
  int verbosity = solParams.getIntValue("verbosity");
  if (type == "bcgs")
    return new Dune::BiCGSTABSolver<ISTLMatrix::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "cg")
    return new Dune::CGSolver<ISTLMatrix::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "minres")
    return new Dune::MINRESSolver<ISTLMatrix::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "gmres") {
    int restart = solParams.getIntValue("gmres_restart_iterations");
    return new Dune::RestartedGMResSolver<ISTLMatrix::Vec>(op, pre, rtol,
                                                           restart, maxits, verbosity);
  } else
    std::cerr << "**ISTLMatrix** Unknown solver type " << type << "." << std::endl;

  return nullptr;
}

/*! \brief Helper template for setting up an AMG preconditioner with a given smoother */

template<class Smoother>
static std::pair<Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*,
                 Dune::Preconditioner<ISTLMatrix::Vec, ISTLMatrix::Vec>*> 
  setupAMG(const LinSolParams& params, size_t block,
           Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>& op, int nsd)
{
  // The coupling metric used in the AMG
  typedef Dune::Amg::FirstDiagonal CouplingMetric;
  // The coupling criterion used in the AMG
  typedef Dune::Amg::SymmetricCriterion<ISTLMatrix::Mat, CouplingMetric> CritBase;
  // The coarsening criterion used in the AMG
  typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;
  Criterion crit;
  typedef Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec> Operator;
  typedef typename Dune::Amg::AMG<Operator, ISTLMatrix::Vec, Smoother> AMG;
  typename AMG::SmootherArgs args;
  args.relaxationFactor = 1.0;
  args.iterations = std::max(1, params.getBlock(block).getIntValue("multigrid_no_smooth"));

  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  if (params.getBlock(block).hasValue("multigrid_no_smooth")) {
    int val = params.getBlock(block).getIntValue("multigrid_no_smooth");
    crit.setNoPreSmoothSteps(val);
    crit.setNoPostSmoothSteps(val);
  }

  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  crit.setDefaultValuesIsotropic(nsd);

  std::pair<Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*,
            Dune::Preconditioner<ISTLMatrix::Vec, ISTLMatrix::Vec>*> res;

  res.second = new AMG(op, crit, args);
  res.first = setupWithPreType<AMG>(params, op, *res.second);

  return res;
}


#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
/*! \brief Helper template for setting up a AMG preconditioner
 *!        with fine smoother differing from the other levels */

template<class FineSmoother, class Smoother>
static std::pair<Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*,
                 Dune::Preconditioner<ISTLMatrix::Vec, ISTLMatrix::Vec>*> 
  setupAMG2_full(const LinSolParams& params, size_t block,
                 Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>& op, int nsd,
                 FineSmoother* fsmooth)
{
  typedef Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec> Operator;
  typedef Dune::Amg::FirstDiagonal CouplingMetric;
  typedef Dune::Amg::SymmetricCriterion<ISTLMatrix::Mat,CouplingMetric> CritBase;
  typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;
  typedef Dune::Amg::AggregationLevelTransferPolicy<Operator,Criterion> TransferPolicy;
  typedef Dune::Amg::OneStepAMGCoarseSolverPolicy<Operator,Smoother,Criterion> CoarsePolicy;
  typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
  typedef typename Dune::Amg::TwoLevelMethod<Operator, CoarsePolicy, FineSmoother> AMG2;
  SmootherArgs args;
  args.relaxationFactor = 1.0;

  Criterion crit;
  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  if (params.getBlock(block).hasValue("multigrid_no_smooth")) {
    int val = params.getBlock(block).getIntValue("multigrid_no_smooth");
    crit.setNoPreSmoothSteps(val);
    crit.setNoPostSmoothSteps(val);
  }

  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  crit.setDefaultValuesIsotropic(nsd);

  CoarsePolicy coarsePolicy(args, crit);
  TransferPolicy policy(crit);
  std::pair<Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*,
            Dune::Preconditioner<ISTLMatrix::Vec, ISTLMatrix::Vec>*> res;
  Dune::shared_ptr<FineSmoother> fsp(fsmooth);
  res.second = new AMG2(op, fsp, policy, coarsePolicy);
  res.first = setupWithPreType<AMG2>(params, op, *res.second);

  return res;
}


template<class FineSmoother>
static std::pair<Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*,
                 Dune::Preconditioner<ISTLMatrix::Vec, ISTLMatrix::Vec>*> 
  setupAMG2_smoother(const LinSolParams& params, size_t block,
                     Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>& op, int nsd, FineSmoother* fsmooth)
{
  std::string smoother = params.getBlock(block).getStringValue("multigrid_smoother");
  if (smoother == "ssor")
    return setupAMG2_full<FineSmoother,Dune::SeqSSOR<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>>(params, block, op, nsd, fsmooth);
  else if (smoother == "sor")
    return setupAMG2_full<FineSmoother,Dune::SeqSOR<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>>(params, block, op, nsd, fsmooth);
  else if (smoother == "jacobi")
    return setupAMG2_full<FineSmoother,Dune::SeqJac<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>>(params, block, op, nsd, fsmooth);
  else if (smoother == "ilu") {
    int level = params.getBlock(block).getIntValue("ilu_fill_level");
    if (level > 0)
      std:: cerr << "**ISTLMATRIX ** ILU(n) smoothing currently not available, using ILU(0)." << std::endl;
    return setupAMG2_full<FineSmoother,Dune::SeqILU0<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>>(params, block, op, nsd, fsmooth);
  }
  else {
    std::cerr << "**ISTLMatrix** Invalid smoother " << smoother << "." << std::endl;
    return {nullptr, nullptr};
  }
}


static std::pair<Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*,
                 Dune::Preconditioner<ISTLMatrix::Vec, ISTLMatrix::Vec>*> 
  setupAMG2(const LinSolParams& params, size_t block,
            Dune::MatrixAdapter<ISTLMatrix::Mat,ISTLMatrix::Vec,ISTLMatrix::Vec>& op, int nsd)
{
  std::string fsmoother = params.getBlock(block).getStringValue("multigrid_finesmoother");
  int nosmooth = std::max(1, params.getBlock(block).getIntValue("multigrid_no_fine_smooth"));
  if (fsmoother == "ilu") {
    auto fsmooth = new Dune::SeqILU0<ISTLMatrix::Mat, ISTLMatrix::Vec, ISTLMatrix::Vec>(op.getmat(), 1.0);
    return setupAMG2_smoother(params, block, op, nsd, fsmooth);
  } else if (fsmoother == "ssor") {
    auto fsmooth = new Dune::SeqSSOR<ISTLMatrix::Mat, ISTLMatrix::Vec, ISTLMatrix::Vec>(op.getmat(), nosmooth, 1.0);
    return setupAMG2_smoother(params, block, op, nsd, fsmooth);
  } else {
    std::cerr << "**ISTLMatrix** Invalid fine smoother " << fsmoother << "." << std::endl;
    return {nullptr, nullptr};
  }
}
#endif


void ISTLMatrix::setupSolver()
{
  if (!pre) {
    std::string prec = solParams.getBlock(0).getStringValue("pc");
    if (prec == "ilu") {
      if (solParams.getBlock(0).getIntValue("ilu_fill_level") == 0) {
        pre.reset(new Dune::SeqILU0<Mat,Vec,Vec>(A, 1.0));
        solver.reset(setupWithPreType<Dune::SeqILU0<Mat,Vec,Vec>>(solParams, op, *pre));
      } else {
        pre.reset(new Dune::SeqILUn<Mat,Vec,Vec>(A, solParams.getBlock(0).getIntValue("ilu_fill_level"), 1.0));
        solver.reset(setupWithPreType<Dune::SeqILUn<Mat,Vec,Vec>>(solParams, op, *pre));
      }
    } else if (prec == "sor") {
      pre.reset(new Dune::SeqSOR<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqSOR<Mat,Vec,Vec>>(solParams, op, *pre));
    } else if (prec == "ssor") {
      pre.reset(new Dune::SeqSOR<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqSSOR<Mat,Vec,Vec>>(solParams, op, *pre));
    } else if (prec == "jacobi") {
      pre.reset(new Dune::SeqJac<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqJac<Mat,Vec,Vec>>(solParams, op, *pre));
    } else if (prec == "gs") {
      pre.reset(new Dune::SeqGS<Mat,Vec,Vec>(A, 1, 1.0));
      solver.reset(setupWithPreType<Dune::SeqGS<Mat,Vec,Vec>>(solParams, op, *pre));
    }
    else if (prec == "asm" || prec == "asmlu") {
      size_t nx = solParams.getBlock(0).getIntValue("nx");
      nx = std::max(1ul, nx);
      size_t ny = solParams.getBlock(0).getIntValue("ny");
      ny = std::max(1ul, ny);
      size_t nz = solParams.getBlock(0).getIntValue("nz");
      nz = std::max(1ul, nz);
      int overlap = solParams.getBlock(0).getIntValue("overlap");

      if (!samp)
        return;

      std::vector<std::set<int>>  locSubdDofs(nx*ny*nz*samp->getPatches().size());
      size_t d = 0;
      for (const auto& it : samp->getPatches()) {
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
        Dune::SeqOverlappingSchwarz<Mat, Vec, Dune::AdditiveSchwarzMode, Dune::SuperLU<Mat>>::subdomain_vector ddofs(locSubdDofs.size());
        for (size_t i = 0; i < locSubdDofs.size(); ++i)
          for (const auto& it : locSubdDofs[i])
            ddofs[i].insert(it);

        pre.reset(new Dune::SeqOverlappingSchwarz<Mat, Vec, Dune::AdditiveSchwarzMode, Dune::SuperLU<Mat>>(A, ddofs));
        solver.reset(setupWithPreType<Dune::SeqOverlappingSchwarz<Mat, Vec, Dune::AdditiveSchwarzMode, Dune::SuperLU<Mat>>>(solParams, op, *pre));
      } else {
        Dune::SeqOverlappingSchwarz<Mat, Vec>::subdomain_vector ddofs(locSubdDofs.size());
        for (size_t i = 0; i < locSubdDofs.size(); ++i)
          ddofs[i].insert(locSubdDofs[i].begin(), locSubdDofs[i].end());

        for (auto& it : ddofs) {
          for (auto& it2 : it)
            std::cout << it2 << " ";
          std::cout << std::endl;
        }

        pre.reset(new Dune::SeqOverlappingSchwarz<Mat, Vec>(A, ddofs));
        solver.reset(setupWithPreType<Dune::SeqOverlappingSchwarz<Mat, Vec>>(solParams, op, *pre));
      }
    } else if (prec == "amg") {
      std::pair<Dune::InverseOperator<ISTLMatrix::Vec,ISTLMatrix::Vec>*,
                Dune::Preconditioner<ISTLMatrix::Vec, ISTLMatrix::Vec>*> res {nullptr, nullptr};
      if (solParams.getBlock(0).getStringValue("multigrid_smoother") !=
          solParams.getBlock(0).getStringValue("multigrid_finesmoother")) {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
        auto res = setupAMG2(solParams, 0, op, adm.dd.getNoSpaceDim());
        solver.reset(res.first);
        pre.reset(res.second);
#else
        std::cerr << "**ISTLMatrix** Separate fine smoother not implemented." << std::endl;
        return;
#endif
      } else {
        std::string smoother = solParams.getBlock(0).getStringValue("multigrid_smoother");
        if (smoother.empty()) {
          std::cerr << "** ISTLMatrix ** No smoother defined for AMG, defaulting to ILU" << std::endl;
          smoother = "ilu";
        }
        if (smoother == "ilu")
          res = setupAMG<Dune::SeqILU0<Mat,Vec,Vec>>(solParams, 0, op, adm.dd.getNoSpaceDim());
        else if (smoother == "sor")
          res = setupAMG<Dune::SeqSOR<Mat,Vec,Vec>>(solParams, 0, op, adm.dd.getNoSpaceDim());
        else if (smoother == "ssor")
          res = setupAMG<Dune::SeqSSOR<Mat,Vec,Vec>>(solParams, 0, op, adm.dd.getNoSpaceDim());
        else if (smoother == "jacobi")
          res = setupAMG<Dune::SeqJac<Mat,Vec,Vec>>(solParams, 0, op, adm.dd.getNoSpaceDim());

        solver.reset(res.first);
        pre.reset(res.second);
      }
    }
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
    else if (prec == "fastamg") {
      // The coupling metric used in the AMG
      typedef Dune::Amg::FirstDiagonal CouplingMetric;
      // The coupling criterion used in the AMG
      typedef Dune::Amg::SymmetricCriterion<ISTLMatrix::Mat, CouplingMetric> CritBase;
      // The coarsening criterion used in the AMG
      typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;
      Criterion crit;

      pre.reset(new Dune::Amg::FastAMG<Operator,Vec>(op, crit));
      solver.reset(setupWithPreType<Dune::Amg::FastAMG<Operator,Vec>>(solParams, op, *pre));
    }
#endif
    else
      std::cerr << "**ISTLMatrix** Unknown preconditioner " << prec << "." << std::endl;
  }
}


bool ISTLMatrix::solve (SystemVector& B, bool newLHS, Real*)
{
  setupSolver();

  ISTLVector* Bptr = dynamic_cast<ISTLVector*>(&B);
  if (!Bptr || !solver)
    return false;

  try {
    Dune::InverseOperatorResult r;
    Vec b(Bptr->getVector());
    Bptr->getVector() = 0;
    solver->apply(Bptr->getVector(), b, r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Bptr)(i+1) = Bptr->getVector()[i];

  return true;
}


bool ISTLMatrix::solve (const SystemVector& b, SystemVector& x, bool newLHS)
{
  setupSolver();

  const ISTLVector* Bptr = dynamic_cast<const ISTLVector*>(&b);
  if (!Bptr || ! solver)
    return false;

  ISTLVector* Xptr = dynamic_cast<ISTLVector*>(&x);
  if (!Xptr)
    return false;

  try {
    Dune::InverseOperatorResult r;
    solver->apply(Xptr->getVector(),
                  const_cast<Vec&>(Bptr->getVector()), r);
  } catch (Dune::ISTLError& e) {
    std::cerr << "ISTL exception " << e << std::endl;
    return false;
  }

  for (size_t i = 0; i < rows(); ++i)
    (*Xptr)(i+1) = Xptr->getVector()[i];

  return true;
}


Real ISTLMatrix::Linfnorm () const
{
  return A.infinity_norm();
}
