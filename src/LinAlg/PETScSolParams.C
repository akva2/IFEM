// $Id$
//==============================================================================
//!
//! \file PETScSolParams.C
//!
//! \date Mar 10 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//!
//==============================================================================

#include "PETScSolParams.h"
#include "PCPerm.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <fstream>
#include <sstream>
#include <utility>
#include <iterator>


void PETScSolParams::setParams(const LinSolParams& params,
                               KSP& ksp, PetscIntMat& locSubdDofs,
			       PetscIntMat& subdDofs, PetscRealVec& coords,
			       ISMat& dirIndexSet, int nsd)
{
  // Set linear solver method
  KSPSetType(ksp,params.getStringValue("type").c_str());
  KSPSetTolerances(ksp,params.getDoubleValue("rtol"),
                   params.getDoubleValue("atol"),
                   params.getDoubleValue("dtol"),
                   params.getIntValue("maxits"));

  // Set preconditioner
  PC pc;
  KSPGetPC(ksp,&pc);
  std::string prec = params.getStringValue("pc");
  if (!strncasecmp(prec.c_str(),"compositedir",12)) {
    Mat mat;
    Mat Pmat;
#if PETSC_VERSION_MINOR < 5
    MatStructure flag;
    PCGetOperators(pc,&mat,&Pmat,&flag);
#else
    PCGetOperators(pc,&mat,&Pmat);
#endif
    addDirSmoother(pc,Pmat,params.getBlock(0),0,dirIndexSet);
  }
  else
    PCSetType(pc,prec.c_str());

#if PETSC_HAVE_HYPRE
  if (!strncasecmp(prec.c_str(),"hypre",5)) {
    PCHYPRESetType(pc,params.getBlock(0).getStringValue("hypre_type").c_str());
    setHypreOptions("", params.getBlock(0));
  }
#endif

  if (!strncasecmp(prec.c_str(),"gamg",4) || !strncasecmp(prec.c_str(),"ml",2) ) {
    PetscInt nloc = coords.size()/nsd;
    PCSetCoordinates(pc,nsd,nloc,&coords[0]);
    PCGAMGSetType(pc, PCGAMGAGG); // TODO?
  }

  if (!strncasecmp(prec.c_str(),"asm",3) || !strncasecmp(prec.c_str(),"gasm",4))
    setupAdditiveSchwarz(pc, params.getBlock(0).getIntValue("overlap"), !strncasecmp(prec.c_str(),"asmlu",5),
                         locSubdDofs, subdDofs, false);
  else if (!strncasecmp(prec.c_str(),"ml",2) || !strncasecmp(prec.c_str(),"gamg",4)) {
    if (!strncasecmp(prec.c_str(),"ml",2))
      setMLOptions("", params.getBlock(0));
    else if (!strncasecmp(prec.c_str(),"gamg",4))
      setGAMGOptions("", params.getBlock(0));

    PCSetFromOptions(pc);
    PCSetUp(pc);

    // Settings for coarse solver
    if (!params.getBlock(0).hasValue("ml_coarse_solver"))
      setupCoarseSolver(pc, "", params.getBlock(0));

    setupSmoothers(pc, params, 0, dirIndexSet, locSubdDofs, subdDofs);
  }
  else {
    PCSetFromOptions(pc);
    PCSetUp(pc);
  }

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  MPI_Comm comm;
  PetscObjectGetComm((PetscObject)ksp,&comm);
  PCView(pc,PETSC_VIEWER_STDOUT_(comm));
}


bool PETScSolParams::addDirSmoother(PC pc, Mat P, const LinSolParams::BlockParams& block,
                                    int iBlock, ISMat& dirIndexSet)
{
  PCSetType(pc,"composite");
  PCCompositeSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);
  for (size_t k = 0;k < block.dirSmoother.size();k++)
    PCCompositeAddPC(pc,"shell");
  for (size_t k = 0;k < block.dirSmoother.size();k++) {
    PC dirpc;
    PCCompositeGetPC(pc,k,&dirpc);
    PCPerm *pcperm;
    PCPermCreate(&pcperm);
    PCShellSetApply(dirpc,PCPermApply);
    PCShellSetContext(dirpc,pcperm);
    PCShellSetDestroy(dirpc,PCPermDestroy);
    PCShellSetName(dirpc,"dir");
    PCPermSetUp(dirpc,&dirIndexSet[iBlock][k],P,block.dirSmoother[k].type.c_str());
  }

  return true;
}


/*! \brief Static helper to optionally add a prefix to a PETSc parameter */
static std::string AddPrefix(const std::string& prefix, const std::string& data)
{
  if (prefix.empty())
    return "-"+data;

  return "-"+prefix+"_"+data;
}


static void condSetup(const std::string& prefix, const std::string& petsc_option,
                      const std::string& map_option, const SettingMap& map)
{
  if (map.hasValue(map_option))
    PetscOptionsSetValue(AddPrefix(prefix,petsc_option).c_str(), map.getStringValue(map_option).c_str());
}


void PETScSolParams::setMLOptions(const std::string& prefix, const SettingMap& map)
{
  condSetup(prefix, "pc_ml_maxNLevels", "multigrid_levels", map);
  condSetup(prefix, "pc_ml_maxCoarseSize", "ml_max_coarse_size", map);
  condSetup(prefix, "pc_ml_maxCoarsenScheme", "ml_coarsen_scheme", map);
  condSetup(prefix, "pc_ml_Threshold", "ml_threshold", map);
  condSetup(prefix, "pc_ml_DampingFactor", "ml_damping_factor", map);
  condSetup(prefix, "pc_ml_repartitionMaxMinRatio", "ml_repartition_max_min_ratio", map);
  condSetup(prefix, "pc_ml_Symmetrize", "ml_symmetrize", map);
  condSetup(prefix, "pc_ml_repartition", "ml_repartition", map);
  condSetup(prefix, "pc_ml_BlockScaling", "ml_block_scaling", map);
  condSetup(prefix, "pc_ml_repartitionPutOnSingleProc", "ml_put_on_single_proc", map);
  condSetup(prefix, "pc_ml_reuse_interpolation", "ml_reuse_interpolation", map);
  condSetup(prefix, "pc_ml_KeepAggInfo", "ml_keep_agg_info", map);
  condSetup(prefix, "pc_ml_Reusable", "ml_reusable", map);
  condSetup(prefix, "pc_ml_Aux", "ml_aux", map);
  condSetup(prefix, "pc_ml_AuxThreshold", "ml_aux_threshold", map);
}


void PETScSolParams::setGAMGOptions(const std::string& prefix, const SettingMap& map)
{
  condSetup(prefix, "pc_gamg_type", "gamg_type", map);
  condSetup(prefix, "pc_gamg_coarse_eq_limit", "gamg_coarse_eq_limit", map);
  condSetup(prefix, "pc_gamg_process_eq_limit", "gamg_process_eq_limit", map);
  condSetup(prefix, "pc_gamg_repartition", "gamg_repartition", map);
  condSetup(prefix, "pc_gamg_use_agg_gasm", "gamg_use_agg_gasm", map);
  condSetup(prefix, "pc_gamg_reuse_interpolation", "gamg_reuse_interpolation", map);
  condSetup(prefix, "pc_gamg_threshold", "gamg_threshold", map);
  condSetup(prefix, "pc_mg_levels", "multigrid_levels", map);
}


void PETScSolParams::setHypreOptions(const std::string& prefix, const SettingMap& map)
{
  condSetup(prefix,"pc_hypre_type", "hypre_type", map);
  condSetup(prefix, "pc_hypre_boomeramg_max_levels", "multigrid_levels", map);
  // TODO: Investigate why these are hardcoded.
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_max_iter").c_str(), "1");
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_tol").c_str(), "0.0");;
  condSetup(prefix, "pc_hypre_boomeramg_strong_threshold", "hypre_threshold", map);
  condSetup(prefix, "pc_hypre_boomeramg_coarsen_type", "hypre_boomeramg_coarsen_type", map);
  condSetup(prefix, "pc_hypre_boomeramg_agg_nl", "hypre_boomeramg_agg_nl", map);
  condSetup(prefix, "pc_hypre_boomeramg_agg_num_paths", "hypre_boomeramg_agg_num_paths", map);
  condSetup(prefix, "pc_hypre_boomeramg_truncfactor", "hypre_boomeramg_truncation_factor", map);
}


void PETScSolParams::setupCoarseSolver(PC& pc, const std::string& prefix, const SettingMap& map)
{
  std::string coarseSolver = map.getStringValue("ml_coarse_solver");
  std::string coarsePackage = map.getStringValue("ml_coarse_package");
  if (coarseSolver == "OneLevelSchwarz" ||
      coarseSolver == "TwoLevelSchwarz") {
    KSP cksp;
    PC  cpc;
    PCMGGetCoarseSolve(pc,&cksp);
    KSPSetType(cksp,"preonly");
    KSPGetPC(cksp,&cpc);
    PCSetType(cpc,"redistribute");
    PCSetUp(cpc);

    KSP sksp;
    PC  spc;
    PCRedistributeGetKSP(cpc,&sksp);
    KSPSetTolerances(sksp,1.0e-2,PETSC_DEFAULT,PETSC_DEFAULT,10);
    KSPGetPC(sksp,&spc);
    if (coarseSolver == "OneLevelSchwarz") {
      PCSetType(spc,PCGASM);
      PCSetUp(spc);
    } else {
      PCSetType(spc,PCML);
      PetscOptionsSetValue(AddPrefix(prefix,"mg_coarse_pc_ml_maxNLevels").c_str(),"2");	
      PCSetFromOptions(spc);
      PCSetUp(spc);

      KSP csksp;
      PC cspc;
      PCMGGetSmoother(spc,1,&csksp);
      KSPSetType(csksp,"richardson");
      KSPGetPC(csksp,&cspc);
      PCSetType(cspc,"asm");
      PCSetUp(cspc);
    }

    KSP* subsksp;
    PC   subspc;
    PetscInt first, nlocal;
    PCGASMGetSubKSP(spc,&nlocal,&first,&subsksp);
    for (int k = 0; k < nlocal; k++) {
      KSPGetPC(subsksp[k],&subspc);
      PCSetType(subspc,PCLU);
      KSPSetType(subsksp[k],KSPPREONLY);
    }
  }

  if (!coarsePackage.empty()) {
    KSP cksp;
    PC  cpc;
    PCMGGetCoarseSolve(pc,&cksp);
    KSPGetPC(cksp,&cpc);
    PCSetType(cpc,PCLU);
    PCFactorSetMatSolverPackage(cpc,coarsePackage.c_str());
    PCSetUp(cpc);
  }
}


void PETScSolParams::setupSmoothers(PC& pc, const LinSolParams& params, size_t iBlock,
                                    ISMat& dirIndexSet,
                                    const PetscIntMat& locSubdDofs,
                                    const PetscIntMat& subdDofs)
{
  PetscInt n;
  PCMGGetLevels(pc,&n);

  std::string mgKSP = params.getBlock(iBlock).getStringValue("multigrid_ksp");

  // Presmoother settings
  for (int i = 1;i < n;i++) {
    KSP preksp;
    PC  prepc;

    // Set smoother
    std::string smoother;
    PetscInt noSmooth;

    PCMGGetSmoother(pc,i,&preksp);

    // warn that richardson might break symmetry if the KSP is CG
    if (mgKSP == "defrichardson" && params.getStringValue("type") == KSPCG)
      std::cerr << "WARNING: Using multigrid with Richardson on sublevels.\n"
                << "If you get divergence with KSP_DIVERGED_INDEFINITE_PC, try\n"
                << "adding <mgksp>chebyshev</mgksp. Add <mgksp>richardson</mgksp>\n"
                << "to quell this warning." << std::endl;

    if (mgKSP == "richardson" || mgKSP == "defrichardson")
      KSPSetType(preksp,KSPRICHARDSON);
    else if (mgKSP == "chebyshev")
      KSPSetType(preksp,KSPCHEBYSHEV);

    std::string finesmoother = params.getBlock(iBlock).getStringValue("multigrid_finesmoother");
    if ((i == n-1) && !finesmoother.empty()) {
      smoother = finesmoother;
      noSmooth = params.getBlock(iBlock).getIntValue("multigrid_no_fine_smooth");
    }
    else {
      smoother = params.getBlock(iBlock).getStringValue("multigrid_smoother");
      noSmooth = params.getBlock(iBlock).getIntValue("multigrid_no_smooth");;
    }

    KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noSmooth);
    KSPGetPC(preksp,&prepc);

    if (smoother == "asm" || smoother == "asmlu")
      setupAdditiveSchwarz(prepc, params.getBlock(iBlock).getIntValue("overlap"), smoother == "asmlu",
                           locSubdDofs, subdDofs, true);
    else if (smoother == "compositedir" && (i==n-1)) {
      Mat mat;
      Mat Pmat;
#if PETSC_VERSION_MINOR < 5
      MatStructure flag;
      PCGetOperators(prepc,&mat,&Pmat,&flag);
#else
      PCGetOperators(prepc,&mat,&Pmat);
#endif

      addDirSmoother(prepc,Pmat,params.getBlock(iBlock),iBlock,dirIndexSet);
    }
    else
      PCSetType(prepc,smoother.c_str());

    PCFactorSetLevels(prepc,params.getBlock(0).getIntValue("ilu_fill_level"));
    KSPSetUp(preksp);
  }
}


void PETScSolParams::setupAdditiveSchwarz(PC& pc, int overlap, bool asmlu,
                                          const PetscIntMat& locSubdDofs,
                                          const PetscIntMat& subdDofs, bool smoother)
{
  PCSetType(pc, PCASM);
  if (!smoother)
    PCASMSetType(pc,PC_ASM_BASIC);
  PCASMSetOverlap(pc,overlap);

  if (!locSubdDofs.empty() && !subdDofs.empty()) {
    const size_t nsubds = subdDofs.size();

    IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
    for (size_t i = 0;i < nsubds;i++) {
      ISCreateGeneral(PETSC_COMM_SELF,locSubdDofs[i].size(),
                      &(const_cast<PetscIntMat&>(locSubdDofs)[i][0]),
                      PETSC_USE_POINTER,&(isLocSubdDofs[i]));
      ISCreateGeneral(PETSC_COMM_SELF,subdDofs[i].size(),
                      &(const_cast<PetscIntMat&>(subdDofs)[i][0]),
                      PETSC_USE_POINTER,&(isSubdDofs[i]));
    }
    PCASMSetLocalSubdomains(pc,nsubds,isSubdDofs,isLocSubdDofs);
  }

  PCSetFromOptions(pc);
  PCSetUp(pc);

  if (asmlu) {
    KSP* subksp;
    PC   subpc;
    PetscInt first, nlocal;
    PCASMGetSubKSP(pc,&nlocal,&first,&subksp);

    for (int i = 0; i < nlocal; i++) {
      KSPGetPC(subksp[i],&subpc);
      PCSetType(subpc,PCLU);
      KSPSetType(subksp[i],KSPPREONLY);
    }
  }
}
