// $Id$
//==============================================================================
//!
//! \file PETScSchurPC.C
//!
//! \date Jun 4 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Schur-complement preconditioner using PETSc.
//!
//==============================================================================

#include "PETScSchurPC.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"

#include <iostream>


PETScSchurPC::PETScSchurPC (PC& pc_init, const std::vector<Mat>& blocks,
                            const LinSolParams::BlockParams& params, const ProcessAdm& adm,
                            const PETScMGLevels* mg, int verbosity)
  : m_blocks(&blocks)
{
  PCSetType(pc_init, PCSHELL);
  PCShellSetContext(pc_init, this);
  PCShellSetName(pc_init, "Schur complement preconditioner");
  PCShellSetApply(pc_init, PETScSchurPC::Apply_Outer);
  PCShellSetDestroy(pc_init, PETScSchurPC::Destroy);

  KSPCreate(*adm.getCommunicator(), &inner_ksp);
  KSPSetType(inner_ksp, KSPPREONLY);
  KSPSetOperators(inner_ksp, blocks[0], blocks[0]);
  PC pc;
  KSPGetPC(inner_ksp, &pc);
  // The <schur> container of the input gives its settings the prefix schur_,
  // so the preconditioner of the inner solve arrives as schur_pc. The
  // unprefixed spelling is honoured as well, since that is the key this used
  // to look for, although no input file in the tree sets it.
  std::string schurpc = params.getStringValue("schur_pc");
  if (schurpc.empty())
    schurpc = params.getStringValue("schurpc");
  if (schurpc.empty())
    schurpc = PCGAMG;

  // The inner solve is where the approximation of the inverse of the momentum
  // operator enters the Schur complement, and it is the expensive part of the
  // preconditioner. A geometric multigrid hierarchy built from the adaptive
  // mesh sequence goes in here in place of the algebraic default.
  if (schurpc == "gmg") {
    if (mg && mg->size() > 1) {
      LinSolParams dummy;
      PETScSolParams(dummy,adm).setupGeometricMG(pc,*mg,params);
    }
    else {
      std::cerr <<"  ** PETScSchurPC: No geometric multigrid hierarchy for the"
                <<" inner operator,\n     falling back on "<< PCGAMG <<"."
                << std::endl;
      PCSetType(pc,PCGAMG);
    }
  }
  else
    PCSetType(pc, schurpc.c_str());

  KSPSetFromOptions(inner_ksp);
  KSPSetUp(inner_ksp);
  if (verbosity > 1)
    KSPView(inner_ksp, PETSC_VIEWER_STDOUT_WORLD);

  KSPCreate(*adm.getCommunicator(), &outer_ksp);
  KSPGetPC(outer_ksp, &pc);
  PCSetType(pc, PCNONE);
  int maxits = params.getIntValue("schur_maxits");
  if (maxits < 1)
    maxits = 1000;
  double atol = params.getDoubleValue("schur_atol");
  if (atol == 0.0)
    atol = 1e-16;
  double dtol = params.getDoubleValue("schur_dtol");
  if (dtol == 0.0)
    dtol = 1e100;
  double rtol = params.getDoubleValue("schur_rtol");
  if (rtol == 0.0)
    rtol = 1e-12;
  std::string type = params.getStringValue("schur_type");
  if (type.empty())
    type = KSPGMRES;

  // The outer solver has no preconditioner, so preonly makes it the identity
  // and the Schur complement operator is never applied at all, inner solve
  // included. That is rarely what is wanted, and is silent otherwise.
  if (type == KSPPREONLY)
    std::cerr <<"  ** PETScSchurPC: The Schur complement solver is preonly"
              <<" and has no preconditioner,\n     so it reduces to the"
              <<" identity and the Schur operator is never applied."
              << std::endl;

  KSPSetType(outer_ksp, type.c_str());
  KSPSetTolerances(outer_ksp, rtol, atol, dtol, maxits);

  MatCreate(*adm.getCommunicator(), &outer_mat);
  PetscInt r;
  MatGetSize(blocks[3], &r, &r);
  MatSetSizes(outer_mat, r, r, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetFromOptions(outer_mat);
  MatSetType(outer_mat, MATSHELL);
  MatShellSetContext(outer_mat, this);
  MatShellSetOperation(outer_mat, MATOP_MULT,
                       (void(*)(void))&PETScSchurPC::Apply_Schur);
  MatSetUp(outer_mat);

  KSPSetOperators(outer_ksp, outer_mat, outer_mat);
  KSPSetFromOptions(outer_ksp);
  KSPSetUp(outer_ksp);
  if (verbosity > 1)
    KSPView(outer_ksp, PETSC_VIEWER_STDOUT_WORLD);
  PCSetUp(pc_init);

  VecCreate(*adm.getCommunicator(), &tmp);
  VecSetFromOptions(tmp);
  MatGetSize(blocks[0], &r, &r);
  VecSetSizes(tmp, r, PETSC_DETERMINE);
  VecSetUp(tmp);
}


PETScSchurPC::~PETScSchurPC ()
{
  KSPDestroy(&inner_ksp);
  KSPDestroy(&outer_ksp);
  MatDestroy(&outer_mat);
  VecDestroy(&tmp);
}


PetscErrorCode PETScSchurPC::Apply_Schur (Mat A, Vec x, Vec y)
{
  void* p;
  MatShellGetContext(A, &p);
  PETScSchurPC* spc = static_cast<PETScSchurPC*>(p);
  MatMult(spc->m_blocks->at(1), x, spc->tmp);
  KSPSolve(spc->inner_ksp, spc->tmp, spc->tmp);
  MatMult(spc->m_blocks->at(2), spc->tmp, y);

  return 0;
}


PetscErrorCode PETScSchurPC::Apply_Outer (PC pc, Vec x, Vec y)
{
  void* p;
  PCShellGetContext(pc, &p);
  PETScSchurPC* spc = static_cast<PETScSchurPC*>(p);
  KSPSolve(spc->outer_ksp, x, y);

  return 0;
}


PetscErrorCode PETScSchurPC::Destroy (PC pc)
{
  PETScSchurPC *shell;

  PCShellGetContext(pc,(void**)&shell);
  delete shell;

  return 0;
}
