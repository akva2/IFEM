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

#include <cstring>
#include <iostream>


PETScSchurPC::PETScSchurPC (PC& pc_init, const std::vector<Mat>& blocks,
                            const LinSolParams::BlockParams& params, const ProcessAdm& adm,
                            const PETScMGLevels* mg, int verbosity,
                            const SchurOperators& ops)
  : pressure(ops),
    m_blocks(&blocks),
    hasC(params.hasValue("has_C") ? params.getIntValue("has_C") > 0
                                  : ops.hasC)
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

  // Applying the Schur complement means a Krylov solve against it, and what
  // that costs decides whether the preconditioner is worth having at all.
  // Unpreconditioned it takes hundreds of iterations, every one of them an
  // inner solve of the momentum operator, so it is cut short long before it
  // converges and the outer solver is handed something inconsistent.
  //
  // Whoever set this preconditioner up may have left an approximation of the
  // Schur complement on it. It is a matrix over the same equations, and
  // whatever it is worth as an approximation it is worth more as something
  // to precondition the solve against the real thing with.
  Mat sA = nullptr, sP = nullptr;
  MatType ptype = nullptr;
  PCGetOperators(pc_init, &sA, &sP);
  if (sP)
    MatGetType(sP, &ptype);

  // Only an assembled one is of any use: where the Schur complement is
  // preconditioned with itself, what comes back is the operator again and
  // there is nothing in it to factorize.
  const bool preconditioned = !pressure.empty() ||
                              (sP && ptype &&
                               strcmp(ptype, MATSHELL) &&
                               strcmp(ptype, MATSCHURCOMPLEMENT));

  // What is applied above is the negative of the Schur complement, that
  // being the positive definite one of the two, while the approximation
  // left here is of the Schur complement itself. Preconditioning the one
  // with the other takes a change of sign, which an operator with a sign
  // to speak of does not survive being without.
  if (!pressure.empty())
  {
    // The pressure operators are what the Schur complement is like, rather
    // than an approximation of what it is, so they are applied as a
    // preconditioner of their own rather than factorized.
    PCSetType(pc, PCSHELL);
    PCShellSetContext(pc, this);
    PCShellSetName(pc, "Cahouet-Chabard");
    PCShellSetApply(pc, PETScSchurPC::Apply_CahouetChabard);

    auto&& solverFor = [&adm](Mat A, KSP& ksp)
    {
      if (!A) return;
      KSPCreate(*adm.getCommunicator(), &ksp);
      KSPSetOperators(ksp, A, A);
      KSPSetType(ksp, KSPPREONLY);
      PC sub;
      KSPGetPC(ksp, &sub);
      PCSetType(sub, PCILU);
      KSPSetUp(ksp);
    };
    solverFor(pressure.mass, massKsp);
    solverFor(pressure.laplacian, lapKsp);
    MatCreateVecs(pressure.mass ? pressure.mass : pressure.laplacian,
                  nullptr, &ctmp);
  }
  else if (preconditioned)
  {
    MatDuplicate(sP, MAT_COPY_VALUES, &prec);
    MatScale(prec, -1.0);
    PCSetType(pc, PCILU);
    // A constraint on the integrated pressure leaves its multiplier in these
    // equations, and the row it sits in has nothing on the diagonal for a
    // factorization to work with.
    PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
  }
  else
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

  // The Schur complement operator is laid out over the processes the way the
  // block it acts on is, so it is the local sizes of that block which say how
  // much of it belongs here. Taking the global ones instead makes every
  // process claim the whole of it, which is the same thing on one process
  // and a matrix as many times too large as there are processes on more.
  MatCreate(*adm.getCommunicator(), &outer_mat);
  PetscInt m, n;
  MatGetLocalSize(blocks[3], &m, &n);
  MatSetSizes(outer_mat, m, n, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetFromOptions(outer_mat);
  MatSetType(outer_mat, MATSHELL);
  MatShellSetContext(outer_mat, this);
  MatShellSetOperation(outer_mat, MATOP_MULT,
                       (void(*)(void))&PETScSchurPC::Apply_Schur);
  MatSetUp(outer_mat);

  KSPSetOperators(outer_ksp, outer_mat,
                  prec ? prec : outer_mat);
  KSPSetFromOptions(outer_ksp);
  KSPSetUp(outer_ksp);
  if (verbosity > 1)
    KSPView(outer_ksp, PETSC_VIEWER_STDOUT_WORLD);
  PCSetUp(pc_init);

  // The inner solve works on the rows of the momentum operator, and the
  // pressure-pressure term on those of the block being preconditioned, so
  // the two vectors take their layout from those rather than being sized by
  // hand.
  MatCreateVecs(blocks[0], nullptr, &tmp);
  if (hasC)
    MatCreateVecs(blocks[3], nullptr, &ptmp);
}


PETScSchurPC::~PETScSchurPC ()
{
  KSPDestroy(&inner_ksp);
  KSPDestroy(&outer_ksp);
  MatDestroy(&outer_mat);
  if (prec) MatDestroy(&prec);
  if (massKsp) KSPDestroy(&massKsp);
  if (lapKsp) KSPDestroy(&lapKsp);
  if (ctmp) VecDestroy(&ctmp);
  VecDestroy(&tmp);
  if (ptmp) VecDestroy(&ptmp);
}


/*!
  The Schur complement of the block being preconditioned is

  \f[ {\bf S} = {\bf A}_{11} - {\bf A}_{10}{\bf A}_{00}^{-1}{\bf A}_{01} \f]

  and what is applied here is its negative, which is the positive definite
  one of the two and so the one an outer solver asking for that can be given.

  The first of the two terms is nothing for a mixed discretization, which has
  no pressure-pressure coupling, and is the stabilization for one which is
  stabilized. Leaving it out costs the latter the very term which makes its
  system solvable, and the preconditioner is then not one.

  It is not enough to look at whether that block holds anything, since a
  constraint on the integrated pressure puts its multiplier there and a mixed
  discretization is then no longer empty in the corner without coupling the
  pressure to itself in any way that belongs in a Schur complement. Which of
  the two it is, is asked of the discretization.
*/

PetscErrorCode PETScSchurPC::Apply_Schur (Mat A, Vec x, Vec y)
{
  void* p;
  MatShellGetContext(A, &p);
  PETScSchurPC* spc = static_cast<PETScSchurPC*>(p);
  MatMult(spc->m_blocks->at(1), x, spc->tmp);
  KSPSolve(spc->inner_ksp, spc->tmp, spc->tmp);
  MatMult(spc->m_blocks->at(2), spc->tmp, y);

  if (spc->hasC)
  {
    MatMult(spc->m_blocks->at(3), x, spc->ptmp);
    VecAXPY(y, -1.0, spc->ptmp);
  }

  return 0;
}


/*!
  The Schur complement of a generalised Stokes system is like the pressure
  mass matrix where the problem is steady and like the pressure Laplacian
  where the time increment is what dominates it, and

  \f[ {\bf S}^{-1} \approx \mu{\bf M}_p^{-1} + \alpha{\bf L}_p^{-1} \f]

  holds across the two, which is the preconditioner of Cahouet and Chabard.
  It is a sum of the two inverses and not the inverse of a sum: as the time
  increment grows the second term falls away and leaves the first, which a
  sum of the operators themselves would not do.
*/

PetscErrorCode PETScSchurPC::Apply_CahouetChabard (PC pc, Vec x, Vec y)
{
  void* p;
  PCShellGetContext(pc, &p);
  PETScSchurPC* spc = static_cast<PETScSchurPC*>(p);

  VecZeroEntries(y);

  if (spc->massKsp && spc->pressure.viscosity != 0.0)
  {
    KSPSolve(spc->massKsp, x, spc->ctmp);
    VecAXPY(y, spc->pressure.viscosity, spc->ctmp);
  }

  if (spc->lapKsp && spc->pressure.transient != 0.0)
  {
    KSPSolve(spc->lapKsp, x, spc->ctmp);
    VecAXPY(y, spc->pressure.transient, spc->ctmp);
  }

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
