//==============================================================================
//!
//! \file OptControlOps.C
//!
//! \date Sep 19 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Reduced-space operators for linear-quadratic optimal control.
//!
//==============================================================================

#include "OptControlOps.h"

#include "SAM.h"
#include "SIMbase.h"
#include "SIMoptions.h"
#include "Profiler.h"
#include "SparseMatrix.h"
#include "SystemMatrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>


OptControlOps::OptControlOps (SIMbase& sim) : mySim(sim)
{
}


OptControlOps::~OptControlOps () = default;


bool OptControlOps::initSystem ()
{
  return mySim.initSystem(mySim.opt.solver,1,1);
}


bool OptControlOps::setAdjointTranspose (bool yes)
{
  adjTranspose = yes;
  if (!yes)
    return true;

  SystemMatrix* A = mySim.getLHSmatrix();
  if (A && A->canSolveTranspose())
    return true;

  std::cerr <<" *** OptControlOps: The equation solver in use can not solve"
            <<" the transposed system,\n     which the adjoint equation of a"
            <<" non-symmetric state operator needs."<< std::endl;
  return false;
}


bool OptControlOps::captureFixedLoad ()
{
  const SystemVector* b = mySim.getRHSvector(0);
  if (!b)
    return false;

  myFixedLoad.resize(b->dim(),true);
  if (b->dim() > 0)
    myFixedLoad.fill(b->getRef(),b->dim());

  return true;
}


SparseMatrix* OptControlOps::initMassMatrix ()
{
  // This copy of the mass matrix is only multiplied with, never factorized,
  // hence no equation solver is associated with it
  myMass = std::make_unique<SparseMatrix>(SparseMatrix::NONE);

  return myMass.get();
}


void OptControlOps::initLoadOperator (SparseMatrix*& C, SparseMatrix*& Ct,
                                      SparseMatrix*& Mr)
{
  myCoup  = std::make_unique<SparseMatrix>(SparseMatrix::NONE);
  myCoupT = std::make_unique<SparseMatrix>(SparseMatrix::NONE);
  // Unlike the other two, this one is factorized rather than multiplied with
  myRiesz = std::make_unique<SparseMatrix>(SparseMatrix::SUPERLU);

  C  = myCoup.get();
  Ct = myCoupT.get();
  Mr = myRiesz.get();
}


double OptControlOps::getTargetNorm () const
{
  return normTarget2 > 0.0 ? std::sqrt(normTarget2) : 0.0;
}


bool OptControlOps::setControlDofs (const std::vector<char>& nodeTypes,
                                    unsigned char nComp)
{
  const SAM* sam = mySim.getSAM();
  if (!sam)
    return false;

  ctrlDof.clear();
  if (nodeTypes.empty() && nComp == 0)
    return true; // The control lives on all degrees of freedom

  ctrlDof.assign(sam->getNoDOFs(),false);
  for (int inod = 1; inod <= sam->getNoNodes(); inod++)
  {
    if (!nodeTypes.empty() &&
        std::find(nodeTypes.begin(),nodeTypes.end(),
                  sam->getNodeType(inod)) == nodeTypes.end())
      continue;

    std::pair<int,int> dofs = sam->getNodeDOFs(inod);
    if (nComp > 0)
    {
      // A node with fewer degrees of freedom than the control has components
      // is not one the control lives on, a global Lagrange multiplier say
      if (dofs.second - dofs.first + 1 < nComp)
        continue;

      dofs.second = dofs.first + nComp - 1;
    }

    for (int idof = dofs.first; idof <= dofs.second; idof++)
      ctrlDof[idof-1] = true;
  }

  return this->getNoControls() > 0;
}


size_t OptControlOps::getNoControls () const
{
  const SAM* sam = mySim.getSAM();
  if (ctrlDof.empty())
    return sam ? sam->getNoDOFs() : 0;

  return std::count(ctrlDof.begin(),ctrlDof.end(),true);
}


void OptControlOps::maskControl (Vector& X) const
{
  for (size_t i = 0; i < ctrlDof.size() && i < X.size(); i++)
    if (!ctrlDof[i])
      X[i] = 0.0;
}


bool OptControlOps::applyMass (const Vector& X, Vector& MX) const
{
  if (!myMass || X.size() != myMass->dim(2))
    return false;

  StdVector x(X.ptr(),X.size());
  StdVector y(X.size());
  if (!myMass->multiply(x,y))
    return false;

  MX = y;
  return true;
}


bool OptControlOps::applyLoad (const Vector& F, Vector& LF) const
{
  if (!this->applyMass(F,LF))
    return false;

  if (massScale != 1.0)
    LF *= massScale;

  if (myCoup)
  {
    StdVector x(F.ptr(),F.size());
    StdVector y(F.size());
    if (!myCoup->multiply(x,y))
      return false;

    LF.add(y,coupScale);
  }

  return true;
}


bool OptControlOps::applyRiesz (Vector& X) const
{
  if (!myRiesz)
    return false;

  if (!rieszReady)
  {
    // The control space mass matrix is singular on the degrees of freedom the
    // control does not live on, so those are given a unit diagonal, which
    // makes the Riesz map leave them untouched
    for (size_t i = 1; i <= myRiesz->dim(1); i++)
      if (i > ctrlDof.size() || !ctrlDof[i-1])
        (*myRiesz)(i,i) = 1.0;
    rieszReady = true;
  }

  StdVector b(X.ptr(),X.size());
  if (!myRiesz->solve(b))
    return false;

  X = b;
  return true;
}


bool OptControlOps::applyLoadTranspose (const Vector& Y, Vector& G) const
{
  // The mass matrix part of the load operator cancels against the Riesz map
  G = Y;
  this->maskControl(G);
  if (massScale != 1.0)
    G *= massScale;

  if (myCoupT)
  {
    StdVector x(Y.ptr(),Y.size());
    StdVector y(Y.size());
    if (!myCoupT->multiply(x,y))
      return false;

    Vector Cy(y);
    if (!this->applyRiesz(Cy))
      return false;

    G.add(Cy,coupScale);
    this->maskControl(G);
  }

  return true;
}


bool OptControlOps::solveForward (const Vector& F, Vector& U)
{
  const SAM* sam = mySim.getSAM();
  SystemVector* b = mySim.getRHSvector(0);
  if (!sam || !b)
    return false;

  Vector LF;
  if (!this->applyLoad(F,LF))
    return false;

  // b = b0 + R*L*f, where b0 is the constant part of the load vector
  b->init();
  Real* bptr = b->getPtr();
  for (size_t i = 0; i < myFixedLoad.size(); i++)
    bptr[i] = myFixedLoad[i];
  sam->addToRHS(*b,LF);

  return this->solveSystem(U,1.0);
}


bool OptControlOps::solveHomogeneous (const Vector& nodalRHS, Vector& sol,
                                      bool adjoint)
{
  const SAM* sam = mySim.getSAM();
  SystemVector* b = mySim.getRHSvector(0);
  if (!sam || !b)
    return false;

  b->init();
  sam->addToRHS(*b,nodalRHS);

  return this->solveSystem(sol,0.0,adjoint && adjTranspose);
}


bool OptControlOps::solveSystem (Vector& sol, double scaleSD, bool transposed)
{
  const SAM* sam = mySim.getSAM();
  SystemMatrix* A = mySim.getLHSmatrix();
  SystemVector* b = mySim.getRHSvector(0);
  if (!sam || !A || !b)
    return false;

  utl::profiler->start("Equation solving");
  bool ok = transposed ? A->solveTranspose(*b) : A->solve(*b);
  utl::profiler->stop("Equation solving");

  return ok && sam->expandSolution(*b,sol,scaleSD);
}


bool OptControlOps::solveAdjoint (const Vector& U, Vector& P)
{
  Vector rhs;
  if (!this->applyMass(U,rhs))
    return false;

  rhs.add(myTarget,-1.0);

  return this->solveHomogeneous(rhs,P,true);
}


bool OptControlOps::applyHessian (const Vector& D, Vector& HD)
{
  // The reduced Hessian is applied by linearizing the control-to-state map,
  // which for a linear-quadratic problem means repeating the forward and
  // adjoint solves with homogeneous data. Note that only the second of the
  // two is an adjoint solve
  Vector LD, dU, MdU, dP, dG;
  if (!this->applyLoad(D,LD) || !this->solveHomogeneous(LD,dU,false))
    return false;
  if (!this->applyMass(dU,MdU) || !this->solveHomogeneous(MdU,dP,true))
    return false;
  if (!this->applyLoadTranspose(dP,dG))
    return false;

  HD = D;
  HD *= alpha;
  HD.add(dG);
  this->maskControl(HD);

  return true;
}


double OptControlOps::computeObjective (const Vector& F, const Vector& U) const
{
  Vector MU, MF;
  if (!this->applyMass(U,MU) || !this->applyMass(F,MF))
    return std::numeric_limits<double>::max();

  double misfit = U*MU - 2.0*(U*myTarget) + normTarget2;

  return 0.5*(misfit < 0.0 ? 0.0 : misfit) + 0.5*alpha*(F*MF);
}


bool OptControlOps::computeGradient (const Vector& F, const Vector& P,
                                     Vector& G) const
{
  if (F.size() != P.size())
    return false;

  Vector LtP;
  if (!this->applyLoadTranspose(P,LtP))
    return false;

  G = F;
  G *= alpha;
  G.add(LtP);
  this->maskControl(G);

  return true;
}
