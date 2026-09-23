//==============================================================================
//!
//! \file GlbNodalSystem.C
//!
//! \date Sep 19 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Global integral for assembly in the unconstrained nodal ordering.
//!
//==============================================================================

#include "GlbNodalSystem.h"

#include "ASMbase.h" // for IntMat
#include "ElmMats.h"
#include "SAM.h"
#include "SparseMatrix.h"

#include <algorithm>
#include <iostream>


GlbNodalSystem::GlbNodalSystem (const SAM& sam,
                                const std::vector<MatrixTarget>& mats,
                                Vector* vec, RealArray* scl)
  : mySam(sam), myMats(mats), myVec(vec), myScl(scl)
{
  this->preAssemble();
}


GlbNodalSystem::GlbNodalSystem (const SAM& sam, SparseMatrix* mat,
                                Vector* vec, RealArray* scl)
  : mySam(sam), myVec(vec), myScl(scl)
{
  if (mat)
    myMats.push_back({mat,0,false});

  this->preAssemble();
}


void GlbNodalSystem::preAssemble ()
{
  if (myMats.empty())
    return;

  // Bake the sparsity patterns, such that the element loop below only updates
  // existing entries and therefore can run in parallel
  IntMat mmnpc(mySam.getNoElms());
  for (size_t iel = 0; iel < mmnpc.size(); iel++)
    if (this->elementDofs(mmnpc[iel],1+iel))
      for (int& idof : mmnpc[iel])
        --idof; // SparseMatrix::preAssemble expects zero-based indices

  for (const MatrixTarget& target : myMats)
    if (target.mat)
    {
      target.mat->resize(mySam.getNoDOFs(),mySam.getNoDOFs(),true);
      target.mat->preAssemble(mmnpc,mmnpc.size());
    }
}


bool GlbNodalSystem::elementDofs (IntVec& meen, int elmId) const
{
  IntVec mnpc;
  if (!mySam.getElmNodes(mnpc,elmId))
    return false;

  meen.clear();
  meen.reserve(mnpc.size());
  for (int inod : mnpc)
  {
    // A negative node number means the DOFs of that node are not present in
    // the element matrices, but they still occupy their place in the ordering
    std::pair<int,int> dofs = mySam.getNodeDOFs(inod > 0 ? inod : -inod);
    for (int idof = dofs.first; idof <= dofs.second; idof++)
      meen.push_back(inod > 0 ? idof : 0);
  }

  return true;
}


void GlbNodalSystem::initialize (char)
{
  for (const MatrixTarget& target : myMats)
    if (target.mat)
      target.mat->init(); // Zero the values, but retain the sparsity pattern

  if (myVec)
    myVec->resize(mySam.getNoDOFs(),true);

  if (myScl)
    std::fill(myScl->begin(),myScl->end(),0.0);
}


bool GlbNodalSystem::sizeError (const char* what, int elmId, size_t nedof)
{
  std::cerr <<" *** GlbNodalSystem::assemble: Element "<< what <<" of element "
            << elmId <<" is inconsistent with its "<< nedof
            <<" degrees of freedom."<< std::endl;
  return false;
}


bool GlbNodalSystem::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmMats* elMat = dynamic_cast<const ElmMats*>(elmObj);
  if (!elMat)
    return false;

  IntVec meen;
  if (!this->elementDofs(meen,elmId))
    return false;

  const size_t nedof = meen.size();

  // An element quantity may be smaller than the element it belongs to, as
  // global Lagrange multipliers are appended at the end of the element
  // ordering and the integrands here do not contribute to them
  for (const MatrixTarget& target : myMats)
  {
    if (!target.mat || target.idx >= elMat->A.size())
      continue;

    const Matrix& eM = elMat->A[target.idx];
    if (eM.rows() > nedof || eM.cols() != eM.rows())
      return sizeError("matrix",elmId,nedof);

    for (size_t j = 1; j <= eM.cols(); j++)
      if (meen[j-1] > 0)
        for (size_t i = 1; i <= eM.rows(); i++)
          if (meen[i-1] > 0)
          {
            if (target.transposed)
              (*target.mat)(meen[j-1],meen[i-1]) += eM(i,j);
            else
              (*target.mat)(meen[i-1],meen[j-1]) += eM(i,j);
          }
  }

  if (myVec && !elMat->b.empty())
  {
    const Vector& eS = elMat->b.front();
    if (eS.size() > nedof)
      return sizeError("vector",elmId,nedof);

    for (size_t i = 1; i <= eS.size(); i++)
      if (meen[i-1] > 0)
        (*myVec)(meen[i-1]) += eS(i);
  }

  // All elements contribute to each of the global scalars, so these are not
  // covered by the node-disjointness of the thread groups
  if (myScl)
    for (size_t i = 0; i < elMat->c.size() && i < myScl->size(); i++)
    {
      double& sum = (*myScl)[i];
#pragma omp atomic
      sum += elMat->c[i];
    }

  return true;
}
