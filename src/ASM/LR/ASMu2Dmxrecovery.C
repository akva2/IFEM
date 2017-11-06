// $Id$
//==============================================================================
//!
//! \file ASMu2Dmxrecovery.C
//!
//! \date May 11 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Recovery techniques for unstructured mixed LR B-splines.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"

#include "ASMu2Dmx.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "SplineUtils.h"
#include "Profiler.h"
#include <numeric>


/*!
  \brief Expands a tensor parametrization point to an unstructured one.
*/

static void expandTensorGrid (const RealArray* in, RealArray* out)
{
  out[0].resize(in[0].size()*in[1].size());
  out[1].resize(in[0].size()*in[1].size());

  size_t i, j, ip = 0;
  for (j = 0; j < in[1].size(); j++)
    for (i = 0; i < in[0].size(); i++, ip++) {
      out[0][ip] = in[0][i];
      out[1][ip] = in[1][j];
    }
}


bool ASMu2Dmx::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                   const IntegrandBase& integrand,
                                   bool continuous) const
{
  const int p1 = projBasis->order(0);
  const int p2 = projBasis->order(1);

  // Get Gaussian quadrature points
  const int ng1 = continuous ? nGauss : p1 - 1;
  const int ng2 = continuous ? nGauss : p2 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* wg = continuous ? GaussQuadrature::getWeight(nGauss) : nullptr;
  if (!xg || !yg) return false;
  if (continuous && !wg) return false;

  double dA = 0.0;
  Vectors phi(2);
  Matrices dNdu(2);
  Matrix sField, Xnod, Jac;
  std::vector<Go::BasisDerivsSf> spl1(2);
  std::vector<Go::BasisPtsSf> spl0(2);


  // === Assembly loop over all elements in the patch ==========================

  std::vector<LR::Element*>::iterator el1 = m_basis[geoBasis-1]->elementBegin();
  for (int iel = 1; el1 != m_basis[geoBasis-1]->elementEnd(); ++el1, ++iel)
  {
    double uh = ((*el1)->umin()+(*el1)->umax())/2.0;
    double vh = ((*el1)->vmin()+(*el1)->vmax())/2.0;
    std::vector<size_t> els;
    els.push_back(projBasis->getElementContaining(uh, vh)+1);
    els.push_back(m_basis[geoBasis-1]->getElementContaining(uh,vh)+1);

    int geoEl = els[1];

    if (continuous)
    {
      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,geoEl))
        return false;
      else if ((dA = 0.25*this->getParametricArea(geoEl)) < 0.0)
        return false; // topology error (probably logic error)
    }

    // Compute parameter values of the Gauss points over this element
    RealArray gpar[2], unstrGpar[2];
    this->getGaussPointParameters(gpar[0],0,ng1,geoEl,xg);
    this->getGaussPointParameters(gpar[1],1,ng2,geoEl,yg);

    // convert to unstructred mesh representation
    expandTensorGrid(gpar, unstrGpar);

    // Evaluate the secondary solution at all integration points
    if (!this->evalSolution(sField,integrand,unstrGpar))
      return false;

    // set up basis function size (for extractBasis subroutine)
    const LR::Element* elm = projBasis->getElement(els[0]-1);
    phi[0].resize(elm->nBasisFunctions());
    phi[1].resize(m_basis[geoBasis-1]->getElement(els[1]-1)->nBasisFunctions());
    IntVec lmnpc;
    if (projBasis != m_basis[0]) {
      lmnpc.reserve(phi[0].size());
      for (const LR::Basisfunction* f : elm->support())
        lmnpc.push_back(f->getId());
    }
    const IntVec& mnpc = projBasis == m_basis[0] ? MNPC[els[0]-1] : lmnpc;

    // --- Integration loop over all Gauss points in each direction ----------
    int ip = 0;
    for (int j = 0; j < ng2; j++)
      for (int i = 0; i < ng1; i++, ip++)
      {
        if (continuous)
        {
          projBasis->computeBasis(gpar[0][i], gpar[1][j], spl1[0], els[0]-1);
          SplineUtils::extractBasis(spl1[0],phi[0],dNdu[0]);
          m_basis[geoBasis-1]->computeBasis(gpar[0][i], gpar[1][j],
                                            spl1[1], els[1]-1);
          SplineUtils::extractBasis(spl1[1], phi[1], dNdu[1]);
        }
        else
        {
          projBasis->computeBasis(gpar[0][i], gpar[1][j], spl0[0], els[0]-1);
          phi[0] = spl0[0].basisValues;
        }

        // Compute the Jacobian inverse and derivatives
        double dJw = 1.0;
        if (continuous)
        {
          dJw = dA*wg[i]*wg[j]*utl::Jacobian(Jac,dNdu[1],Xnod,dNdu[1],false);
          if (dJw == 0.0) continue; // skip singular points
        }

        // Integrate the linear system A*x=B
        size_t ncmp = sField.rows();
        for (size_t ii = 0; ii < phi[0].size(); ii++)
        {
          int inod = mnpc[ii];
          for (size_t jj = 0; jj < phi[0].size(); jj++)
          {
            int jnod = mnpc[jj];
            for (size_t k = 1; k <= ncmp; ++k)
              A(inod*ncmp+k, jnod*ncmp+k) += phi[0][ii]*phi[0][jj]*dJw;
          }
          for (size_t k = 1; k <= ncmp; ++k)
            B(inod*ncmp+k) += phi[0][ii]*sField(k,ip+1)*dJw;
        }
      }
  }

  return true;
}


bool ASMu2Dmx::globalL2projection (Matrix& sField,
                                   const IntegrandBase& integrand,
                                   bool continuous) const
{
  for (size_t b = 0; b < m_basis.size(); b++)
    if (!m_basis[b]) return true; // silently ignore empty patches

  PROFILE2("ASMu2Dmx::globalL2");

  // Assemble the projection matrices
  size_t nnod = this->getNoProjectionNodes() * integrand.getNoFields(2);
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(nnod);
  A.redim(nnod,nnod);

  if (!this->assembleL2matrices(A,B,integrand,continuous))
    return false;

#if SP_DEBUG > 1
  std::cout <<"---- Matrix A -----\n"<< A
            <<"-------------------"<< std::endl;
  std::cout <<"---- Vector B -----\n"<< B
            <<"-------------------"<< std::endl;
#endif

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the control-point values of the projected field
  sField.resize(integrand.getNoFields(2), nb[0]);

  size_t inod = 1, jnod = 1;
  for (size_t i = 1; i <= this->getNoProjectionNodes(); i++, inod++)
    for (size_t j = 1; j <= integrand.getNoFields(2); j++, jnod++)
      sField(j,inod) = B(jnod);

#if SP_DEBUG > 1
  std::cout <<"- Solution Vector -"<< sField
            <<"-------------------"<< std::endl;
#endif
  return true;
}
