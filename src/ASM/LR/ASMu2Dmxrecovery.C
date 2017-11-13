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
#include "matrix.h"
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

  size_t nnod = this->getNoProjectionNodes();
  double dA = 0.0;
  Vectors phi(2);
  Matrices dNdu(2);
  Matrix sField, Xnod, Jac;
  std::vector<Go::BasisDerivsSf> spl1(2);
  std::vector<Go::BasisPtsSf> spl0(2);


  // === Assembly loop over all elements in the patch ==========================

  for (auto el1 : m_basis[geoBasis-1]->getAllElements())
  {
    double uh = (el1->umin()+el1->umax())/2.0;
    double vh = (el1->vmin()+el1->vmax())/2.0;
    std::vector<size_t> els;
    els.push_back(projBasis->getElementContaining(uh, vh)+1);
    els.push_back(m_basis[geoBasis-1]->getElementContaining(uh,vh)+1);

    if (continuous)
    {
      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,els[1]))
        return false;
      else if ((dA = 0.25*this->getParametricArea(els[1])) < 0.0)
        return false; // topology error (probably logic error)
    }

    // Compute parameter values of the Gauss points over this element
    RealArray gpar[2], unstrGpar[2];
    this->getGaussPointParameters(gpar[0],0,ng1,els[1],xg);
    this->getGaussPointParameters(gpar[1],1,ng2,els[1],yg);

    // convert to unstructred mesh representation
    expandTensorGrid(gpar, unstrGpar);

    // Evaluate the secondary solution at all integration points
    if (!this->evalSolution(sField,integrand,unstrGpar))
      return false;

    // set up basis function size (for extractBasis subroutine)
    const LR::Element* elm = projBasis->getElement(els[0]-1);
    phi[0].resize(elm->nBasisFunctions());
    phi[1].resize(el1->nBasisFunctions());
    IntVec lmnpc;
    if (projBasis != m_basis[0]) {
      lmnpc.reserve(phi[0].size());
      for (const LR::Basisfunction* f : elm->support())
        lmnpc.push_back(f->getId());
    }
    const IntVec& mnpc = projBasis == m_basis[0] ? MNPC[els[1]-1] : lmnpc;

    // --- Integration loop over all Gauss points in each direction ----------
    Matrix eA(phi[0].size(), phi[0].size());
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

        // Integrate the mass matrix
        eA.outer_product(phi[0], phi[0], true, dJw);

        // Integrate the rhs vector B
        size_t ncmp = sField.rows();
        for (size_t ii = 0; ii < phi[0].size(); ii++)
        {
          int inod = mnpc[ii]+1;
          for (size_t r = 1; r <= ncmp; ++r)
            B(inod + (r-1)*nnod) += phi[0][ii]*sField(r,ip+1)*dJw;
        }
      }

    for (size_t j = 1; j <= eA.cols(); ++j)
      for (size_t i = 1; i <= eA.rows(); ++i)
        A(mnpc[i-1]+1, mnpc[j-1]+1) += eA(i,j);
  }

  return true;
}
