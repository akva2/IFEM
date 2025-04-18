// $Id$
//==============================================================================
//!
//! \file ASMs3Drecovery.C
//!
//! \date Mar 06 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Recovery of secondary solutions for structured 3D spline FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

#include "ASMs3D.h"
#include "ItgPoint.h"
#include "Field.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "GlbL2projector.h"
#include "SparseMatrix.h"
#include "SplineUtils.h"
#include "Profiler.h"
#include <array>


bool ASMs3D::getGrevilleParameters (RealArray& prm, int dir, int basisNum) const
{
  if (dir < 0 || dir > 2) return false;

  const Go::BsplineBasis& basis = this->getBasis(basisNum)->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs3D::getQuasiInterplParameters (RealArray& prm, int dir) const
{
  if (!svol || dir < 0 || dir > 2) return false;

  const Go::BsplineBasis& basis = svol->basis(dir);

  RealArray knots_simple;
  basis.knotsSimple(knots_simple);

  prm.clear();
  for (size_t i = 0; i+1 < knots_simple.size(); i++)
  {
    prm.push_back(knots_simple[i]);
    prm.push_back(0.5*(knots_simple[i]+knots_simple[i+1]));
  }
  prm.push_back(knots_simple.back());

  return true;
}


bool ASMs3D::evaluate (const ASMbase* basis, const Vector& locVec,
                       RealArray& vec, int basisNum) const
{
  const ASMs3D* pch = dynamic_cast<const ASMs3D*>(basis);
  if (!pch) return false;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Matrix sValues;
  if (!pch->evalSolution(sValues,locVec,gpar.data()))
    return false;

  const Go::SplineVolume* svol = this->getBasis(basisNum);

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  Go::SplineVolume* vol_new =
    Go::VolumeInterpolator::regularInterpolation(svol->basis(0),
                                                 svol->basis(1),
                                                 svol->basis(2),
                                                 gpar[0], gpar[1], gpar[2],
                                                 sValues, sValues.rows(),
                                                 svol->rational(), weights);

  vec.assign(vol_new->coefs_begin(),vol_new->coefs_end());
  delete vol_new;

  return true;
}


Go::SplineVolume* ASMs3D::projectSolution (const IntegrandBase& integrand) const
{
  PROFILE2("ASMs3D::projectSolution");

  const int basis = 1;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basis))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()) || sValues.rows() == 0)
    return nullptr;

  const Go::SplineVolume* pvol = this->getBasis(basis);

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (pvol->rational())
    pvol->getWeights(weights);

  return Go::VolumeInterpolator::regularInterpolation(pvol->basis(0),
                                                      pvol->basis(1),
                                                      pvol->basis(2),
                                                      gpar[0], gpar[1], gpar[2],
                                                      sValues, sValues.rows(),
                                                      pvol->rational(),
                                                      weights);
}


Go::GeomObject* ASMs3D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMs3D::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                 const L2Integrand& integrand,
                                 bool continuous) const
{
  const size_t nnod = this->getNoProjectionNodes();

  const Go::SplineVolume* geo = this->getBasis(ASM::GEOMETRY_BASIS);
  const Go::SplineVolume* proj = this->getBasis(ASM::PROJECTION_BASIS);
  const bool separateProjBasis = proj != geo;
  const bool singleBasis = !separateProjBasis && this->getNoBasis() == 1;

  const int p1 = proj->order(0);
  const int p2 = proj->order(1);
  const int p3 = proj->order(2);
  const int n1 = proj->numCoefs(0);
  const int n2 = proj->numCoefs(1);
  const int nel1 = proj->numCoefs(0) - p1 + 1;
  const int nel2 = proj->numCoefs(1) - p2 + 1;
  const int nel3 = proj->numCoefs(2) - p3 + 1;

  int pmax = p1 > p2 ? p1 : p2;
  if (pmax < p3) pmax = p3;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? this->getNoGaussPt(pmax,true) : p1 - 1;
  const int ng2 = continuous ? ng1 : p2 - 1;
  const int ng3 = continuous ? ng2 : p3 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* zg = GaussQuadrature::getCoord(ng3);
  const double* wg = continuous ? GaussQuadrature::getWeight(ng1) : 0;
  if (!xg || !yg || !zg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  std::array<RealArray,3> gpar;
  SplineUtils::getGaussParameters(gpar[0],ng1,xg,proj->basis(0));
  SplineUtils::getGaussParameters(gpar[1],ng2,yg,proj->basis(1));
  SplineUtils::getGaussParameters(gpar[2],ng3,zg,proj->basis(2));

  // Evaluate basis functions at all integration points
  std::vector<Go::BasisPts>    spl1;
  std::vector<Go::BasisDerivs> spl2;
  if (continuous)
    geo->computeBasisGrid(gpar[0],gpar[1],gpar[2],spl2);

  if (!continuous || separateProjBasis)
    proj->computeBasisGrid(gpar[0],gpar[1],gpar[2],spl1);

  // Evaluate the secondary solution at all integration points
  Matrix sField;
  if (!integrand.evaluate(sField,gpar.data()))
  {
    std::cerr <<" *** ASMs3D::assembleL2matrices: Failed for patch "<< idx+1
              <<" nPoints="<< gpar[0].size()*gpar[1].size()*gpar[2].size()
              << std::endl;
    return false;
  }

  double dV = 1.0;
  Vector phi(p1*p2*p3);
  Matrix dNdu, Xnod, J;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i3 = 0; i3 < nel3; i3++)
    for (int i2 = 0; i2 < nel2; i2++)
      for (int i1 = 0; i1 < nel1; i1++, iel++)
      {
        int ip = ((i3*ng2*nel2 + i2)*ng1*nel1 + i1)*ng3;
        IntVec lmnpc;
        if (!singleBasis && proj->knotSpan(0,i1+p1-1) > 0.0 &&
                            proj->knotSpan(1,i2+p2-1) > 0.0 &&
                            proj->knotSpan(2,i3+p3-1) > 0.0)
        {
          // Establish nodal point correspondance for the projection element
          int vidx, widx;
          lmnpc.reserve(phi.size());
          if (separateProjBasis)
            widx = ((spl1[ip].left_idx[2]-p3+1)*n1*n2 +
                    (spl1[ip].left_idx[1]-p2+1)*n1    +
                    (spl1[ip].left_idx[0]-p1+1));
          else
            widx = ((spl2[ip].left_idx[2]-p3+1)*n1*n2 +
                    (spl2[ip].left_idx[1]-p2+1)*n1    +
                    (spl2[ip].left_idx[0]-p1+1));

          for (int k = 0; k < p3; k++, widx += n1*n2)
            for (int j = vidx = 0; j < p2; j++, vidx += n1)
              for (int i = 0; i < p1; i++)
                lmnpc.push_back(widx+vidx+i);
        }
        const IntVec& mnpc = singleBasis ? MNPC[iel] : lmnpc;
        if (mnpc.empty())
          continue;

        if (continuous)
        {
          // Set up control point (nodal) coordinates for current element
          if (singleBasis) {
            if (!this->getElementCoordinates(Xnod,1+iel))
              return false;
            else if ((dV = 0.125*this->getParametricVolume(1+iel)) < 0.0)
              return false; // topology error (probably logic error)
          } else {
            if (!this->getElementCoordinatesPrm(Xnod,gpar[0][i1*ng1],
                                                gpar[1][i2*ng2],gpar[2][i3*ng3]))
              return false;
            else if ((dV = 0.125 * proj->knotSpan(0, mnpc.back() % n1)
                                 * proj->knotSpan(1,(mnpc.back() / n1) % n2)
                                 * proj->knotSpan(2, mnpc.back() / (n1*n2))) < 0.0)
              return false;
          }
        }

        // --- Integration loop over all Gauss points in each direction --------

        Matrix eA(p1*p2*p3, p1*p2*p3);
        Vectors eB(sField.rows(), Vector(p1*p2*p3));
        for (int k = 0; k < ng3; k++, ip += ng2*(nel2-1)*ng1*nel1)
          for (int j = 0; j < ng2; j++, ip += ng1*(nel1-1))
            for (int i = 0; i < ng1; i++, ip++)
            {
              if (continuous)
                SplineUtils::extractBasis(spl2[ip],phi,dNdu);

              if (!continuous || separateProjBasis)
                phi = spl1[ip].basisValues;

              // Compute the Jacobian inverse and derivatives
              double dJw = dV;
              if (continuous)
              {
                dJw *= wg[i]*wg[j]*wg[k]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
                if (dJw == 0.0) continue; // skip singular points
              }

              // Integrate the mass matrix
              eA.outer_product(phi, phi, true, dJw);

              // Integrate the rhs vector B
              for (size_t r = 1; r <= sField.rows(); r++)
                eB[r-1].add(phi,sField(r,ip+1)*dJw);
            }

        for (int i = 0; i < p1*p2*p3; ++i) {
          for (int j = 0; j < p1*p2*p3; ++j)
            A(mnpc[i]+1, mnpc[j]+1) += eA(i+1, j+1);

          int jp = mnpc[i]+1;
          for (size_t r = 0; r < sField.rows(); r++, jp += nnod)
            B(jp) += eB[r](1+i);
        }
      }

  return true;
}


#include "ASMs3DInterpolate.C" // TODO: inline these methods instead...


/*!
  \note A Variation Diminishing Spline Approximation is used here as the
  regular interpolation method in GoTools only works with uniform knots.
*/

bool ASMs3D::evaluate (const Field* field, RealArray& vec, int basisNum) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Vector sValues;
  sValues.reserve(gpar[0].size()*gpar[1].size()*gpar[2].size());
  for (double w : gpar[2])
    for (double v : gpar[1])
      for (double u : gpar[0])
        sValues.push_back(field->valueFE(ItgPoint(u,v,w)));

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  Go::SplineVolume* vol_new =
    VariationDiminishingSplineApproximation(this->getBasis(basisNum),sValues,1);

  vec.assign(vol_new->coefs_begin(),vol_new->coefs_end());
  delete vol_new;

  return true;
}


Go::SplineVolume* ASMs3D::projectSolutionLeastSquare (const IntegrandBase& integrand) const
{
  if (!svol) return nullptr;

  // Get Gaussian quadrature points and weights
  int p = svol->order(0);
  for (int d = 1; d < 3; d++)
    if (p < svol->order(d)) p = svol->order(d);
  const int ng = this->getNoGaussPt(p,true);
  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = GaussQuadrature::getWeight(ng);
  if (!xg || !wg) return nullptr;

  std::array<RealArray,3> gpar, wgpar;
  for (int dir = 0; dir < 3; dir++)
  {
    // Get Gauss parameter values and associated weights
    const Go::BsplineBasis& basis = svol->basis(dir);
    SplineUtils::getGaussParameters(gpar[dir],ng,xg,basis);
    RealArray::const_iterator knotit = basis.begin();
    RealArray& tmp = wgpar[dir];
    tmp.reserve(ng*(basis.numCoefs()-basis.order()));
    for (int i = basis.order(); i <= basis.numCoefs(); i++)
    {
      double d = knotit[i] - knotit[i-1];
      for (int j = 0; j < ng; j++)
        tmp.push_back(d > 0.0 ? wg[j]*d*0.5 : 0.0);
    }
  }

  // Evaluate the secondary solution at all sampling points (Gauss points)
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  return leastsquare_approximation(svol->basis(0),
                                   svol->basis(1),
                                   svol->basis(2),
                                   gpar[0], gpar[1], gpar[2],
                                   wgpar[0], wgpar[1], wgpar[2],
                                   sValues,
                                   sValues.rows(),
                                   svol->rational(),
                                   weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLocal (const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getQuasiInterplParameters(gpar[dir],dir))
      return nullptr;

  for (int dir = 2; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  return quasiInterpolation(svol->basis(0),
                            svol->basis(1),
                            svol->basis(2),
                            gpar[0], gpar[1], gpar[2],
                            sValues,
                            sValues.rows(),
                            svol->rational(),
                            weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLocalApprox(const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  // Project onto the geometry basis
  return VariationDiminishingSplineApproximation(svol.get(),sValues,sValues.rows());
}
