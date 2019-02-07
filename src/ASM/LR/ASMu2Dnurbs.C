// $Id$
//==============================================================================
//!
//! \file ASMu2Dnurbs.C
//!
//! \date Nov 28 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D NURBS FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "ASMu2Dnurbs.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "LagrangeInterpolator.h"
#include "ElementBlock.h"
#include "MPC.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Function.h"
#include "Vec3Oper.h"
#include <array>
#include <fstream>


ASMu2Dnurbs::ASMu2Dnurbs (unsigned char n_s, unsigned char n_f)
  : ASMu2D(n_s, n_f)
{
  isNurbs = true;
}


ASMu2Dnurbs::ASMu2Dnurbs (const ASMu2Dnurbs& patch, unsigned char n_f)
  : ASMu2D(patch, n_f)
{
  isNurbs = patch.isNurbs;
}

/*!
  \brief Create a dim+1 dimensional LRSplineSurface from a tensor NURBS surface.
*/
static LR::LRSplineSurface* createFromTensor(const Go::SplineSurface* srf)
{
  return new LR::LRSplineSurface(srf->numCoefs_u(), srf->numCoefs_v(),
                                 srf->order_u(), srf->order_v(),
                                 srf->basis_u().begin(),
                                 srf->basis_v().begin(),
                                 srf->rcoefs_begin(), srf->dimension()+1);
}


bool ASMu2Dnurbs::read (std::istream& is)
{
  if (shareFE) return true;

  // read inputfile as either an LRSpline file directly
  // or a tensor product B-spline and convert
  char firstline[256];
  is.getline(firstline, 256);
  if (strncmp(firstline, "# LRSPLINE", 10) == 0) {
    lrspline.reset(new LR::LRSplineSurface());
    is >> *lrspline;
  } else { // probably a SplineSurface, so we'll read that and convert
    tensorspline = new Go::SplineSurface();
    is >> *tensorspline;
    if (tensorspline->rational())
      lrspline.reset(createFromTensor(tensorspline));
    else {
      lrspline.reset(new LR::LRSplineSurface(tensorspline));
      std::cout << "LR-nurbs requested but input is a spline" << std::endl;
      isNurbs = false;
    }
  }

  // Eat white-space characters to see if there is more data to read
  char c;
  while (is.get(c))
    if (!isspace(c)) {
      is.putback(c);
      break;
    }

  if (!is.good() && !is.eof())
  {
    std::cerr <<" *** ASMu2Dnurbs::read: Failure reading spline data"<< std::endl;
    lrspline.reset();
    return false;
  }
  else if (lrspline->dimension() < 2)
  {
    std::cerr <<" *** ASMu2Dnurbs::read: Invalid lrspline patch, dim="
              << lrspline->dimension() << std::endl;
    lrspline.reset();
    return false;
  }
  else if (lrspline->dimension() < nsd)
  {
    std::cout <<"  ** ASMu2Dnurbs::read: The dimension of this lrsplineace patch "
        << lrspline->dimension() <<" is less than nsd="<< nsd
        <<".\n                   Resetting nsd to "<< lrspline->dimension()
        <<" for this patch."<< std::endl;
    nsd = lrspline->dimension();
  }

  geo = lrspline.get();
  return true;
}


bool ASMu2Dnurbs::uniformRefine (int dir, int nInsert)
{
  if (!this->ASMu2D::uniformRefine(dir,nInsert))
    return false;

  if (!isNurbs)
    return true;

  lrspline.reset(createFromTensor(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu2Dnurbs::refine (int dir, const RealArray& xi, double scale)
{
  if (!this->ASMu2D::refine(dir,xi,scale))
    return false;

  if (!isNurbs)
    return true;

  lrspline.reset(createFromTensor(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu2Dnurbs::raiseOrder (int ru, int rv)
{
  if (!this->ASMu2D::raiseOrder(ru,rv))
    return false;

  lrspline.reset(createFromTensor(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu2Dnurbs::evaluateBasis (int iel, FiniteElement& fe,
                                 int derivs) const
{
  if (!isNurbs)
    return this->ASMu2D::evaluateBasis(iel,fe,derivs);

  const LR::Element* el = lrspline->getElement(iel);
  if (!el) return false;

  fe.xi  = 2.0*(fe.u - el->umin()) / (el->umax() - el->umin()) - 1.0;
  fe.eta = 2.0*(fe.v - el->vmin()) / (el->vmax() - el->vmin()) - 1.0;
  RealArray Nu = bezier_u.computeBasisValues(fe.xi, derivs);
  RealArray Nv = bezier_v.computeBasisValues(fe.eta,derivs);
  const Matrix& C = bezierExtract[iel];

  size_t wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;

  if (derivs < 1) {
    Matrix B;
    B.outer_product(Nu,Nv);
    Vector N  = C*static_cast<const Vector&>(B);

    double W = 0.0;
    for (size_t i = 1; i <= N.size(); ++i)
      W += N(i)*lrspline->getBasisfunction(MNPC[iel][i-1])->cp(wi);

    fe.N.resize(N.size(),true);
    for (size_t i = 1; i <= N.size(); ++i)
      fe.N(i) = N(i)*lrspline->getBasisfunction(MNPC[iel][i-1])->cp(wi)/W;

#if SP_DEBUG > 2
    if (fabs(fe.N.sum()-1.0) > 1.0e-10)
      std::cerr <<"fe.N do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(static_cast<const Vector&>(B).sum()-1.0) > 1.0e-10)
      std::cerr <<"Bezier basis do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else
      return true; // The basis is OK

    return false;
#endif
  }
  else {
    int p = lrspline->order(0)*lrspline->order(1);

    Vector B(p);
    Vector Bu(p); // Bezier basis functions differentiated wrt u
    Vector Bv(p); // Bezier basis functions differentiated wrt v

    size_t i, j, k = 0;
    for (j = 0; j < Nv.size(); j+=(derivs+1))
      for (i = 0; i < Nu.size(); i+=(derivs+1), k++) {
        B[k]  = Nu[i  ]*Nv[j  ];
        Bu[k] = Nu[i+1]*Nv[j  ];
        Bv[k] = Nu[i  ]*Nv[j+1];
      }

    Vector N = C*B;
    Vector dNxi = C*Bu;
    Vector dNeta = C*Bv;

    double W = 0.0;
    double Wderxi = 0.0;
    double Wdereta = 0.0;
    for (size_t i = 1; i <= N.size(); ++i) {
      double w = lrspline->getBasisfunction(MNPC[iel][i-1])->cp(wi);
      W += N(i)*w;
      Wderxi += dNxi(i)*w;
      Wdereta += dNeta(i)*w;
    }

    fe.N.resize(N.size(),true);
    Matrix dNdu(el->nBasisFunctions(),2);
    for (size_t i = 1; i <= N.size(); ++i) {
      double w = lrspline->getBasisfunction(MNPC[iel][i-1])->cp(wi);
      fe.N(i) = N(i)*w/W;
      dNdu(i,1) = (dNxi(i)*W - N(i)*Wderxi)*w/(W*W);
      dNdu(i,2) = (dNeta(i)*W - N(i)*Wdereta)*w/(W*W);
    }

    Matrix Xnod, Jac;
    this->getElementCoordinates(Xnod,iel+1);
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

#if SP_DEBUG > 2
    if (fabs(fe.N.sum()-1.0) > 1.0e-10)
      std::cerr <<"fe.N do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(B.sum()-1.0) > 1.0e-10)
      std::cerr <<"Bezier basis do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(dNdu.getColumn(1).sum()) > 1.0e-10)
      std::cerr <<"dNdu not sums to zero at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(dNdu.getColumn(1).sum()) > 1.0e-10)
      std::cerr <<"dNdv not sums to zero at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(Bu.sum()) > 1.0e-10 || fabs(Bv.sum()) > 1.0e-10)
      std::cerr <<"Bezier derivatives do not sum to zero at integration point #"
                << fe.iGP << std::endl;
    else
      return true; // The basis is OK

    return false;
#endif
  }
  return true;
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf& bas, int iel) const
{
  if (!isNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  Go::BasisDerivsSf tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  size_t wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;

  double W = 0.0;
  double Wderxi = 0.0;
  double Wdereta = 0.0;
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);
    W += tmp.basisValues[i]*w;
    Wderxi += tmp.basisDerivs_u[i]*w;
    Wdereta += tmp.basisDerivs_v[i]*w;
  }

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);
    bas.basisValues[i] = tmp.basisValues[i]*w/W;
    bas.basisDerivs_u[i] = (tmp.basisDerivs_u[i]*W-tmp.basisValues[i]*Wderxi)*w/(W*W);
    bas.basisDerivs_v[i] = (tmp.basisDerivs_v[i]*W-tmp.basisValues[i]*Wdereta)*w/(W*W);
  }
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf2& bas, int iel) const
{
  if (!isNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  Go::BasisDerivsSf2 tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  size_t wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;

  double W = 0.0;
  double Wderxi = 0.0;
  double Wdereta = 0.0;
  double Wderxixi = 0.0;
  double Wderetaeta = 0.0;
  double Wderxieta = 0.0;
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);
    W += tmp.basisValues[i]*w;
    Wderxi += tmp.basisDerivs_u[i]*w;
    Wdereta += tmp.basisDerivs_v[i]*w;
    Wderxixi += tmp.basisDerivs_uu[i]*w;
    Wderetaeta += tmp.basisDerivs_vv[i]*w;
    Wderxieta += tmp.basisDerivs_uv[i]*w;
  }

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);
    bas.basisValues[i] = tmp.basisValues[i]*w/W;
    double H1 = (tmp.basisDerivs_u[i]*W-tmp.basisValues[i]*Wderxi);
    bas.basisDerivs_u[i] = H1*w/(W*W);
    double H2 = (tmp.basisDerivs_v[i]*W-tmp.basisValues[i]*Wdereta);
    bas.basisDerivs_v[i] = H2*w/(W*W);
    double dH1dx = tmp.basisDerivs_uu[i]*W-tmp.basisValues[i]*Wderxixi;
    bas.basisDerivs_uu[i] = (dH1dx*W-2*H1*Wderxi)*w/(W*W*W);
    double dH2dy = tmp.basisDerivs_vv[i]*W-tmp.basisValues[i]*Wderetaeta;
    bas.basisDerivs_vv[i] = (dH2dy*W-2*H2*Wdereta)*w/(W*W*W);
    double dH1dy = tmp.basisDerivs_uv[i]*W + tmp.basisDerivs_u[i]*Wdereta -
                   tmp.basisDerivs_v[i]*Wderxi - tmp.basisValues[i]*Wderxieta;
    bas.basisDerivs_uv[i] = (dH1dy*W-2*H1*Wdereta)*w/(W*W*W);
  }
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf3& bas, int iel) const
{
  if (!isNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  std::cerr << "Third order derivatives not implemented for NURBS." << std::endl;
  assert(0);
}
