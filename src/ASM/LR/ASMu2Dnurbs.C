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
  noNurbs = false;
}


ASMu2Dnurbs::ASMu2Dnurbs (const ASMu2Dnurbs& patch, unsigned char n_f)
  : ASMu2D(patch, n_f)
{
  noNurbs = patch.noNurbs;
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
    std::cout << "LR-nurbs requested but input is a spline" << std::endl;
    noNurbs = true;
  } else { // probably a SplineSurface, so we'll read that and convert
    tensorspline = new Go::SplineSurface();
    is >> *tensorspline;
    if (tensorspline->rational())
      lrspline.reset(createFromTensor(tensorspline));
    else {
      lrspline.reset(new LR::LRSplineSurface(tensorspline));
      std::cout << "LR-nurbs requested but input is a spline" << std::endl;
      noNurbs = true;
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
    std::cout <<"  ** ASMu2Dnurbs::read: The dimension of this lrspline patch "
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

  if (noNurbs)
    return true;

  lrspline.reset(createFromTensor(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu2Dnurbs::refine (int dir, const RealArray& xi, double scale)
{
  if (!this->ASMu2D::refine(dir,xi,scale))
    return false;

  if (noNurbs)
    return true;

  lrspline.reset(createFromTensor(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu2Dnurbs::raiseOrder (int ru, int rv)
{
  if (!this->ASMu2D::raiseOrder(ru,rv))
    return false;

  if (noNurbs)
    return true;

  lrspline.reset(createFromTensor(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu2Dnurbs::evaluateBasis (int iel, FiniteElement& fe,
                                 int derivs) const
{
  if (noNurbs)
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
                               Go::BasisPtsSf& bas, int iel,
                               const LR::LRSplineSurface* spline) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel,spline);

  if (!spline)
    spline = lrspline.get();

  Go::BasisPtsSf tmp;
  spline->computeBasis(u,v,tmp,iel);
  size_t wi = spline->getBasisfunction(MNPC[iel][0])->dim()-1;

  double W = 0.0;
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = spline->getBasisfunction(MNPC[iel][i])->cp(wi);
    W += tmp.basisValues[i]*w;
  }

  bas.preparePts(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = spline->getBasisfunction(MNPC[iel][i])->cp(wi);
    bas.basisValues[i] = tmp.basisValues[i]*w/W;
  }
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf& bas, int iel,
                               const LR::LRSplineSurface* spline) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel,spline);

  if (!spline)
    spline = lrspline.get();

  Go::BasisDerivsSf tmp;
  spline->computeBasis(u,v,tmp,iel);
  size_t wi = spline->getBasisfunction(MNPC[iel][0])->dim()-1;

  double W = 0.0;
  double Wderxi = 0.0;
  double Wdereta = 0.0;
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = spline->getBasisfunction(MNPC[iel][i])->cp(wi);
    W += tmp.basisValues[i]*w;
    Wderxi += tmp.basisDerivs_u[i]*w;
    Wdereta += tmp.basisDerivs_v[i]*w;
  }

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = spline->getBasisfunction(MNPC[iel][i])->cp(wi);
    bas.basisValues[i] = tmp.basisValues[i]*w/W;
    bas.basisDerivs_u[i] = (tmp.basisDerivs_u[i]*W-tmp.basisValues[i]*Wderxi)*w/(W*W);
    bas.basisDerivs_v[i] = (tmp.basisDerivs_v[i]*W-tmp.basisValues[i]*Wdereta)*w/(W*W);
  }
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf2& bas, int iel) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  Go::BasisDerivsSf2 tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  size_t wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;

  double W = 0.0;
  double Wderxi = 0.0;
  double Wdereta = 0.0;
  double Wder2xi = 0.0;
  double Wder2eta = 0.0;
  double Wderxieta = 0.0;
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);
    W += tmp.basisValues[i]*w;
    Wderxi += tmp.basisDerivs_u[i]*w;
    Wdereta += tmp.basisDerivs_v[i]*w;
    Wder2xi += tmp.basisDerivs_uu[i]*w;
    Wder2eta += tmp.basisDerivs_vv[i]*w;
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
    double dH1dx = tmp.basisDerivs_uu[i]*W-tmp.basisValues[i]*Wder2xi;
    bas.basisDerivs_uu[i] = (dH1dx*W-2*H1*Wderxi)*w/(W*W*W);
    double dH2dy = tmp.basisDerivs_vv[i]*W-tmp.basisValues[i]*Wder2eta;
    bas.basisDerivs_vv[i] = (dH2dy*W-2*H2*Wdereta)*w/(W*W*W);
    double dH1dy = tmp.basisDerivs_uv[i]*W + tmp.basisDerivs_u[i]*Wdereta -
                   tmp.basisDerivs_v[i]*Wderxi - tmp.basisValues[i]*Wderxieta;
    bas.basisDerivs_uv[i] = (dH1dy*W-2*H1*Wdereta)*w/(W*W*W);
  }
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf3& bas, int iel) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  Go::BasisDerivsSf3 tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  size_t wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;

  double W = 0.0;
  double Wderxi = 0.0;
  double Wdereta = 0.0;
  double Wder2xi = 0.0;
  double Wder2eta = 0.0;
  double Wderxieta = 0.0;
  double Wder3xi = 0.0;
  double Wder3eta = 0.0;
  double Wder2xieta = 0.0;
  double Wderxi2eta = 0.0;
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);
    W += tmp.basisValues[i]*w;
    Wderxi += tmp.basisDerivs_u[i]*w;
    Wdereta += tmp.basisDerivs_v[i]*w;
    Wder2xi += tmp.basisDerivs_uu[i]*w;
    Wder2eta += tmp.basisDerivs_vv[i]*w;
    Wderxieta += tmp.basisDerivs_uv[i]*w;
    Wder3xi += tmp.basisDerivs_uuu[i]*w;
    Wder3eta += tmp.basisDerivs_vvv[i]*w;
    Wder2xieta = tmp.basisDerivs_uuv[i]*w;
    Wderxi2eta = tmp.basisDerivs_uvv[i]*w;
  }

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); ++i) {
    double w = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);
    bas.basisValues[i] = tmp.basisValues[i]*w/W;
    double H1 = (tmp.basisDerivs_u[i]*W-tmp.basisValues[i]*Wderxi);
    bas.basisDerivs_u[i] = H1*w/(W*W);
    double H2 = (tmp.basisDerivs_v[i]*W-tmp.basisValues[i]*Wdereta);
    bas.basisDerivs_v[i] = H2*w/(W*W);
    double dH1dx = tmp.basisDerivs_uu[i]*W-tmp.basisValues[i]*Wder2xi;
    double dH2dx = tmp.basisDerivs_vv[i]*W-tmp.basisValues[i]*Wder2eta;
    double G1 = dH1dx*W-2*H1*Wderxi;
    bas.basisDerivs_uu[i] = G1*w/(W*W*W);
    double dH2dy = tmp.basisDerivs_vv[i]*W-tmp.basisValues[i]*Wder2eta;
    double G2 = dH2dy*W-2*H2*Wdereta;
    bas.basisDerivs_vv[i] = G2*w/(W*W*W);
    double dH1dy = tmp.basisDerivs_uv[i]*W + tmp.basisDerivs_u[i]*Wdereta -
                   tmp.basisDerivs_v[i]*Wderxi - tmp.basisValues[i]*Wderxieta;
    bas.basisDerivs_uv[i] = (dH1dy*W-2*H1*Wdereta)*w/(W*W*W);
    double d2H1dx2  =   tmp.basisDerivs_uuu[i]*W + tmp.basisDerivs_uu[i]*Wderxi
                      - tmp.basisDerivs_u[i]*Wder2xi - tmp.basisValues[i]*Wder3xi;
    double d2H1dxdy =   tmp.basisDerivs_uuv[i]*W + tmp.basisDerivs_uu[i]*Wdereta
                      - tmp.basisDerivs_v[i]*Wder2xi - tmp.basisValues[i]*Wder2xieta;
    double d2H2dy2  =   tmp.basisDerivs_vvv[i]*W + tmp.basisDerivs_vv[i]*Wdereta
                      - tmp.basisDerivs_v[i]*Wder2eta - tmp.basisValues[i]*Wder3eta;
    double d2H2dxdy =   tmp.basisDerivs_uvv[i]*W + tmp.basisDerivs_vv[i]*Wderxi
                      - tmp.basisDerivs_u[i]*Wder2eta - tmp.basisValues[i]*Wderxi2eta;

    double dG1dx = d2H1dx2*W + dH1dx*Wderxi - 2*dH1dx*Wderxi -2*H1*Wder2xi;
    double dG1dy = d2H1dxdy*W + dH1dx*Wdereta - 2*dH1dy*Wderxi - 2*H1*Wderxieta;
    double dG2dx = d2H2dxdy*W + dH2dy*Wderxi - 2*dH2dx*Wdereta -2*H2*Wderxieta;
    double dG2dy = d2H2dy2*W + dH2dy*Wdereta - 2*dH2dy*Wdereta - 2*H2*Wder2eta;

    bas.basisDerivs_uuu[i] = (dG1dx*W-3*G1*Wderxi)*w/(W*W*W*W);
    bas.basisDerivs_vvv[i] = (dG2dy*W-3*G2*Wdereta)*w/(W*W*W*W);
    bas.basisDerivs_uuv[i] = (dG1dy*W-3*G1*Wdereta)*w/(W*W*W*W);
    bas.basisDerivs_uvv[i] = (dG2dx*W-3*G2*Wderxi)*w/(W*W*W*W);
  }
}


LR::LRSplineSurface* ASMu2Dnurbs::createLRfromTensor ()
{
  if (tensorspline)
  {
    if (tensorspline->rational())
      lrspline.reset(createFromTensor(tensorspline));
    else
      lrspline.reset(new LR::LRSplineSurface(tensorspline));
    delete tensorspline;
    tensorspline = nullptr;
  }

  return lrspline.get();
}


bool ASMu2Dnurbs::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return false;
  }
#endif

  const LR::Element* el = lrspline->getElement(iel-1);
  X.resize(nsd,el->nBasisFunctions());

  int n = 1;
  for (LR::Basisfunction* b : el->support()) {
    X.fillColumn(n,&(*b->cp()));
    for (int j = 1; j < b->dim(); ++j)
      X(j,n) /= b->cp(b->dim()-1);
    ++n;
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}
