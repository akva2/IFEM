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
#include "FiniteElement.h"
#include "CoordinateMapping.h"


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

  int wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;
  RealArray w(MNPC[iel].size());
  for (size_t i = 0; i < w.size(); i++)
    w[i] = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);

  if (derivs < 1) {
    Matrix B;
    B.outer_product(Nu,Nv);
    fe.N = C*static_cast<const Vector&>(B);
    double W = fe.N.dot(w);
    for (size_t i = 0; i < fe.N.size(); i++)
      fe.N[i] *= w[i]/W;

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

    fe.N = C*B;
    Vector dNxi  = C*Bu;
    Vector dNeta = C*Bv;

    double W       = fe.N.dot(w);
    double Wderxi  = dNxi.dot(w);
    double Wdereta = dNeta.dot(w);

    Matrix dNdu(w.size(),2);
    for (size_t i = 1; i <= fe.N.size(); i++) {
      fe.N(i)  *= w[i-1]/W;
      dNdu(i,1) = (dNxi(i)*W  - fe.N(i)*Wderxi )*w[i-1]/(W*W);
      dNdu(i,2) = (dNeta(i)*W - fe.N(i)*Wdereta)*w[i-1]/(W*W);
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
  int wi = spline->getBasisfunction(MNPC[iel][0])->dim()-1;
  Vector w(tmp.basisValues.size());
  for (size_t i = 0; i < w.size(); i++)
    w[i] = spline->getBasisfunction(MNPC[iel][i])->cp(wi);

  double W = w.dot(tmp.basisValues);

  bas.preparePts(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;
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
  int wi = spline->getBasisfunction(MNPC[iel][0])->dim()-1;
  Vector w(tmp.basisValues.size());
  for (size_t i = 0; i < w.size(); i++)
    w[i] = spline->getBasisfunction(MNPC[iel][i])->cp(wi);

  double W  = w.dot(tmp.basisValues);
  double Wx = w.dot(tmp.basisDerivs_u);
  double Wy = w.dot(tmp.basisDerivs_v);

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
  {
    bas.basisValues[i]   =  tmp.basisValues[i]*w[i]/W;
    bas.basisDerivs_u[i] = (tmp.basisDerivs_u[i]*w[i] - bas.basisValues[i]*Wx)/W;
    bas.basisDerivs_v[i] = (tmp.basisDerivs_v[i]*w[i] - bas.basisValues[i]*Wy)/W;
  }
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf2& bas, int iel) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  Go::BasisDerivsSf2 tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  int wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;
  Vector w(tmp.basisValues.size());
  for (size_t i = 0; i < w.size(); i++)
    w[i] = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);

  double W   = w.dot(tmp.basisValues);
  double Wx  = w.dot(tmp.basisDerivs_u);
  double Wy  = w.dot(tmp.basisDerivs_v);
  double Wxx = w.dot(tmp.basisDerivs_uu);
  double Wyy = w.dot(tmp.basisDerivs_vv);
  double Wxy = w.dot(tmp.basisDerivs_uv);

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); ++i)
  {
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;

    double H1 = (tmp.basisDerivs_u[i]*W - tmp.basisValues[i]*Wx);
    double H2 = (tmp.basisDerivs_v[i]*W - tmp.basisValues[i]*Wy);
    bas.basisDerivs_u[i] = H1*w[i]/(W*W);
    bas.basisDerivs_v[i] = H2*w[i]/(W*W);

    double H1x = tmp.basisDerivs_uu[i]*W - tmp.basisValues[i]*Wxx;
    double H2y = tmp.basisDerivs_vv[i]*W - tmp.basisValues[i]*Wyy;
    double H1y = tmp.basisDerivs_uv[i]*W - tmp.basisValues[i]*Wxy
               + tmp.basisDerivs_u[i]*Wy - tmp.basisDerivs_v[i]*Wx;
    bas.basisDerivs_uu[i] = (H1x*W - 2.0*H1*Wx)*w[i]/(W*W*W);
    bas.basisDerivs_vv[i] = (H2y*W - 2.0*H2*Wy)*w[i]/(W*W*W);
    bas.basisDerivs_uv[i] = (H1y*W - 2.0*H1*Wy)*w[i]/(W*W*W);
  }
}


void ASMu2Dnurbs::computeBasis(double u, double v,
                               Go::BasisDerivsSf3& bas, int iel) const
{
  if (noNurbs)
    return this->ASMu2D::computeBasis(u,v,bas,iel);

  Go::BasisDerivsSf3 tmp;
  lrspline->computeBasis(u,v,tmp,iel);
  int wi = lrspline->getBasisfunction(MNPC[iel][0])->dim()-1;
  Vector w(tmp.basisValues.size());
  for (size_t i = 0; i < w.size(); i++)
    w[i] = lrspline->getBasisfunction(MNPC[iel][i])->cp(wi);

  double W    = w.dot(tmp.basisValues);
  double Wx   = w.dot(tmp.basisDerivs_u);
  double Wy   = w.dot(tmp.basisDerivs_v);
  double Wxx  = w.dot(tmp.basisDerivs_uu);
  double Wyy  = w.dot(tmp.basisDerivs_vv);
  double Wxy  = w.dot(tmp.basisDerivs_uv);
  double Wxxx = w.dot(tmp.basisDerivs_uuu);
  double Wyyy = w.dot(tmp.basisDerivs_vvv);
  double Wxxy = w.dot(tmp.basisDerivs_uuv);
  double Wxyy = w.dot(tmp.basisDerivs_uvv);

  bas.prepareDerivs(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
  {
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;

    double H1 = (tmp.basisDerivs_u[i]*W - tmp.basisValues[i]*Wx);
    double H2 = (tmp.basisDerivs_v[i]*W - tmp.basisValues[i]*Wy);
    bas.basisDerivs_u[i] = H1*w[i]/(W*W);
    bas.basisDerivs_v[i] = H2*w[i]/(W*W);

    double H1x = tmp.basisDerivs_uu[i]*W - tmp.basisValues[i]*Wxx;
    double H2y = tmp.basisDerivs_vv[i]*W - tmp.basisValues[i]*Wyy;
    double H1y = tmp.basisDerivs_uv[i]*W - tmp.basisValues[i]*Wxy
               + tmp.basisDerivs_u[i]*Wy - tmp.basisDerivs_v[i]*Wx;
    double H2x = tmp.basisDerivs_uv[i]*W - tmp.basisValues[i]*Wxy
               + tmp.basisDerivs_v[i]*Wx - tmp.basisDerivs_u[i]*Wy;
    double G1  = H1x*W - 2.0*H1*Wx;
    double G2  = H2y*W - 2.0*H2*Wy;
    bas.basisDerivs_uu[i] = G1*w[i]/(W*W*W);
    bas.basisDerivs_vv[i] = G2*w[i]/(W*W*W);
    bas.basisDerivs_uv[i] = (H1y*W - 2.0*H1*Wy)*w[i]/(W*W*W);

    double H1xx = tmp.basisDerivs_uuu[i]*W + tmp.basisDerivs_uu[i]*Wx
                - tmp.basisDerivs_u[i]*Wxx - tmp.basisValues[i]*Wxxx;
    double H2yy = tmp.basisDerivs_vvv[i]*W + tmp.basisDerivs_vv[i]*Wy
                - tmp.basisDerivs_v[i]*Wyy - tmp.basisValues[i]*Wyyy;
    double H1xy = tmp.basisDerivs_uuv[i]*W + tmp.basisDerivs_uu[i]*Wy
                - tmp.basisDerivs_v[i]*Wxx - tmp.basisValues[i]*Wxxy;
    double H2xy = tmp.basisDerivs_uvv[i]*W + tmp.basisDerivs_vv[i]*Wx
                - tmp.basisDerivs_u[i]*Wyy - tmp.basisValues[i]*Wxyy;

    double G1x = H1xx*W + H1x*Wx - 2.0*H1x*Wx - 2.0*H1*Wxx;
    double G1y = H1xy*W + H1x*Wy - 2.0*H1y*Wx - 2.0*H1*Wxy;
    double G2x = H2xy*W + H2y*Wx - 2.0*H2x*Wy - 2.0*H2*Wxy;
    double G2y = H2yy*W + H2y*Wy - 2.0*H2y*Wy - 2.0*H2*Wyy;

    bas.basisDerivs_uuu[i] = (G1x*W - 3.0*G1*Wx)*w[i]/(W*W*W*W);
    bas.basisDerivs_vvv[i] = (G2y*W - 3.0*G2*Wy)*w[i]/(W*W*W*W);
    bas.basisDerivs_uuv[i] = (G1y*W - 3.0*G1*Wy)*w[i]/(W*W*W*W);
    bas.basisDerivs_uvv[i] = (G2x*W - 3.0*G2*Wx)*w[i]/(W*W*W*W);
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

/* This is not strictly correct as we return rational spline coefficients,
   rather than coordinates from this method. But since these are used as
   coefficients, and not coordinates at caller sites, we do this for simplicity.
   We should consider introducing ASMbase::getElementCoefficients to lessen confusion.
*/
bool ASMu2Dnurbs::getElementCoordinates (Matrix& X, int iel) const
{
  if (noNurbs)
    return this->ASMu2D::getElementCoordinates(X,iel);

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
    double weight = b->cp(b->dim()-1);
#ifdef SP_DEBUG
    if (weight <= 0.0)
    {
      std::cerr <<" *** ASMu2Dnurbs::getElementCoordinates: Zero weight for"
                <<" node "<< n <<" in element "<< iel << std::endl;
      return false;
    }
#endif
    for (int j = 1; j <= nsd; j++)
      X(j,n) = b->cp(j-1) / weight;
    ++n;
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}
