// $Id$
//==============================================================================
//!
//! \file LRNURBSField2D.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element scalar field in 2D.
//!
//==============================================================================

#include "LRNURBSField2D.h"
#include "LRSpline/LRSplineSurface.h"

#include "ASMu2D.h"
#include "ItgPoint.h"
#include "Vec3.h"

namespace {

void computeBasis(double u, double v,
                Go::BasisPtsSf& bas, int iel,
                const LR::LRSplineSurface* spline)
{
  const LR::Element* el = spline->getElement(iel);

  Go::BasisPtsSf tmp;
  spline->computeBasis(u,v,tmp,iel);
  Vector w; w.reserve(tmp.basisValues.size());
  for (const LR::Basisfunction* func : el->support())
    w.push_back(func->cp(func->dim()-1));

  double W = w.dot(tmp.basisValues);

  bas.preparePts(u, v, 0, -1, tmp.basisValues.size());
  for (size_t i = 0; i < tmp.basisValues.size(); i++)
    bas.basisValues[i] = tmp.basisValues[i]*w[i]/W;
}

}


LRNURBSField2D::LRNURBSField2D (const ASMu2D* patch,
                                const RealArray& v, char nbasis,
                                char cmp, const char* name)
  : FieldBase(name), basis(patch->getBasis(nbasis)), surf(patch->getSurface())
{
  nno = basis->nBasisFunctions();
  nelm = basis->nElements();

  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  auto vit = v.begin()+ofs;

  values.resize(nno);
  int nf = patch->getNoFields(nbasis);
  int ndof = nf > 1 && cmp > 0 ? nf*nno : nno;
  auto end = v.size() > ofs+ndof ? vit+ndof : v.end();
  if (nf == 1 || cmp == 0)
    std::copy(vit,end,values.begin());
  else
    for (size_t i = 0; i < nno && vit != end; ++i, vit += nf)
      values[i] = *(vit+cmp-1);

  basis->generateIDs();
}


LRNURBSField2D::LRNURBSField2D (const LR::LRSplineSurface* srf,
                                const RealArray& v, const char* name)
  : FieldBase(name), basis(srf), surf(srf)
{
  values = v;
}


double LRNURBSField2D::valueNode (size_t node) const
{
  return node > 0 && node <= nno ? values(node) : 0.0;
}


double LRNURBSField2D::valueFE (const ItgPoint& x) const
{
  if (!basis) return 0.0;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(x.u,x.v);
  auto elm = basis->getElement(iel);

  Go::BasisPtsSf spline;
  computeBasis(x.u,x.v,spline,iel,basis);

  Vector Vnod;
  Vnod.reserve(elm->nBasisFunctions());
  for (const LR::Basisfunction* f : elm->support())
    Vnod.push_back(values(f->getId()+1));

  return Vnod.dot(spline.basisValues);
}


double LRNURBSField2D::valueCoor (const Vec4& x) const
{
  // Just produce a segmentation fault, if invoked without the parameters
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1]));

  std::cerr << "** LRNURBSField2D::valueCoor: "
            << "not implemented without parameters\n";

  return false;
}


bool LRNURBSField2D::gradFE (const ItgPoint& x, Vector& grad) const
{
  return false;
}


bool LRNURBSField2D::hessianFE (const ItgPoint& x, Matrix& H) const
{
  return false;
}
