// $Id$
//==============================================================================
//!
//! \file LRNURBSFields2D.C
//!
//! \date Apr 23 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR NURBS-based finite element vector fields in 2D.
//!
//==============================================================================

#include "LRNURBSFields2D.h"
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



LRNURBSFields2D::LRNURBSFields2D (const ASMu2D* patch,
                                  const RealArray& v, char nbasis,
                                  int nnf, const char* name)
  : Fields(name), basis(patch->getBasis(nbasis)), surf(patch->getSurface())
{
  nno = basis->nBasisFunctions();
  nelm = basis->nElements();

  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  auto vit = v.begin()+ofs;

  if (nnf == 0)
    nnf = 2;
  nf = nnf;
  int nfc = patch->getNoFields(nbasis);
  values.resize(nno*nf);
  int ndof = nfc*nno;
  auto end = v.size() > ofs+ndof ? vit+ndof : v.end();
  if (nfc == nf)
    std::copy(vit,end,values.begin());
  else
    for (size_t i = 0; i < nno && vit != end; ++i, vit += nfc)
      for (size_t j = 0; j < nf; ++j)
        values[nf*i+j] = *(vit+j);

  basis->generateIDs();
}


LRNURBSFields2D::LRNURBSFields2D (const LR::LRSplineSurface* srf,
                                  const RealArray& v, int cmp, const char* name)
  : Fields(name), basis(srf), surf(srf)
{
  values = v;
  nf = cmp;
}


bool LRNURBSFields2D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool LRNURBSFields2D::valueFE (const ItgPoint& x, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(x.u,x.v);
  auto elm = basis->getElement(iel);

  Go::BasisPtsSf spline;
  computeBasis(x.u,x.v,spline,iel,basis);

  // Evaluate the solution field at the given point
  Matrix Vnod(nf, elm->nBasisFunctions());
  size_t i = 1;
  for (const LR::Basisfunction* f : elm->support()) {
    for (size_t j = 1; j <= nf; ++j)
      Vnod(j,i) = values(f->getId()*nf+j);
    ++i;
  }

  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}


bool LRNURBSFields2D::valueCoor (const Vec4& x, Vector& vals) const
{
  if (x.u)
    return this->valueFE(ItgPoint(x.u[0],x.u[1]),vals);

  return false;
}


bool LRNURBSFields2D::gradFE (const ItgPoint& x, Matrix& grad) const
{
  return false;
}


bool LRNURBSFields2D::hessianFE (const ItgPoint& x, Matrix3D& H) const
{
  return false;
}
