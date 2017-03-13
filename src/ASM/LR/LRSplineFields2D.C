// $Id$
//==============================================================================
//!
//! \file LRSplineFields2D.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for LR spline-based finite element vector fields in 2D.
//!
//==============================================================================

#include "LRSplineFields2D.h"
#include "ASMu2D.h"
#include "FiniteElement.h"
#include "CoordinateMapping.h"
#include "Utilities.h"
#include "Vec3.h"

#include "LRSpline/LRSplineSurface.h"


LRSplineFields2D::LRSplineFields2D (const ASMu2D* patch,
                                    const RealArray& v, char nbasis,
                                    const char* name)
  : Fields(name), basis(patch->getBasis(nbasis)), surf(patch->getSurface())
{
  nno = basis->nBasisFunctions();
  nelm = basis->nElements();

  auto vit = v.begin();
  size_t ofs = 0;
  for (char i = 1; i < nbasis; ++i)
    ofs += patch->getNoNodes(i)*patch->getNoFields(i);
  vit += ofs;

  // Ensure the values array has compatible length, pad with zeros if necessary
  nf = 2;
  values.reserve(nno*nf);
  if (patch->getNoFields(nbasis) != 2) {
    RealArray::const_iterator end = v.size() > nno*nf+ofs ? vit+nno*nf : v.end();
    for (size_t i = 0; i < nno && vit != end; ++i) {
      for (size_t k = 0; k < nf && vit != end; ++k, ++vit)
        values.push_back(*vit);
      vit += patch->getNoFields(nbasis)-nf;
    }
  } else {
    RealArray::const_iterator end = v.size() > nno*nf+ofs ? vit+nno*nf : v.end();
    std::copy(vit,end,std::back_inserter(values));
  }
  values.resize(nno*nf);
}


bool LRSplineFields2D::valueNode (size_t node, Vector& vals) const
{
  if (node < 1 || node > nno) return false;

  vals.resize(nf);
  vals.fill(values.ptr()+(node-1)*nf);
  return true;
}


bool LRSplineFields2D::valueFE (const FiniteElement& fe, Vector& vals) const
{
  if (!basis) return false;

  // Evaluate the basis functions at the given point
  int iel = basis->getElementContaining(fe.u,fe.v);
  auto elm = basis->getElement(basis->getElementContaining(fe.u,fe.v));

  Go::BasisPtsSf spline;
  basis->computeBasis(fe.u,fe.v,spline,iel);

  // Evaluate the solution field at the given point

  Matrix Vnod(2, elm->nBasisFunctions());
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 2; ++j)
      Vnod(j,i) = values((*it)->getId()*2+j);

  Vnod.multiply(spline.basisValues,vals); // vals = Vnod * basisValues

  return true;
}


bool LRSplineFields2D::valueCoor (const Vec3& x, Vector& vals) const
{
  // Not implemented yet
  return false;
}


bool LRSplineFields2D::gradFE (const FiniteElement& fe, Matrix& grad) const
{
  if (!basis) return false;
  if (!surf)  return false;

  // Evaluate the basis functions at the given point
  int iel = surf->getElementContaining(fe.u,fe.v);
  auto elm = surf->getElement(iel);
  Go::BasisDerivsSf spline;
  surf->computeBasis(fe.u,fe.v,spline,iel);

  const size_t nen = elm->nBasisFunctions();

  Matrix dNdu(nen,2), dNdX;
  for (size_t n = 1; n <= nen; n++)
  {
    dNdu(n,1) = spline.basisDerivs_u[n-1];
    dNdu(n,2) = spline.basisDerivs_v[n-1];
  }

  Matrix Xnod(2, nen), Jac;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 2; ++j)
      Xnod(j, i) = (*it)->cp(j-1);

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  Matrix Vnod;
  if (basis != surf)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
    int iel = basis->getElementContaining(fe.u,fe.v);
    auto belm = basis->getElement(iel);
    basis->computeBasis(fe.u,fe.v,spline,iel);

    const size_t nbf = belm->nBasisFunctions();
    dNdu.resize(nbf,2);
    for (size_t n = 1; n <= nbf; n++)
    {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }
    dNdX.multiply(dNdu,Jac); // dNdX = dNdu * Jac

    size_t i = 1;
    Vnod.resize(nbf,2);
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= 2; ++j)
        Vnod(i,j) = values((*it)->getId()*2+j);
  } else {
    size_t i = 1;
    Vnod.resize(2, nen);
    for (auto it  = elm->constSupportBegin();
              it != elm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= 2; ++j)
        Vnod(j,i) = values((*it)->getId()*2+j);
  }

  grad.multiply(Vnod,dNdX); // grad = Xnod * dNdX

  return true;
}


bool LRSplineFields2D::hessianFE(const FiniteElement& fe, Matrix3D& H) const
{
  if (!basis) return false;
  if (!surf)  return false;

  int iel = surf->getElementContaining(fe.u,fe.v);
  auto elm = surf->getElement(iel);
  const size_t nen = elm->nBasisFunctions();

  // Evaluate the basis functions at the given point
  Go::BasisDerivsSf  spline;
  Go::BasisDerivsSf2 spline2;
  Matrix3D d2Ndu2;
  Matrix dNdu, dNdX;
  Matrix Xnod(2, nen), Jac, Vnod;
  size_t i = 1;
  for (auto it  = elm->constSupportBegin();
            it != elm->constSupportEnd(); ++it, ++i)
    for (size_t j = 1; j <= 2; ++j)
      Xnod(j, i) = (*it)->cp(j-1);

  if (surf == basis) {
    surf->computeBasis(fe.u,fe.v,spline2,iel);

    dNdu.resize(nen,2);
    d2Ndu2.resize(nen,2,2);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
    }

    size_t i = 1;
    Vnod.resize(2, nen);
    for (auto it  = elm->constSupportBegin();
        it != elm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= 2; ++j)
        Vnod(i,j) = values((*it)->getId()*2+j);
  }
  else {
    surf->computeBasis(fe.u,fe.v,spline,iel);

    dNdu.resize(nen,2);
    for (size_t n = 1; n <= nen; n++) {
      dNdu(n,1) = spline.basisDerivs_u[n-1];
      dNdu(n,2) = spline.basisDerivs_v[n-1];
    }
  }

  // Evaluate the Jacobian inverse
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Evaluate the gradient of the solution field at the given point
  if (basis != surf)
  {
    // Mixed formulation, the solution uses a different basis than the geometry
    int iel = basis->getElementContaining(fe.u,fe.v);
    auto belm = basis->getElement(iel);
    basis->computeBasis(fe.u,fe.v,spline2,iel);

    const size_t nbf = belm->nBasisFunctions();
    dNdu.resize(nbf,2);
    d2Ndu2.resize(nbf,2,2);
    for (size_t n = 1; n <= nbf; n++) {
      dNdu(n,1) = spline2.basisDerivs_u[n-1];
      dNdu(n,2) = spline2.basisDerivs_v[n-1];
      d2Ndu2(n,1,1) = spline2.basisDerivs_uu[n-1];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline2.basisDerivs_uv[n-1];
      d2Ndu2(n,2,2) = spline2.basisDerivs_vv[n-1];
    }

    size_t i = 1;
    Vnod.resize(2, nbf);
    for (auto it  = belm->constSupportBegin();
              it != belm->constSupportEnd(); ++it, ++i)
      for (size_t j = 1; j <= 2; ++j)
        Vnod(j,i) = values((*it)->getId()*2+j);
  }

  return H.multiply(Vnod,d2Ndu2);
}


bool LRSplineFields2D::gradCoor (const Vec3& x, Matrix& grad) const
{
  // Not implemented yet
  return false;
}
