// $Id$
//==============================================================================
//!
//! \file ASMLRSplineGo.h
//!
//! \date 11 Apr 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for FE assembly drivers using Go LR B-splines.
//!
//==============================================================================

#ifndef _ASM_LR_SPLINE_GO_H
#define _ASM_LR_SPLINE_GO_H

#include "ASMbase.h"
#include "ASMunstruct.h"
#include "ThreadGroups.h"

#include "GoTools/geometry/BsplineBasis.h"

namespace Go {
  class LRSplineSurface;
}

namespace GoLR //! Utilities for LR-splines.
{
  class Basisfunction;

  //! \brief Expands the basis coefficients of an LR-spline object.
  //! \param basis The spline object to extend
  //! \param[in] v The vector to append to the basis coefficients
  //! \param[in] nf Number of fields in the given vector
  //! \param[in] ofs Offset in vector
  int extendControlPoints(Go::LRSplineSurface* basis, const Vector& v,
                          int nf, int ofs = 0);

  //! \brief Contracts the basis coefficients of an LR-spline object.
  //! \param basis The spline object to contact
  //! \param[out] v Vector containing the extracted basis coefficients
  //! \param[in] nf Number of fields in the given vector
  //! \param[in] ofs Offset in vector
  void contractControlPoints(Go::LRSplineSurface* basis, Vector& v,
                             int nf, int ofs = 0);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[in] spline The LR-spline object to get parameter values for
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] d Parameter direction (0,1,2)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] iel Element index
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  void getGaussPointParameters(const Go::LRSplineSurface* spline, RealArray& uGP,
                               int d, int nGauss, int iel, const double* xi);

  //! \brief Generates thread groups for a LR-spline mesh.
  //! \param[out] threadGroups The generated thread groups
  //! \param[in] lr The LR-spline to generate thread groups for
  //! \param[in] addConstraints If given, additional constraint bases
  void generateThreadGroups(ThreadGroups& threadGroups,
                            const Go::LRSplineSurface* lr,
                            const std::vector<Go::LRSplineSurface*>& addConstraints = {});
}


/*!
  \brief Base class for LR B-spline FE assembly drivers.
  \details This class contains methods common for unstructured spline patches.
*/

class ASMLRSplineGo : public ASMbase, public ASMunstruct
{
protected:
  //! \brief The constructor sets the space dimensions.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMLRSplineGo(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Special copy constructor for sharing of FE data.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  ASMLRSplineGo(const ASMLRSplineGo& patch, unsigned char n_f);

public:
  //! \brief Checks if the patch is empty.
  bool empty() const override { return geomB == nullptr; }

  //! \brief Returns parameter values and node numbers of the domain corners.
  bool getParameterDomain(Real2DMat&, IntVec*) const override;

  //! \brief Returns a const pointer to refinement basis.
  const Go::LRSplineSurface* getRefinementBasis() const { return refB.get(); }

  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  bool refine(const LR::RefineData& prm, Vectors& sol) override;

  using ASMbase::evalSolution;
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual Go::LRSplineSurface* evalSolution(const IntegrandBase& integrand) const = 0;

  //! \brief Returns a Bezier basis of order \a p.
  static Go::BsplineBasis getBezierBasis(int p,
                                         double start = -1.0, double end = 1.0);

  //! \brief Returns a list of basis functions having support on given elements.
  IntVec getFunctionsForElements(const IntVec& elements,
                                 bool globalId = false) const;
  //! \brief Returns a list of basis functions having support on given elements.
  void getFunctionsForElements(IntSet& functions, const IntVec& elements,
                               bool globalId = true) const;

  //! \brief Sort basis functions based on local knot vectors.
  // static void Sort(int u, int v, int orient,
  //                  std::vector<Go::Basisfunction*>& functions);

  //! \brief Returns all boundary functions that are covered by the given nodes.
  //! \param[in] nodes Set of (0-based) patch local node IDs
  //! \return 0-based node IDs for boundary functions whose support is
  //! completely covered by the union of the support of the input nodes
  IntVec getBoundaryCovered(const IntSet& nodes) const override;

  //! \brief Returns all functions whose support overlap with the input nodes.
  //! \param[in] nodes List of (0-based) patch local node IDs
  //! (typically requested by adaptive refinement)
  //! \param[in] dir 3-bit binary mask on which parameter directions are allowed
  //! to grow; i.e. bin(011)=dec(3) allows u-direction and v-direction to grow,
  //! default is bin(111)=dec(7) all directions
  //! \return 0-based node IDs for functions with overlapping support with
  //! the ones in boundary
  IntVec getOverlappingNodes(const IntSet& nodes, int dir = 7) const;

  //! \brief Returns all functions whose support overlap with the input node.
  //! \param[in] node 0-based patch local node ID
  //! \param[in] dir 3-bit binary mask for which parameter directions can grow
  //! \return 0-based node IDs for functions with overlapping support
  IntVec getOverlappingNodes(int node, int dir = 7) const
  {
    return this->getOverlappingNodes(IntSet(&node,(&node)+1),dir);
  }

  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVars(const Go::LRSplineSurface* oldBasis,
                                   const RealArray& oldVar, RealArray& newVar,
                                   int nGauss) const = 0;
  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVarsN(const Go::LRSplineSurface* oldBasis,
                                    const RealArray& oldVar, RealArray& newVar,
                                    int nGauss) const = 0;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferCntrlPtVars(const Go::LRSplineSurface* oldBasis,
                                   RealArray& newVar, int nGauss) const = 0;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Control point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] nf Number of field components
  bool transferCntrlPtVars(Go::LRSplineSurface* oldBasis,
                           const RealArray& oldVar, RealArray& newVar,
                           int nGauss, int nf = 1) const;

  //! \brief Finds the node that is closest to the given point \b X.
  std::pair<size_t,double> findClosestNode(const Vec3& X) const override;

  //! \brief Returns the coordinates of the element center.
  Vec3 getElementCenter(int iel) const override;

  //! \brief Computes the total number of integration points in this patch.
  void getNoIntPoints(size_t& nPt, size_t& nIPt) override;

  //! \brief Swaps between the first and second projection basis.
  void swapProjectionBasis() override;

protected:
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param lrspline The spline to perform adaptation for
  bool doRefine(const LR::RefineData& prm, Go::LRSplineSurface* lrspline);

  using ASMbase::evalPoint;
  //! \brief Evaluates the geometry at a specified point.
  virtual int evalPoint(int iel, const double* param, Vec3& X) const = 0;

  //! \brief Santity check thread groups.
  //! \param[in] groups The generated thread groups
  //! \param[in] bases The bases to check for
  //! \param[in] threadBasis The LRSpline the element groups are derived from
  static bool checkThreadGroups(const IntMat& groups,
                                const std::vector<const Go::LRSplineSurface*>& bases,
                                const Go::LRSplineSurface* threadBasis);

  //! \brief Analyze and print thread group statistics.
  //! \param[in] groups The generated thread groups
  static void analyzeThreadGroups(const IntMat& groups);

  std::shared_ptr<Go::LRSplineSurface> geomB;  //!< Pointer to spline object of the geometry basis
  std::shared_ptr<Go::LRSplineSurface> projB;  //!< Pointer to spline object of the projection basis
  std::shared_ptr<Go::LRSplineSurface> projB2; //!< Pointer to spline object of the secondary projection basis
  std::shared_ptr<Go::LRSplineSurface> refB;   //!< Pointer to spline object of the refinement basis

  ThreadGroups threadGroups; //!< Element groups for multi-threaded assembly
  ThreadGroups projThreadGroups; //!< Element groups for multi-threaded assembly - projection basis
  ThreadGroups proj2ThreadGroups; //!< Element groups for multi-threaded assembly - second projection basis
};

#endif
