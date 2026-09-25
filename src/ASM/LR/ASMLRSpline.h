// $Id$
//==============================================================================
//!
//! \file ASMLRSpline.h
//!
//! \date December 2010
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for FE assembly drivers using LR B-splines.
//!
//==============================================================================

#ifndef _ASM_LR_SPLINE_H
#define _ASM_LR_SPLINE_H

#include "ASMbase.h"
#include "ASMunstruct.h"
#include "ThreadGroups.h"

#include "GoTools/geometry/BsplineBasis.h"


namespace LR //! Utilities for LR-splines.
{
  class Basisfunction;
  class LRSpline;

  //! \brief Expands the basis coefficients of an LR-spline object.
  //! \param basis The spline object to extend
  //! \param[in] v The vector to append to the basis coefficients
  //! \param[in] nf Number of fields in the given vector
  //! \param[in] ofs Offset in vector
  int extendControlPoints(LRSpline* basis, const Vector& v,
                          int nf, int ofs = 0);

  //! \brief Contracts the basis coefficients of an LR-spline object.
  //! \param basis The spline object to contact
  //! \param[out] v Vector containing the extracted basis coefficients
  //! \param[in] nf Number of fields in the given vector
  //! \param[in] ofs Offset in vector
  void contractControlPoints(LRSpline* basis, Vector& v,
                             int nf, int ofs = 0);

  //! \brief Extracts parameter values of the Gauss points in one direction.
  //! \param[in] spline The LR-spline object to get parameter values for
  //! \param[out] uGP Parameter values in given direction for all points
  //! \param[in] d Parameter direction (0,1,2)
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] iel Element index
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  void getGaussPointParameters(const LRSpline* spline, RealArray& uGP,
                               int d, int nGauss, int iel, const double* xi);

  //! \brief Generates thread groups for a LR-spline mesh.
  //! \param[out] threadGroups The generated thread groups
  //! \param[in] lr The LR-spline to generate thread groups for
  //! \param[in] addConstraints If given, additional constraint bases
  void generateThreadGroups(ThreadGroups& threadGroups,
                            const LRSpline* lr,
                            const std::vector<LRSpline*>& addConstraints = {});

  //! \brief Createss the matrix of nodal point correspondance for a LR-spline.
  //! \param[in] basis LR-spline to get nodal point correspondance for
  //! \param[out] result Matrix of nodal point correspondance for the elements
  void createMNPC(const LR::LRSpline* basis, IntMat& result);
}


/*!
  \brief Base class for LR B-spline FE assembly drivers.
  \details This class contains methods common for unstructured spline patches.
*/

class ASMLRSpline : public ASMbase, public ASMunstruct
{
protected:
  //! \brief The constructor sets the space dimensions.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMLRSpline(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Special copy constructor for sharing of FE data.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  ASMLRSpline(const ASMLRSpline& patch, unsigned char n_f);

public:
  //! \brief Empty destructor.
  virtual ~ASMLRSpline() {}

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geomB == nullptr; }

  //! \brief Returns parameter values and node numbers of the domain corners.
  virtual bool getParameterDomain(Real2DMat&, IntVec*) const;

  //! \brief Returns a const pointer to refinement basis.
  const LR::LRSpline* getRefinementBasis() const { return refB.get(); }

  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol);

  //! \brief Stores solution vectors as extra control point dimensions.
  //! \param[in] sol Solution vectors to store
  //! \param[out] nf Number of field components stored for each vector
  //!
  //! \details Refining an LR-spline is knot insertion, which interpolates the
  //! control points onto the refined mesh. Solution vectors are therefore
  //! carried across a refinement by storing them as extra control point
  //! dimensions, and extracting them again afterwards. The two halves are
  //! exposed separately such that a multi-patch driver can keep the vectors
  //! stored while it makes the patch meshes conform with each other.
  virtual bool packSolution(const Vectors& sol, IntVec& nf);
  //! \brief Extracts solution vectors from the extra control point dimensions.
  //! \param sol Solution vectors to extract, resized to the current mesh
  //! \param[in] nf Number of field components stored for each vector
  virtual void unpackSolution(Vectors& sol, const IntVec& nf);

  //! \brief Refines the mesh, leaving any stored solution vectors alone.
  //! \param[in] prm Input data used to control the mesh refinement
  virtual bool refineMesh(const LR::RefineData& prm);

  //! \brief Propagates the current mesh onto a separate projection basis.
  //! \details Does nothing unless this patch has a projection basis of its
  //! own. Insertion of an existing knot line is a no-op, so this may be
  //! invoked whenever the mesh has settled, also more than once.
  virtual bool refineProjectionBasis() { return true; }

  //! \brief Groups the basis functions into lines along a parameter direction.
  //! \param[in] dir Parameter direction the lines run along, 1-based
  //! \param[out] lines Patch-local node numbers of each line, 1-based
  //!
  //! \details Two functions belong to the same line when their local knot
  //! vectors agree in every direction but \a dir, and a line comes out
  //! ordered along the direction it runs in. On a tensor mesh that is a mesh
  //! line, which is what a smoother wants to solve along on an anisotropic
  //! mesh; on a locally refined one it is what is left of a mesh line after
  //! the refinement, so the grouping loses accuracy rather than meaning.
  bool getLineDofs(int dir, std::vector<IntVec>& lines) const;

  //! \brief Checks the basis of this patch for linear independence.
  //! \return \e false if the basis is linearly dependent, or inconclusive
  bool checkLinearIndependence() const;

  using ASMbase::evalSolution;
  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual LR::LRSpline* evalSolution(const IntegrandBase& integrand) const = 0;

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
  static void Sort(int u, int v, int orient,
                   std::vector<LR::Basisfunction*>& functions);

  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVars(const LR::LRSpline* oldBasis,
                                   const RealArray& oldVar, RealArray& newVar,
                                   int nGauss) const = 0;
  //! \brief Transfers Gauss point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Gauss point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferGaussPtVarsN(const LR::LRSpline* oldBasis,
                                    const RealArray& oldVar, RealArray& newVar,
                                    int nGauss) const = 0;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  virtual bool transferCntrlPtVars(const LR::LRSpline* oldBasis,
                                   RealArray& newVar, int nGauss) const = 0;
  //! \brief Transfers control point variables from old basis to this patch.
  //! \param[in] oldBasis The LR-spline basis to transfer from
  //! \param[in] oldVar Control point variables associated with \a oldBasis
  //! \param[out] newVar Gauss point variables associated with this patch
  //! \param[in] nGauss Number of Gauss points along a knot-span
  //! \param[in] nf Number of field components
  bool transferCntrlPtVars(LR::LRSpline* oldBasis,
                           const RealArray& oldVar, RealArray& newVar,
                           int nGauss, int nf = 1) const;

  //! \brief Finds the node that is closest to the given point \b X.
  virtual std::pair<size_t,double> findClosestNode(const Vec3& X) const;

  //! \brief Returns the coordinates of the element center.
  virtual Vec3 getElementCenter(int iel) const;

  //! \brief Computes the total number of integration points in this patch.
  virtual void getNoIntPoints(size_t& nPt, size_t& nIPt);

  //! \brief Swaps between the first and second projection basis.
  virtual void swapProjectionBasis();

protected:
  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param lrspline The spline to perform adaptation for
  bool doRefine(const LR::RefineData& prm, LR::LRSpline* lrspline);

  //! \brief Updates the patch after its mesh has been changed.
  //! \details Regenerates the spline function IDs and discards the FE data
  //! established for the mesh as it was before the change. A patch with more
  //! than one basis also has to bring the other bases along here, since only
  //! the refinement basis is matched with that of a neighbouring patch.
  virtual void meshUpdated();

  using ASMbase::evalPoint;
  //! \brief Evaluates the geometry at a specified point.
  virtual int evalPoint(int iel, const double* param, Vec3& X) const = 0;

  //! \brief Santity check thread groups.
  //! \param[in] groups The generated thread groups
  //! \param[in] bases The bases to check for
  //! \param[in] threadBasis The LRSpline the element groups are derived from
  static bool checkThreadGroups(const IntMat& groups,
                                const std::vector<const LR::LRSpline*>& bases,
                                const LR::LRSpline* threadBasis);

  std::shared_ptr<LR::LRSpline> geomB;  //!< Pointer to spline object of the geometry basis
  std::shared_ptr<LR::LRSpline> projB;  //!< Pointer to spline object of the projection basis
  std::shared_ptr<LR::LRSpline> projB2; //!< Pointer to spline object of the secondary projection basis
  std::shared_ptr<LR::LRSpline> refB;   //!< Pointer to spline object of the refinement basis

  ThreadGroups threadGroups; //!< Element groups for multi-threaded assembly
  ThreadGroups projThreadGroups; //!< Element groups for multi-threaded assembly - projection basis
  ThreadGroups proj2ThreadGroups; //!< Element groups for multi-threaded assembly - second projection basis
};

#endif
