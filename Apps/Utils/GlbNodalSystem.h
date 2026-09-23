//==============================================================================
//!
//! \file GlbNodalSystem.h
//!
//! \date Sep 19 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Global integral for assembly in the unconstrained nodal ordering.
//!
//==============================================================================

#ifndef _GLB_NODAL_SYSTEM_H_
#define _GLB_NODAL_SYSTEM_H_

#include "GlobalIntegral.h"
#include "MatVec.h"

#include <vector>

using IntVec = std::vector<int>; //!< General integer vector

class SAM;
class SparseMatrix;


/*!
  \brief Global integral assembling into the unconstrained nodal system.

  \details Optimal control formulations need \f$L_2\f$ inner products over a
  control space which is not subject to the Dirichlet conditions of the state
  problem. Those cannot be obtained from the constrained system matrices held
  by the AlgEqSystem of a simulator, as the constrained degrees of freedom have
  been eliminated from them.

  This class therefore assembles element matrices, element vectors and element
  scalars directly in the degree-of-freedom ordering of the model, bypassing
  the constraint handling performed by SAM. It handles any number of degrees of
  freedom per node, including the mixed case where the number varies from node
  to node, and orders the element degrees of freedom exactly as SAM does.

  The constructor pre-computes the sparsity pattern of the coefficient matrix,
  so that the element loop only updates entries that already exist. Since the
  ASM classes integrate the elements in node-disjoint thread groups whenever
  the global integral is not thread safe, the matrix and vector assembly is
  then race free without any locking. This is the same approach as in
  ASMbase::globalL2projection(), which assembles the patch-local
  \f$L_2\f$-projection matrix. The global scalars are the exception, as every
  element contributes to all of them, and they are updated atomically.
*/

class GlbNodalSystem : public GlobalIntegral
{
public:
  //! \brief Where one of the element matrices is to be assembled.
  struct MatrixTarget
  {
    SparseMatrix* mat;        //!< The matrix to assemble into
    size_t        idx = 0;    //!< Index of the element matrix to assemble
    bool   transposed = false; //!< If \e true, assemble the transpose
  };

  //! \brief The constructor pre-computes the sparsity patterns.
  //! \param[in] sam Data for FE assembly management
  //! \param[in] mats Where each of the element matrices is to be assembled
  //! \param vec Vector to assemble the first element vector into, if any
  //! \param scl Scalars to accumulate the element scalars into, if any
  GlbNodalSystem(const SAM& sam, const std::vector<MatrixTarget>& mats,
                 Vector* vec, RealArray* scl);
  //! \brief Convenience constructor for a single coefficient matrix.
  //! \param[in] sam Data for FE assembly management
  //! \param mat Matrix to assemble the first element matrix into, if any
  //! \param vec Vector to assemble the first element vector into, if any
  //! \param scl Scalars to accumulate the element scalars into, if any
  GlbNodalSystem(const SAM& sam, SparseMatrix* mat,
                 Vector* vec, RealArray* scl);
  //! \brief Empty destructor.
  virtual ~GlbNodalSystem() {}

  //! \brief Initializes the integrated quantities to zero.
  void initialize(char) override;

  //! \brief Adds an element contribution into the global quantities.
  //! \param[in] elmObj The element quantities to add
  //! \param[in] elmId Global number of the element associated with \a elmObj
  bool assemble(const LocalIntegral* elmObj, int elmId) override;

private:
  //! \brief Finds the degree-of-freedom numbers of an element.
  //! \param[out] meen One-based global DOF numbers, zero for absent nodes
  //! \param[in] elmId Global element number
  bool elementDofs(IntVec& meen, int elmId) const;

  //! \brief Prints an element size mismatch error message.
  //! \param[in] what Which element quantity that did not match
  //! \param[in] elmId Global element number
  //! \param[in] nedof Number of degrees of freedom on the element
  static bool sizeError(const char* what, int elmId, size_t nedof);

  //! \brief Pre-computes the sparsity patterns of the target matrices.
  void preAssemble();

  const SAM& mySam; //!< Data for FE assembly management

  std::vector<MatrixTarget> myMats; //!< Nodal coefficient matrices
  Vector*    myVec; //!< Nodal right-hand-side vector, or \e nullptr
  RealArray* myScl; //!< Global scalar quantities, or \e nullptr
};

#endif
