// $Id$
//==============================================================================
//!
//! \file SAMpatchPara.h
//!
//! \date Sep 6 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for distributed models.
//!
//==============================================================================

#ifndef _SAM_PATCH_PARA_H
#define _SAM_PATCH_PARA_H

#include "SAMpatch.h"
#include "PETScSupport.h"
#include <map>

class ProcessAdm;


/*!
  \brief This is a sub-class of SAMpatch with support for parallel processing.
*/

class SAMpatchPETSc : public SAMpatch
{
  typedef std::vector<IntVec> IntMat; //!< General integer matrix

public:
  //! \brief The constructor initializes the \a l2gn array.
  //! \param[in] g2ln Global-to-local node numbers for this processor
  //! \param[in] padm Parallel process administrator
  SAMpatchPETSc(const std::map<int,int>& g2ln, const ProcessAdm& padm);
  //! \brief The destructor destroys the index set arrays.
  virtual ~SAMpatchPETSc();

  //! \brief Allocates the dynamic arrays and populates them with data.
  //! \param[in] model All spline patches in the model
  //! \param[in] numNod Total number of unique nodes in the model
  virtual bool init(const std::vector<ASMbase*>& model, int numNod = 0);

  //! \brief Computes the dot-product of two vectors.
  //! \param[in] x First vector in dot-product
  //! \param[in] y Second vector in dot-product
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of dot-product
  //!
  //! \details Both vectors are defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real dot(const Vector& x, const Vector& y, char dofType = 'D') const;

  //! \brief Computes the L2-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of L2-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normL2(const Vector& x, char dofType = 'D') const;

  //! \brief Computes the inf-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param comp Local nodal DOF on input, index of the largest value on output
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of inf-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normInf(const Vector& x, size_t& comp, char dofType = 'D') const;

  //! \brief Returns the patch vector.
  const std::vector<ASMbase*> getPatches() const { return patch; }

private:
  //! \brief Setup a parallel index set for a given dofType
  void setupIS(char dofType) const;

  // Parameters for parallel computing
  int    nProc;      //!< Number of processes
  int    nleq;       //!< Number of equations for this processor
  int    nnodGlob;   //!< Number of global nodes;
  IntVec ghostNodes; //!< Indices for the ghost nodes
  IntVec l2gn;       //!< Local-to-global node numbers for this processor
  int    ieqmin;     //!< Minium equation number
  int    ieqmax;     //!< Maximun equation number

  const ProcessAdm& adm; //!< Parallel process administrator

  // For domain decomposition preconditioner
  std::vector<ASMbase*> patch; //!< The spline patches

  //! \brief Struct holding dof vector info
  struct DofIS {
    IS local;        //!< Local index set for dof type
    IS global;       //!< Global index set for dof type
    bool scatterCreated = false;
    VecScatter ctx;  //!< Scatterer
  };
  mutable std::map<char, DofIS> dofIS; //!< Map of dof type scatter info
};

#endif
