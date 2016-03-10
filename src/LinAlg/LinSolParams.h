// $Id$
//==============================================================================
//!
//! \file LinSolParams.h
//!
//! \date Jan 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//! \details Includes linear solver method, preconditioner
//! and convergence criteria.
//!
//==============================================================================

#ifndef _LINSOLPARAMS_H
#define _LINSOLPARAMS_H

#include "PETScSupport.h"

#include <array>
#include <iostream>
#include <string>
#include <vector>


class ParamMap {
public:
  void addValue(const std::string& key, const std::string& value);
  std::string getStringValue(const std::string& key) const;
  int getIntValue(const std::string& key) const;
  double getDoubleValue(const std::string& key) const;

  bool hasValue(const std::string& key) const;
private:
  std::map<std::string, std::string> values;
};

#ifdef HAS_PETSC

typedef std::vector<int>         IntVec;       //!< Integer vector
typedef std::vector<IntVec>      IntMat;       //!< Integer matrix
typedef std::vector<std::string> StringVec;    //!< String vector
typedef std::vector<StringVec>   StringMat;    //!< String matrix
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix

//! \brief Null-space for matrix operator
enum NullSpace { NONE, CONSTANT, RIGID_BODY };

//! \brief Schur preconditioner methods
enum SchurPrec { SIMPLE, MSIMPLER, PCD };

#endif

class TiXmlElement;


#ifdef HAS_PETSC
  #define BLANK_IF_NO_PETSC(a) a
#else
  #define BLANK_IF_NO_PETSC(a) ""
#endif


/*!
  \brief Class for linear solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class LinSolParams
{
public:
  //! \brief Default constructor.
  LinSolParams();

  //! \brief Read linear solver parameters from XML document.
  bool read(const TiXmlElement* elem);

  //! \brief Linear solver settings for a block of the linear system
  struct BlockParams {

    BlockParams();
      basis(1), comps(0)
    {}

    //! \brief Settings for a directional smoother
    struct DirSmoother {
      int order; //!< Ordering of DOFs
      std::string type; //!< Directional smoother types

      DirSmoother(int o, const std::string& t) : order(o), type(t) {}
    };

    //! \brief Read settings from XML block
    //! \param[in] child XML block
    bool read(const TiXmlElement* child, const std::string& prefix="");

    size_t basis; //!< Basis for block
    size_t comps; //!< Components from basis (1, 2, 3, 12, 13, 23, 123, ..., 0 = all)
    ValueMap settings;
  };

  //! \brief Number of blocks in matrix system
  size_t getNoBlocks() const { return blocks.size(); }

  //! \brief Obtain settings for a given block
  const BlockParams& getBlock(size_t i) const { return blocks[i]; }

private:
  std::vector<BlockParams> blocks; //!< Parameters for each block
  ValueMap settings; //!< Generic settings
};

#endif
