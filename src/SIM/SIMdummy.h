// $Id$
//==============================================================================
//!
//! \file SIMdummy.h
//!
//! \date May 27 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dummy %SIM template class for drivers not associated with a FE model.
//!
//==============================================================================

#ifndef _SIM_DUMMY_H_
#define _SIM_DUMMY_H_

#include <vector>
#include <iostream>

class ASMbase;
class IntegrandBase;
class ModelGenerator;
class TiXmlElement;
namespace ASM { struct Interface; }


/*!
  \brief Template %SIM class with some dummy implementations.
  \details This class only implements dummy versions for the pure virtual
  virtual methods of the base class SIMbase, and can be used as a base for
  simulator drivers that do not require any FE model.
*/

template<class Base> class SIMdummy : public Base
{
protected:
  //! \brief Default constructor.
  explicit SIMdummy(IntegrandBase* p = nullptr) : Base(p) {}
  //! \brief Empty destructor.
  virtual ~SIMdummy() {}
public:
  //! \brief Returns the number of parameter dimensions in the model.
  virtual unsigned short int getNoParamDim() const { return 0; }
  //! \brief Creates the computational FEM model from the spline patches.
  virtual bool createFEMmodel(char) { return false; }
  //! \brief Element-element connectivities.
  virtual std::vector<std::vector<int>> getElmConnectivities() const
  { return std::vector<std::vector<int>>(); }
protected:
  //! \brief Parses a dimension-specific subelement of the \a geometry XML-tag.
  virtual bool parseGeometryDimTag(const TiXmlElement*) { return false; }
  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  virtual bool addConstraint(int,int,int,int,int,int&,char) { return false; }
  //! \brief Preprocesses the result sampling points.
  virtual void preprocessResultPoints() {}
  //! \brief Creates a model generator.
  virtual ModelGenerator* getModelGenerator(const TiXmlElement*) const
  { return nullptr; }
  //! \brief Reads a patch from given input stream.
  virtual ASMbase* readPatch(std::istream&,int,
                             const std::vector<unsigned char>&,
                             const char*) const { return nullptr; }
  //! \brief Connects two patches.
  virtual bool connectPatches(const ASM::Interface&,bool) { return false; }
};

#endif
