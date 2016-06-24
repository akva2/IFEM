// $Id$
//==============================================================================
//!
//! \file GlobalIntegral.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for classes representing integrated quantities.
//!
//==============================================================================

#ifndef _GLOBAL_INTEGRAL_H
#define _GLOBAL_INTEGRAL_H

class LocalIntegral;


/*!
  \brief Abstract base class representing a system level integrated quantity.
*/

class GlobalIntegral
{
public:
  //! \brief The default constructor.
  GlobalIntegral() {}
  //! \brief Empty destructor.
  virtual ~GlobalIntegral() {}

  //! \brief Initializes the integrated quantity to zero.
  virtual void initialize(bool) {}
  //! \brief Finalizes the integrated quantity after element assembly.
  virtual bool finalize(bool) { return true; }

  //! \brief Adds a LocalIntegral object into a corresponding global object.
  //! \param[in] elmObj The local integral object to add into \a *this.
  //! \param[in] elmId Global number of the element associated with elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId) { return true; }

  //! \brief Adds a LocalIntegral object into a corresponding global object.
  //! \param[in] elmObj The local integral object to add into \a *this.
  //! \param[in] elmId1 Global number of the first element associated with elmObj
  //! \param[in] elmId2 Global number of the second element associated with elmObj
  virtual bool assemble(const LocalIntegral* elmObj,
                        int elmId1, int elmId2) { return true; }
};

#endif
