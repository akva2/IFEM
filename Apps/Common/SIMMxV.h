//==============================================================================
//!
//! \file SIMMxV.h
//!
//! \date Nov 28 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Use a simulator as a matrix-vector product for iterative solvers.
//!
//==============================================================================
#ifndef SIM_MXV_H_
#define SIM_MXV_H_

#ifdef HAS_PETSC

#include "PETScMatrix.h"


/*! \brief Template for applying the matrix assembled in a SIM to a vector 
           for use in iterative solvers.
*/

template<class Sim>
class SIMMxV : public PETScMxV {
public:
  //! \brief Default constructor.
  //! \param sim The simulator to wrap
  SIMMxV(Sim& sim) : S1(sim)
  {
    this->S1.setMxV(this);
  }

  //! \brief Evaluate the matrix-vector product y = A*y
  bool evalMxV(Vec& x, Vec& y) override
  {
    PETScMatrix* A = static_cast<PETScMatrix*>(S1.getAlgEqSystem()->getMatrix(0));
    MatMult(A->getBlockMatrices()[0], x, y);
    return true;
  }

  //! \brief Set a matrix-free preconditioner.
  void setPC(PETScPC* pc)
  {
    static_cast<PETScMatrix*>(S1.getAlgEqSystem()->getMatrix(0))->setPC(pc);
  }

  //! \brief Solve at time level.
  bool solveStep(TimeStep& tp)
  {
    return this->S1.solveStep(tp);
  }

  //! \brief The matrix used for building the preconditioner.
  //! \details Defaults to the system matrix
  bool evalPC(Mat& P) override
  {
    PETScMatrix*A = static_cast<PETScMatrix*>(S1.getAlgEqSystem()->getMatrix(0));
    MatDuplicate(A->getMatrix(), MAT_COPY_VALUES, &P);
    return true;
  }
protected:
  Sim& S1; //!< Reference to wrapped simulator
};

#endif

#endif
