// $Id$
//==============================================================================
//!
//! \file PETScMatrix.h
//!
//! \date Jan 15 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the system matrix in PETSc format with interface
//! to PETSc routines for assembling and solving linear equation systems.
//!
//==============================================================================

#ifndef _PETSC_MATRIX_H
#define _PETSC_MATRIX_H

#include "SystemMatrix.h"
#include "SparseMatrix.h"
#include "PETScSupport.h"
#include "LinAlgenums.h"

class LinSolParams;

typedef std::vector<PetscInt>    PetscIntVec;  //!< PETSc integer vector
typedef std::vector<PetscIntVec> PetscIntMat;  //!< PETSc integer matrix
typedef std::vector<PetscReal>   PetscRealVec; //!< PETSc real vector
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix


/*!
  \brief Class for representing the system vector in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScVector : public StdVector
{
public:
#ifdef HAS_PETSC
  //! \brief Constructor creating an empty vector.
  PETScVector(const ProcessAdm& padm);
  //! \brief Constructor creating a vector of length \a n.
  PETScVector(const ProcessAdm& padm, size_t n);
  //! \brief Constructor creating a vector from an array.
  PETScVector(const ProcessAdm& padm, const Real* values, size_t n);
  //! \brief Copy constructor.
  PETScVector(const PETScVector& vec);
  //! \brief Destructor.
  virtual ~PETScVector();
#endif

  //! \brief Returns the vector type.
  virtual Type getType() const { return PETSC; }

#ifdef HAS_PETSC
  //! \brief Initializes the vector to a given scalar value.
  virtual void init(Real value = Real(0));

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n);

  //! \brief Begins communication step needed in parallel vector assembly.
  //! \details Must be called together with endAssembly after vector assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool beginAssembly();

  //! \brief Ends communication step needed in parallel vector assembly.
  //! \details Must be called together with beginAssembly after vector assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool endAssembly();

  //! \brief L1-norm of vector.
  virtual Real L1norm() const;

  //! \brief L2-norm of vector.
  virtual Real L2norm() const;

  //! \brief Linfinity-norm of vector.
  virtual Real Linfnorm() const;

  //! \brief Returns the PETSc vector (for assignment).
  virtual Vec& getVector() { return x; }
  //! \brief Returns the PETSc vector (for read access).
  virtual const Vec& getVector() const { return x; }

  //! \brief Return associated process administrator
  const ProcessAdm& getAdm() const { return adm; }
protected:
  Vec x;                  //!< The actual PETSc vector
  const ProcessAdm& adm;  //!< Process administrator

#else // dummy implementation when PETSc is not included
  virtual SystemVector* copy() const { return 0; }
  virtual size_t dim() const { return 0; }
  virtual void redim(size_t) {}
  virtual Real* getPtr() { return 0; }
  virtual const Real* getRef() const { return 0; }
  virtual void init(Real = Real(0)) {}
  virtual void mult(Real) {}
  virtual Real L1norm() const { return Real(0); }
  virtual Real L2norm() const { return Real(0); }
  virtual Real Linfnorm() const { return Real(0); }
#endif
};


/*!
  \brief Class for representing the system matrix in PETSc format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScMatrix : public SparseMatrix
{
public:
#ifdef HAS_PETSC
  //! \brief Constructor.
  PETScMatrix(const ProcessAdm& padm, const LinSolParams& spar,
              LinAlg::LinearSystemType ltype);
  //! \brief Copy constructor.
  PETScMatrix(const PETScMatrix& A);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~PETScMatrix();
#else
  //! \brief Constructor.
  PETScMatrix(const ProcessAdm&, const LinSolParams&) {}
#endif

  //! \brief Returns the matrix type.
  virtual Type getType() const { return PETSC; }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int = 1) const { return 0; }

#ifdef HAS_PETSC
  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new PETScMatrix(*this); }

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! The PETSc data structures are initialized and the all symbolic operations
  //! that are needed before the actual assembly can start are performed.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam, bool);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  //! \brief Begins communication step needed in parallel matrix assembly.
  //! \details Must be called together with endAssembly after matrix assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool beginAssembly();
  //! \brief Ends communication step needed in parallel matrix assembly.
  //! \details Must be called together with beginAssembly after matrix assembly
  //! is completed on each processor and before the linear system is solved.
  virtual bool endAssembly();

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C) const;

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector& B, bool newLHS, Real*);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] B Right-hand-side vector
  //! \param[out] x Solution vector
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(const SystemVector& B, SystemVector& x, bool newLHS);

  //! \brief Solves a generalized symmetric-definite eigenproblem.
  //! \details The eigenproblem is assumed to be on the form
  //! \b A \b x = \f$\lambda\f$ \b B \b x where \b A ( = \a *this ) and \b B
  //! both are assumed to be symmetric and \b B also to be positive definite.
  //! The eigenproblem is solved by the SLEPc library subroutine \a EPSSolve.
  //! \sa SLEPc library documentation.
  //! \param B Symmetric and positive definite mass matrix.
  //! \param[out] eigVal Computed eigenvalues
  //! \param[out] eigVec Computed eigenvectors stored column by column
  //! \param[in] nev The number of eigenvalues and eigenvectors to compute
  //! \param[in] shift Eigenvalue shift
  //! \param[in] iop Option telling whether to factorize matrix \a A or \b B.
  virtual bool solveEig(PETScMatrix& B, RealArray& eigVal, Matrix& eigVec,
			int nev, Real shift = Real(0), int iop = 1);

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const;

  //! \brief Returns the PETSc matrix (for assignment).
  virtual Mat& getMatrix() { return A; }
  //! \brief Returns the PETSc matrix (for read access).
  virtual const Mat& getMatrix() const { return A; }

  //! \brief Set the linear solver parameters (solver type, preconditioner, tolerances).
  //! \param[in] P Preconditioner  matrix (ignored here)
  //! \param[in] Pb Preconditioner vector (ignored here)
  //! \return True on success
  virtual bool setParameters(PETScMatrix* P = nullptr, PETScVector* Pb = nullptr);
protected:
  //! \brief Solve a linear system
  bool solve(const Vec& b, Vec& x, bool newLHS, bool knoll);

  Mat                 A;               //!< The actual PETSc matrix
  KSP                 ksp;             //!< Linear equation solver
  MatNullSpace*       nsp;             //!< Null-space of linear operator
  const ProcessAdm&   adm;             //!< Process administrator
  const LinSolParams& solParams;       //!< Linear solver parameters
  bool                setParams;       //!< If linear solver parameters are set
  PetscInt            ISsize;          //!< Number of index sets/elements
  PetscIntMat         locSubdDofs;     //!< Degrees of freedom for unique subdomains
  PetscIntMat         subdDofs;        //!< Degrees of freedom for subdomains
  PetscRealVec        coords;          //!< Coordinates of local nodes (x0,y0,z0,x1,y1,...)
  ISMat               dirIndexSet;     //!< Direction ordering
  int                 nLinSolves;      //!< Number of linear solves
  LinAlg::LinearSystemType linsysType; //!< Linear system type
  IS glob2LocEq = nullptr; //!< Index set for global-to-local equations.

#else // dummy implementation when PETSc is not included
  virtual SystemMatrix* copy() const { return 0; }
  virtual void init() {}
  virtual void initAssembly(const SAM&, bool) {}
  virtual bool assemble(const Matrix&, const SAM&, int) { return false; }
  virtual bool assemble(const Matrix&, const SAM&,
			SystemVector&, int) { return false; }
#endif
};


class PETScDOFVectorOps : public DOFVectorOps
{
public:
  //! \brief The constructor initializes the \a l2gn array.
  //! \param[in] g2ln Global-to-local node numbers for this processor
  //! \param[in] padm Parallel process administrator
  PETScDOFVectorOps(const ProcessAdm& padm);
  //! \brief The destructor destroys the index set arrays.
  virtual ~PETScDOFVectorOps();

  //! \brief Computes the dot-product of two vectors.
  //! \param[in] x First vector in dot-product
  //! \param[in] y Second vector in dot-product
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of dot-product
  //!
  //! \details Both vectors are defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real dot(const Vector& x, const Vector& y, char dofType, const SAM& sam) const;

  //! \brief Computes the L2-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of L2-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normL2(const Vector& x, char dofType, const SAM& sam) const;

  //! \brief Computes the inf-norm of a vector.
  //! \param[in] x Vector for norm computation
  //! \param comp Local nodal DOF on input, index of the largest value on output
  //! \param[in] dofType Only consider nodes of this type (for mixed methods)
  //! \return Value of inf-norm
  //!
  //! \details The vector is defined over all nodes in the patch, i.e.
  //! for parallel vectors the ghost entries are also included.
  virtual Real normInf(Real value) const;

private:
  //! \brief Setup a parallel index set for a given dofType
  void setupIS(char dofType, const SAM& sam) const;

#ifdef HAS_PETSC
  //! \brief Struct holding dof vector info
  struct DofIS {
    IS local;        //!< Local index set for dof type
    IS global;       //!< Global index set for dof type
    bool scatterCreated = false;
    VecScatter ctx;  //!< Scatterer
  };
  mutable std::map<char, DofIS> dofIS; //!< Map of dof type scatter info
#endif

  const ProcessAdm& adm; //<! Process administrator and domain decomposition
};


#ifdef HAS_PETSC
//! \brief Matrix-vector product
PETScVector operator*(const SystemMatrix& A, const PETScVector& b);

//! \brief Solve linear system
PETScVector operator/(const SystemMatrix& A, const PETScVector& b);
#endif

#endif
