// $Id$
//==============================================================================
//!
//! \file OptimizationDriver.h
//!
//! \brief Reduced-space optimization driver for adjoint-based optimal control.
//!
//==============================================================================

#ifndef _OPTIMIZATION_DRIVER_H_
#define _OPTIMIZATION_DRIVER_H_

#include "IFEM.h"
#include "MatVec.h"
#include "Utilities.h"

#include "tinyxml2.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


/*!
  \brief Reduced-space optimizer for PDE-constrained optimal control problems.

  \details This class holds the optimization algorithm alone, decoupled from
  how the sequence of discretizations is arrived at. It solves the problem in
  the reduced space, using adjoint sensitivity analysis to evaluate the
  gradient of the objective with respect to the control, and can be
  instantiated over any simulator implementing the ISolverOpt interface. It
  supports

  - Inexact Newton-CG, which is mesh-independent as long as the simulator
    provides the reduced Hessian through ISolverOpt::applyHessian(),
  - L-BFGS with a limited memory of secant pairs,
  - steepest descent,

  all three with an Armijo backtracking line search and optional projection
  onto box constraints on the control.

  All inner products and norms are taken in the metric supplied by the
  simulator through ISolverOpt::applyMass(), which for a finite element control
  space is the \f$L_2\f$ inner product. Using this metric rather than the
  Euclidean one on the DOF values is what makes the convergence history
  independent of the mesh resolution.

  Separating the algorithm from the driver is what lets the same optimizer be
  run once on a fixed mesh, or repeatedly inside an adaptive loop where each
  cycle is warm started from the control transferred off the previous mesh.
*/

template<class T1>
class OptimizationDriver
{
public:
  //! \brief Optimization method.
  enum class Method
  {
    NEWTON_CG,       //!< Inexact Newton with an inner conjugate gradient solve
    LBFGS,           //!< Limited memory BFGS
    GRADIENT_DESCENT //!< Steepest descent
  };

  //! \brief The constructor initializes the reference to the simulator.
  explicit OptimizationDriver(T1& s1) : S1(s1) {}

  //! \brief Defines the optimization method to use.
  void setMethod(Method m) { method = m; }
  //! \brief Defines the absolute gradient norm tolerance.
  void setTolerance(double t) { tol = t; }
  //! \brief Defines the maximum number of optimization iterations.
  void setMaxIterations(int it) { maxIter = it; }
  //! \brief Clears the reference gradient norm of the relative criterion.
  //! \details Call this to have the next optimize() re-anchor the relative
  //! convergence criterion to its own starting point.
  void resetReference() { gNormRef = 0.0; }

  //! \brief Defines box constraints on the control.
  void setBounds(double low, double high)
  {
    lowerBound = low;
    upperBound = high;
    hasBounds = true;
  }

  //! \brief Runs the optimization on the current discretization.
  //! \param F Control vector, holding the initial guess on entry and the
  //! optimal control on return
  //! \param[out] U Optimal state
  //! \param[out] P Adjoint solution at the optimum
  //! \return Zero if the convergence criterion was met, 12 if the iteration
  //! ran out of steps or the line search gave up, and a code in the range 4-9
  //! if one of the simulator operations failed
  int optimize(Vector& F, Vector& U, Vector& P)
  {
    if (hasBounds)
      this->projectBounds(F);

    nStagnant = 0;

    // Initial forward, adjoint and gradient evaluation
    Vector G;
    if (!this->S1.solveForward(F,U))
      return this->error("Initial forward solve failed",4);

    double J = this->S1.computeObjective(F,U);

    if (!this->S1.solveAdjoint(U,F,P))
      return this->error("Initial adjoint solve failed",5);

    if (!this->S1.computeGradient(F,P,G))
      return this->error("Initial gradient evaluation failed",6);

    double gNorm = this->criticality(F,G);

    // The relative criterion is anchored to the gradient norm of the first
    // discretization the optimizer is run on. In an adaptive loop each cycle
    // after the first starts from the control transferred off the previous
    // mesh, whose gradient is already small, and anchoring to that would ask
    // for a further eight orders of reduction below the attainable floor
    if (gNormRef <= 0.0)
      gNormRef = gNorm;

    const double gNorm0 = gNormRef;

    IFEM::cout <<"\n >>> Optimization history <<<"
               <<"\n Method: "<< this->methodName()
               <<"\n Convergence criterion: ||dJ/df|| < "<< tol
               <<" or ||dJ/df||/||dJ/df||_0 < "<< rtol <<"\n\n"
               <<"  Iter        Objective       ||dJ/df||   Step size   CG its"
               <<"\n-------------------------------------------------"
               <<"-----------------"<< std::endl;

    this->printIteration(0,J,gNorm,0.0,-1);

    bool converged = gNorm <= tol;
    int  totalCGIts = 0;
    int  iter = 0;

    // Secant pairs for the L-BFGS two-loop recursion
    std::deque<Vector> sHist, yHist;
    std::deque<double> rhoHist;

    while (!converged && iter < maxIter)
    {
      ++iter;

      // Variables that are held at a bound are excluded from the second order
      // model and take a steepest descent step instead, which the projection
      // then clips back onto the bound. This is Bertsekas' two-metric
      // projection, and is what keeps a Newton or quasi-Newton direction a
      // descent direction once some of the bounds are active
      std::vector<bool> atBound;
      this->activeSet(F,G,gNorm,atBound);

      Vector Gfree(G);
      this->freeOnly(Gfree,atBound);

      Vector d;
      int innerIts = 0;
      switch (method) {
      case Method::NEWTON_CG:
        this->newtonStep(Gfree,gNorm,atBound,d,innerIts);
        totalCGIts += innerIts;
        break;
      case Method::LBFGS:
        this->lbfgsStep(Gfree,sHist,yHist,rhoHist,d);
        break;
      default:
        d = Gfree;
        d *= -1.0;
      }

      for (size_t i = 0; i < atBound.size(); i++)
        if (atBound[i])
          d[i] = -G[i];

      // Ensure we have a descent direction, fall back to steepest descent
      double dirDeriv = this->dot(G,d);
      if (dirDeriv >= 0.0)
      {
        d = G;
        d *= -1.0;
        dirDeriv = this->dot(G,d);
      }

      // Armijo backtracking line search. The sufficient decrease is measured
      // against the increment that is actually taken, which for a projected
      // step is shorter than the trial step, and reduces to step*dirDeriv
      // when the control is unconstrained
      Vector Fnew, Unew, sk;
      double step = 1.0, Jnew = J;
      bool haveStep = false;
      for (int ls = 0; ls < maxLineSearch && step >= stol; ls++, step *= 0.5)
      {
        Fnew = F;
        Fnew.add(d,step);
        if (hasBounds)
          this->projectBounds(Fnew);

        if (!this->S1.solveForward(Fnew,Unew))
          continue; // The state equation could not be solved for this step

        sk = Fnew;
        sk.add(F,-1.0);

        Jnew = this->S1.computeObjective(Fnew,Unew);
        if (Jnew <= J + armijoC1*this->dot(G,sk))
        {
          haveStep = true;
          break;
        }
      }

      if (!haveStep)
      {
        IFEM::cout <<"  ** Line search failed, the step size dropped below "
                   << stol <<".\n     Stopping with the last accepted iterate."
                   << std::endl;
        --iter;
        break;
      }

      Vector yk;
      if (method == Method::LBFGS)
      {
        yk = G;
        yk *= -1.0;
      }

      F = Fnew;
      U = Unew;
      J = Jnew;

      if (!this->S1.solveAdjoint(U,F,P))
        return this->error("Adjoint solve failed",8);

      if (!this->S1.computeGradient(F,P,G))
        return this->error("Gradient evaluation failed",9);

      if (method == Method::LBFGS)
      {
        yk.add(G);
        double sy = this->dot(sk,yk);
        if (sy > curvTol) // skip the update if the curvature condition fails
        {
          if (sHist.size() >= static_cast<size_t>(lbfgsM))
          {
            sHist.pop_front();
            yHist.pop_front();
            rhoHist.pop_front();
          }
          sHist.push_back(sk);
          yHist.push_back(yk);
          rhoHist.push_back(1.0/sy);
        }
      }

      double gPrev = gNorm;
      gNorm = this->criticality(F,G);
      this->printIteration(iter,J,gNorm,step,innerIts);

      converged = gNorm <= tol || (gNorm0 > 0.0 && gNorm <= rtol*gNorm0);

      // Stop if the gradient has stopped responding to the iteration. This is
      // the roundoff floor of the reduced gradient, below which further steps
      // only consume line search trials
      if (converged)
        nStagnant = 0;
      else if (gNorm < gPrev*(1.0-stagTol))
        nStagnant = 0;
      else if (++nStagnant >= maxStagnant)
      {
        IFEM::cout <<"  ** The gradient norm stagnated at "<< gNorm
                   <<", which is the attainable floor for this problem."
                   <<"\n     Stopping with the last accepted iterate."
                   << std::endl;
        break;
      }
    }

    IFEM::cout <<"--------------------------------------------------"
               <<"----------------\n"
               <<" Optimization "
               << (converged ? "converged in " : "stopped after ")
               << iter <<" iteration"<< (iter == 1 ? "" : "s") <<".";
    if (method == Method::NEWTON_CG)
      IFEM::cout <<"\n Total number of inner CG iterations: "<< totalCGIts;
    IFEM::cout <<"\n Final objective value: "<< J
               <<"\n Final gradient norm:   "<< gNorm << std::endl;

    return converged ? 0 : 12;
  }

  //! \brief Parses the optimization parameters from an XML element.
  //! \param[in] elem The \a optimalcontrol or \a adjoint element to parse
  bool parse(const tinyxml2::XMLElement* elem)
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement("optimizer");
    if (!child)
      return true; // The remaining tags are for the simulator, not the driver

    std::string mStr;
    if (utl::getAttribute(child,"method",mStr,true))
      this->setMethodFromString(mStr);

    utl::getAttribute(child,"tolerance",tol);
    utl::getAttribute(child,"rtolerance",rtol);
    utl::getAttribute(child,"max_iterations",maxIter);
    utl::getAttribute(child,"max_cg_iterations",maxCGIter);
    utl::getAttribute(child,"memory",lbfgsM);

    if (const tinyxml2::XMLElement* ls = child->FirstChildElement("linesearch"))
    {
      utl::getAttribute(ls,"c1",armijoC1);
      utl::getAttribute(ls,"min_step",stol);
    }

    if (const tinyxml2::XMLElement* bc = child->FirstChildElement("bounds"))
      if (utl::getAttribute(bc,"lower",lowerBound) |
          utl::getAttribute(bc,"upper",upperBound))
        hasBounds = true;

    IFEM::cout <<"\tOptimization method: "<< this->methodName()
               <<"\n\tGradient tolerance: "<< tol <<" (absolute), "
               << rtol <<" (relative)"
               <<"\n\tMaximum number of iterations: "<< maxIter;
    if (method == Method::NEWTON_CG)
      IFEM::cout <<"\n\tMaximum number of inner CG iterations: "<< maxCGIter;
    else if (method == Method::LBFGS)
      IFEM::cout <<"\n\tNumber of stored secant pairs: "<< lbfgsM;
    if (hasBounds)
      IFEM::cout <<"\n\tControl bounds: ["<< lowerBound <<","<< upperBound <<"]";
    IFEM::cout << std::endl;

    return true;
  }

private:
  T1& S1; //!< Reference to the simulator being optimized

  //! \brief Computes an inexact Newton step by an inner CG iteration.
  //! \param[in] G Current gradient, zeroed on the active set
  //! \param[in] gNorm Norm of the current gradient
  //! \param[in] atBound Marks the control variables held at a bound
  //! \param[out] d The resulting search direction
  //! \param[out] nIt Number of inner CG iterations performed
  //! \return Always \e true, the method degrades gracefully to a steepest
  //! descent step if the simulator provides no reduced Hessian
  //!
  //! \details The inner iteration solves \f${\bf H}{\bf d} = -{\bf G}\f$ by
  //! Steihaug-style conjugate gradients in the mass matrix metric. It is
  //! truncated when a direction of non-positive curvature is encountered, and
  //! the tolerance follows the Eisenstat-Walker rule to obtain superlinear
  //! convergence close to the optimum without over-solving far away from it.
  bool newtonStep(const Vector& G, double gNorm,
                  const std::vector<bool>& atBound, Vector& d, int& nIt)
  {
    d.resize(G.size(),true);

    Vector r(G);
    r *= -1.0;
    Vector p(r);

    double rr = this->dot(r,r);
    double eta = std::min(0.5,std::sqrt(std::max(gNorm,0.0)));
    double cgTol2 = std::max(eta*eta*rr,1.0e-30);

    for (nIt = 0; nIt < maxCGIter; )
    {
      Vector Hp;
      if (!this->S1.applyHessian(p,Hp))
      {
        if (nIt == 0)
        {
          // No second order information, take a steepest descent step instead
          IFEM::cout <<"  ** No reduced Hessian available,"
                     <<" using the steepest descent direction."<< std::endl;
          d = p;
        }
        return true;
      }

      ++nIt;
      this->freeOnly(Hp,atBound);
      double pHp = this->dot(p,Hp);
      if (pHp <= 0.0)
      {
        if (nIt == 1) d = p; // Steepest descent on the first iteration
        return true;
      }

      double alpha = rr/pHp;
      d.add(p,alpha);
      r.add(Hp,-alpha);

      double rrNew = this->dot(r,r);
      if (rrNew <= cgTol2)
        return true;

      p *= rrNew/rr;
      p.add(r);
      rr = rrNew;
    }

    return true;
  }

  //! \brief Computes an L-BFGS step through the two-loop recursion.
  //! \param[in] G Current gradient
  //! \param[in] sHist Stored control increments
  //! \param[in] yHist Stored gradient increments
  //! \param[in] rhoHist Stored inverse curvatures
  //! \param[out] d The resulting search direction
  void lbfgsStep(const Vector& G,
                 const std::deque<Vector>& sHist,
                 const std::deque<Vector>& yHist,
                 const std::deque<double>& rhoHist,
                 Vector& d)
  {
    const size_t k = sHist.size();

    Vector q(G);
    std::vector<double> alpha(k,0.0);
    for (size_t i = k; i-- > 0;)
    {
      alpha[i] = rhoHist[i]*this->dot(sHist[i],q);
      q.add(yHist[i],-alpha[i]);
    }

    if (k > 0)
    {
      // Scale the initial Hessian approximation by the Barzilai-Borwein factor
      double yy = this->dot(yHist.back(),yHist.back());
      if (yy > 0.0)
        q *= this->dot(sHist.back(),yHist.back())/yy;
    }

    for (size_t i = 0; i < k; i++)
      q.add(sHist[i],alpha[i] - rhoHist[i]*this->dot(yHist[i],q));

    d = q;
    d *= -1.0;
  }

  //! \brief Inner product in the metric defined by the simulator.
  double dot(const Vector& x, const Vector& y) const
  {
    Vector My;
    if (!this->S1.applyMass(y,My))
      return x*y; // No metric available, use the Euclidean inner product

    return x*My;
  }

  //! \brief Flags the control variables that are held at a bound.
  //! \param[in] F Current control
  //! \param[in] G Gradient at \a F
  //! \param[in] gNorm Criticality measure at \a F
  //! \param[out] atBound Marks the variables belonging to the active set
  //!
  //! \details A variable belongs to the active set when it sits on (or very
  //! near) a bound and the negative gradient pushes it further outside. The
  //! width of the "very near" band shrinks with the criticality measure, so
  //! that the active set is identified exactly in the limit.
  void activeSet(const Vector& F, const Vector& G, double gNorm,
                 std::vector<bool>& atBound) const
  {
    atBound.assign(F.size(),false);
    if (!hasBounds)
      return;

    double eps = std::min(bndTol,gNorm);
    for (size_t i = 0; i < F.size(); i++)
      atBound[i] = (F[i] <= lowerBound + eps && G[i] > 0.0) ||
                   (F[i] >= upperBound - eps && G[i] < 0.0);
  }

  //! \brief Zeroes the entries of a vector belonging to the active set.
  static void freeOnly(Vector& x, const std::vector<bool>& atBound)
  {
    for (size_t i = 0; i < x.size() && i < atBound.size(); i++)
      if (atBound[i])
        x[i] = 0.0;
  }

  //! \brief Returns the first order criticality measure of a control state.
  //! \param[in] F Current control
  //! \param[in] G Gradient at \a F
  //!
  //! \details Without box constraints this is simply the norm of the
  //! gradient. With box constraints it is the norm of the projected gradient
  //! \f${\bf f}-P({\bf f}-\nabla J)\f$, which vanishes exactly at a
  //! first order critical point of the constrained problem, whereas the
  //! gradient itself does not vanish where a bound is active.
  double criticality(const Vector& F, const Vector& G) const
  {
    if (!hasBounds)
      return this->norm(G);

    Vector R(F);
    R.add(G,-1.0);
    this->projectBounds(R);
    R.add(F,-1.0);

    return this->norm(R);
  }

  //! \brief Norm in the metric defined by the simulator.
  double norm(const Vector& x) const
  {
    double x2 = this->dot(x,x);
    return x2 > 0.0 ? std::sqrt(x2) : 0.0;
  }

  //! \brief Projects a control vector onto the box constraints.
  void projectBounds(Vector& x) const
  {
    for (double& v : x)
      v = std::clamp(v,lowerBound,upperBound);
  }

  //! \brief Prints an error message and returns the given error code.
  int error(const char* msg, int code) const
  {
    std::cerr <<" *** OptimizationDriver: "<< msg <<"."<< std::endl;
    return code;
  }

  //! \brief Assigns the optimization method from a string.
  void setMethodFromString(const std::string& str)
  {
    if (str == "newton-cg" || str == "newton_cg" || str == "ncg")
      method = Method::NEWTON_CG;
    else if (str == "lbfgs" || str == "l-bfgs" || str == "bfgs")
      method = Method::LBFGS;
    else if (str == "gradient_descent" || str == "steepest_descent" ||
             str == "gd")
      method = Method::GRADIENT_DESCENT;
    else
      std::cerr <<"  ** SIMSolverOpt: Unknown optimization method \""<< str
                <<"\" ignored."<< std::endl;
  }

  //! \brief Returns the name of the current optimization method.
  const char* methodName() const
  {
    switch (method) {
    case Method::NEWTON_CG:        return "Inexact Newton-CG";
    case Method::LBFGS:            return "L-BFGS";
    case Method::GRADIENT_DESCENT: return "Steepest descent";
    }
    return "Unknown";
  }

  //! \brief Prints one line of the convergence history.
  //! \param[in] it Iteration counter
  //! \param[in] obj Current objective value
  //! \param[in] gnorm Current gradient norm
  //! \param[in] step Accepted line search step size
  //! \param[in] inner Number of inner iterations, negative if not applicable
  void printIteration(int it, double obj, double gnorm,
                      double step, int inner) const
  {
    std::streamsize oldPrec = IFEM::cout.precision();
    std::ios::fmtflags oldFlags = IFEM::cout.flags();

    IFEM::cout << std::setw(6) << it << std::scientific
               <<"  "<< std::setprecision(8) << obj
               <<"  "<< std::setprecision(6) << gnorm;
    if (it > 0)
      IFEM::cout <<"  "<< std::setprecision(2) << step;
    if (inner >= 0 && method == Method::NEWTON_CG)
      IFEM::cout << std::setw(9) << inner;
    IFEM::cout << std::endl;

    IFEM::cout.flags(oldFlags);
    IFEM::cout.precision(oldPrec);
  }

  Method method = Method::NEWTON_CG; //!< Optimization method to use
  int    maxIter = 50;      //!< Maximum number of optimization iterations
  int    maxCGIter = 100;   //!< Maximum number of inner CG iterations
  int    lbfgsM = 10;       //!< Number of secant pairs stored by L-BFGS
  double tol = 1.0e-9;      //!< Absolute gradient norm tolerance
  double rtol = 1.0e-8;     //!< Relative gradient norm tolerance
  double stol = 1.0e-12;    //!< Smallest acceptable line search step size
  double armijoC1 = 1.0e-4; //!< Sufficient decrease parameter

  double gNormRef = 0.0;    //!< Reference gradient norm of the relative test
  int    nStagnant = 0;     //!< Number of iterations without gradient decrease

  bool   hasBounds = false;      //!< If \e true, the control is box-constrained
  double lowerBound = -1.0e30;   //!< Lower bound on the control values
  double upperBound =  1.0e30;   //!< Upper bound on the control values

  //! Maximum number of line search trials per optimization iteration
  static constexpr int maxLineSearch = 40;
  //! Smallest curvature accepted for an L-BFGS secant pair
  static constexpr double curvTol = 1.0e-14;
  //! Widest band around a bound within which a variable may be declared active
  static constexpr double bndTol = 1.0e-8;
  //! Smallest relative gradient decrease counted as progress
  static constexpr double stagTol = 1.0e-3;
  //! Number of consecutive stagnant iterations accepted before giving up
  static constexpr int maxStagnant = 3;
};


#endif
