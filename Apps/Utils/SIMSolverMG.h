// $Id$
//==============================================================================
//!
//! \file SIMSolverMG.h
//!
//! \date Sep 26 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Stationary SIM solver class template with geometric multigrid.
//!
//==============================================================================

#ifndef _SIM_SOLVER_MG_H_
#define _SIM_SOLVER_MG_H_

#include "SIMSolver.h"

#include "MultigridProvider.h"
#include "MultigridTransfer.h"

#include "ASMbase.h"
#include "IFEM.h"
#include "LogStream.h"
#include "PETScMatrix.h"
#include "DomainDecomposition.h"
#include "ProcessAdm.h"
#include "SAM.h"
#include "SparseMatrix.h"
#include "TopologySet.h"
#include "Utilities.h"
#include "tinyxml2.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <vector>


/*!
  \brief Simulator driver which solves on the finest of a mesh hierarchy.

  \details This is the machinery on its own, added to whichever driver runs
  the simulation: SIMSolverStatMG solves once on the finest mesh, SIMSolverMG
  steps it through time, and SIMSolverAdapMG refines it and adds the meshes it
  walks through to the hierarchy.

  The levels are meshes of one and the same geometry, each nested in
  the next, which the driver turns into the levels of a PCMG preconditioner:
  it reads them, assembles the operators on each of them, builds the transfer
  operators between them and hands the result to PETSc. The mesh the input
  file names is the finest level and the one the problem is solved on, as it
  is in any other simulation.

  The coarse levels are geometry files, which is what a mesh generator writing
  one file per level gives:
  \code
    <geometry><patchfile>square-L3.g2</patchfile></geometry>
    <multigrid lines="2" set="boundary_layer">
      <level>square-L1.g2</level>
      <level>square-L2.g2</level>
    </multigrid>
  \endcode
  They are geometry files, not input files: everything other than the mesh
  keeps coming from the one input file, which each level reads for itself
  after swapping its mesh for the one in the file. The levels are listed
  coarsest first, and the simulator solving the problem never gives up its own
  mesh, so a hierarchy costs the simulation nothing but the memory of the
  levels.

  The levels may be tensor product splines or locally refined ones, the same
  kind on every level, and the transfer operators are built either way.

  The meshes have to be nested, each in the next, since the transfer operators
  are exact changes of basis and nothing else is meaningful. Generate them by
  inserting nested knot sets into one and the same geometry rather than by
  removing knots from the finest mesh, which perturbs the geometry unless the
  knots happen to be removable. MG::prolongation checks that the spaces really
  are nested and refuses to build an operator between meshes which are not.

  A level may lower the polynomial order as well as coarsen the mesh, which
  is not a nesting and has no change of basis; the coarse basis is projected
  onto the fine one instead, and the transfer applies the projection rather
  than storing it.

  The \a lines attribute turns the point smoother of the levels into one
  solving along the mesh lines of a parameter direction, which is what an
  anisotropic mesh needs, and \a set names the topology set of the patches to
  take those lines from, the others being left to the point smoother.

  The operators a hierarchy is built for come from the simulator through the
  MultigridProvider interface, and are selected in the linear solver input
  with \a pc="gmg" on the matching block.
*/

template<class T1, class Base>
class SIMSolverMGImpl : public Base
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //!
  //! \details \a Base is the driver the multigrid machinery is added to. It
  //! is the plain stationary driver here, and an adaptive one in
  //! SIMSolverAdapMGImpl, which also adds to the hierarchy as it refines.
  explicit SIMSolverMGImpl(T1& s1) : Base(s1)
  {
    this->S1.setPreSolveHook([this]() { return this->installHierarchy(); });
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverMGImpl() {}

  //! \brief Reads solver data from the specified input file.
  //! \details The hierarchy is described in the input file, which a driver
  //! with nothing of its own to read there does not open.
  bool read(const char* file) override { return this->SIMadmin::read(file); }

  //! \brief Solves the problem on the finest level of the hierarchy.
  //!
  //! \details The hierarchy is installed here as well as from the hook,
  //! since a simulator may have allocated its equation system while it was
  //! being configured, which is before this driver existed and before there
  //! was a hierarchy to install. Installing twice costs nothing: everything
  //! the second one hands over has been built by the first.
  int solveProblem(char* infile, const char* heading = nullptr) override
  {
    if (!this->buildMeshLevels(infile) || !this->installHierarchy())
      return 5;

    return this->Base::solveProblem(infile,heading);
  }

protected:
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    if (strcasecmp(elem->Value(),"multigrid"))
      return this->Base::parse(elem);

    utl::getAttribute(elem,"galerkin",galerkin);
    utl::getAttribute(elem,"lines",lineDir);
    utl::getAttribute(elem,"set",lineSet);

    for (const tinyxml2::XMLElement* child = elem->FirstChildElement();
         child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"level") && child->FirstChild())
        levelMeshes.push_back(child->FirstChild()->Value());

    if (!this->parseMore(elem))
      return false;

    IFEM::cout <<"\tGeometric multigrid: ";
    const std::string more = this->moreLevels();
    if (!levelMeshes.empty()) {
      IFEM::cout << levelMeshes.size() <<" mesh level(s)";
      for (const std::string& f : levelMeshes) IFEM::cout <<" "<< f;
      if (!more.empty())
        IFEM::cout <<", then ";
    }
    if (more.empty())
      IFEM::cout <<(levelMeshes.empty() ? "no levels" : " and no more");
    else
      IFEM::cout << more;
    if (galerkin)
      IFEM::cout <<", Galerkin coarse operators";
    IFEM::cout << std::endl;

    return true;
  }

  //! \brief Parses the settings of a driver which adds to the hierarchy.
  //! \details A driver which does not refine adds nothing, and has nothing
  //! to read beyond the levels read here.
  virtual bool parseMore(const tinyxml2::XMLElement*) { return true; }

  //! \brief Returns what the hierarchy gains past the levels it is built as.
  //! \details The empty string, which is what a driver that does not refine
  //! leaves it at, makes the log say that the hierarchy is the whole of it.
  virtual std::string moreLevels() const { return ""; }

  //! \brief Returns whether the mesh being solved on changes as it is.
  //! \details It does not unless something refines it, and what is collected
  //! from a mesh which stays as it is need only be collected once.
  virtual bool refinesMesh() const { return false; }

  /*!
    \brief Builds the levels whose meshes are given as geometry files.

    \details A level swaps in its own mesh before it reads the input file, so
    that the topology and the boundary conditions are established for that
    mesh and the geometry the input file names is skipped, a model which is
    already there being kept. This is the sequence an adaptive simulation goes
    through after a refinement, with the refinement replaced by reading a
    mesh, and it leaves the simulator solving the problem untouched.
  */
  bool buildMeshLevels(char* infile)
  {
    for (const std::string& mesh : levelMeshes) {
      std::unique_ptr<MultigridProvider> level = this->S1.createMGLevel();
      T1* sim = dynamic_cast<T1*>(level.get());
      if (!sim) {
        std::cerr <<" *** SIMSolverMG: The simulator did not create a"
                  <<" level of its own type."<< std::endl;
        return false;
      }

      {
        Quiet quiet;
        sim->opt = this->S1.opt;
        if (!sim->readMesh(mesh) || !sim->read(infile))
          return false;
      }

      if (!this->addLevel(sim,level,"Reading"))
        return false;
    }

    return true;
  }

  /*!
    \brief Collects the mesh lines of a level as sets of equations.

    \details A smoother solving these rather than the points is what an
    anisotropic mesh needs, and the lines have to be taken from the mesh of
    the level they smooth, not from the one being solved on. The lines of a
    patch are the functions whose knot vectors agree in every direction but
    the one the lines run along, which is a mesh line while the mesh is a
    tensor mesh, as the mesh of a boundary layer is.

    \param[in] sim The level the lines are taken from
    \param[in] block Matrix block the equations are numbered within
  */
  std::vector<IntVec> lineSubdomains(const T1& sim, size_t block) const
  {
    std::vector<IntVec> subd;
    const SAM* sam = sim.getSAM();
    if (!sam)
      return subd;

    // A line smoother is only worth its cost where the mesh is anisotropic,
    // which in a boundary layer mesh is the patches next to the wall. The
    // others are left to the point smoother, and saying which is which is a
    // job the topology sets already do.
    std::set<size_t> selected;
    if (!lineSet.empty()) {
      for (const TopItem& item : sim.getEntity(lineSet))
        if (abs(item.idim) == static_cast<int>(sim.getNoParamDim()))
          selected.insert(item.patch);

      if (selected.empty())
        IFEM::cout <<"  ** No patches in the topology set \""<< lineSet
                   <<"\", taking mesh lines from all of them."<< std::endl;
    }

    IntSet covered;
    for (const ASMbase* pch : sim.getFEModel())
    {
      if (!selected.empty() && selected.find(pch->idx+1) == selected.end())
        continue;

      std::vector<IntVec> lines;
      if (!pch->getLineDofs(lineDir,lines))
        continue;

      for (const IntVec& line : lines) {
        IntSet eqs;
        for (int inod : line) {
          IntVec meqn;
          sam->getNodeEqns(meqn,pch->getNodeID(inod));
          for (int eq : meqn)
            if (eq > 0)
              eqs.insert(eq-1); // The lines are zero-based within the level
        }
        if (!eqs.empty()) {
          subd.emplace_back(eqs.begin(),eqs.end());
          covered.insert(eqs.begin(),eqs.end());
        }
      }
    }

    if (subd.empty())
      return subd;

    // Patches which meet within the selection contribute a line each to the
    // equations they share, be it the two halves of a line cut by an
    // interface or the one line two patches both see along a seam. Joining
    // those into a single line makes the lines follow the mesh rather than
    // the patches, and leaves the smoother without the overlap which would
    // otherwise cost it its convergence.
    std::map<int,size_t> owner;
    std::vector<size_t> merge(subd.size());
    std::iota(merge.begin(),merge.end(),0);
    std::function<size_t(size_t)> root = [&merge,&root](size_t i)
    { return merge[i] == i ? i : merge[i] = root(merge[i]); };

    for (size_t i = 0; i < subd.size(); i++)
      for (int eq : subd[i]) {
        std::map<int,size_t>::iterator it = owner.find(eq);
        if (it == owner.end())
          owner[eq] = root(i);
        else
          merge[root(i)] = root(it->second);
      }

    std::map<size_t,IntSet> joined;
    for (size_t i = 0; i < subd.size(); i++)
      joined[root(i)].insert(subd[i].begin(),subd[i].end());

    subd.clear();
    subd.reserve(joined.size());
    for (const std::pair<const size_t,IntSet>& j : joined)
      subd.emplace_back(j.second.begin(),j.second.end());

    // A line of a partitioned mesh runs through the equations of several
    // processes, and each of them keeps the stretch it owns. A line is thus
    // cut where the partition cuts it, which costs the smoother some of its
    // strength but leaves every subdomain solvable without talking to
    // anyone. How much it costs is the number of lines the partition cut,
    // which is reported so that it can be weighed.
    // The equations are those of this level, which has a numbering of its
    // own, and the smoother wants the global numbers of the system it
    // preconditions. Only this level knows how to make that translation, so
    // it is made here rather than where the smoother is set up.
    const DomainDecomposition& dd = sim.getProcessAdm().dd;
    auto&& globalEq = [&dd,block](int eq)
    {
      const int geq = dd.getGlobalEq(eq+1,block);
      return geq >= dd.getMinEq(block) && geq <= dd.getMaxEq(block) ? geq : 0;
    };

    size_t cut = 0;
    std::vector<IntVec> mine;
    mine.reserve(subd.size());
    for (const IntVec& line : subd)
    {
      IntVec keep;
      keep.reserve(line.size());
      for (int eq : line)
        if (int geq = globalEq(eq); geq > 0)
          keep.push_back(geq-1); // PETSc counts equations from zero
      if (keep.empty())
        continue; // the line belongs to another process entirely

      if (keep.size() < line.size())
        ++cut; // the partitioning runs through this one

      std::sort(keep.begin(),keep.end());
      mine.push_back(std::move(keep));
    }
    subd.swap(mine);

    // The subdomains have to span the whole system, or the smoother they
    // define is singular. Whatever the lines did not reach is smoothed one
    // equation at a time, which is what a point smoother does anyway.
    int nOwned = 0;
    for (int eq = 0; eq < sam->getNoEquations(); eq++)
      if (int geq = globalEq(eq); geq > 0)
      {
        ++nOwned;
        if (covered.find(eq) == covered.end())
          subd.emplace_back(1,geq-1);
      }

    if (subd.empty())
      return subd;

    // The counts are of this process, whose stream they are written to.
    size_t mn = subd.front().size(), mx = mn, tot = 0;
    for (const IntVec& d : subd) {
      mn = std::min(mn,d.size());
      mx = std::max(mx,d.size());
      tot += d.size();
    }
    IFEM::cout <<"\tMesh lines: "<< subd.size() <<" subdomains of "<< mn
               <<" to "<< mx <<" equations, covering "<< tot <<" of "
               << nOwned;
    if (cut > 0)
      IFEM::cout <<", "<< cut <<" cut by the partitioning";
    IFEM::cout << std::endl;

    return subd;
  }


  /*!
    \brief Silences the log for as long as it is alive.

    \details A level reads the input file and preprocesses a model of its
    own, which says everything about itself that the simulation being run
    says, and none of it was asked for. What the levels have to say about
    themselves is said by addLevel() in one line each.
  */
  class Quiet
  {
  public:
    //! \brief The constructor silences the log.
    Quiet() : wasMuted(IFEM::cout.mute(true)) {}
    //! \brief The destructor restores it.
    ~Quiet() { IFEM::cout.mute(wasMuted); }

  private:
    bool wasMuted; //!< The setting to restore
  };

  //! \brief Preprocesses a level, assembles its operators and keeps it.
  //! \param[in] sim The level simulator
  //! \param level The level, taken over by this driver
  //! \param[in] what What to call the level in the log
  bool addLevel(T1* sim, std::unique_ptr<MultigridProvider>& level,
                const char* what)
  {
    {
      Quiet quiet;
      if (!sim->preprocess())
        return false;
    }

    IFEM::cout <<"\t"<< what <<" level "<< 1+levels.size() <<" with "
               << sim->getSAM()->getNoEquations() <<" equations for multigrid"
               << std::endl;

    // Assemble the operators on this level while the mesh is current. That
    // has as much to say for itself as an assembly of the system being
    // solved, and none of it was asked for either.
    if (!galerkin)
      for (const MG::Operator& op : this->S1.getMGOperators()) {
        Quiet quiet;
        std::unique_ptr<SystemMatrix> A = level->assembleMGOperator(op);
        if (!A) {
          std::cerr <<" *** SIMSolverMG: Could not assemble the operator"
                    <<" \""<< op.name <<"\" on level "<< 1+levels.size()
                    <<"."<< std::endl;
          return false;
        }
        levelOps[op.name].push_back(std::move(A));
      }

    levels.push_back(std::move(level));
    sims.push_back(sim);

    return true;
  }

  /*!
    \brief Builds the transfer operators and installs them in the solver.

    \details The finest level is always the mesh currently being solved on,
    so the topmost transfer operator is rebuilt on every step while the ones
    between the kept levels stay put once they have been computed.
  */
  bool installHierarchy()
  {
    if (levels.empty())
      return true; // nothing kept yet, solve the first mesh as usual

    PETScMatrix* pA = dynamic_cast<PETScMatrix*>(this->S1.getLHSmatrix());
    if (!pA)
      return true; // not solving with PETSc, nothing to install

    for (const MG::Operator& op : this->S1.getMGOperators()) {
      std::vector<std::unique_ptr<MG::Prolongation>>& P = prolong[op.name];

      // The operators between the kept levels only have to be built once
      for (size_t i = P.size(); i+1 < sims.size(); i++)
        if (!(P.emplace_back(MG::prolongation(*sims[i],*sims[i+1],op))).get())
          return false;
      P.resize(sims.size()-1);

      // The one onto the mesh being solved on changes whenever that mesh is
      // refined, and a driver which does not refine builds it once like the
      // others.
      std::unique_ptr<MG::Prolongation>& top = topProlong[op.name];
      if (!top || this->refinesMesh())
        top = MG::prolongation(*sims.back(),this->S1,op);
      if (!top)
        return false;

      // A level which lowers the order has no operator to hand over, only
      // the two factors of the projection which stands in for one.
      std::vector<const SparseMatrix*> Pptr, massPtr;
      std::vector<std::pair<int,int>> owned;
      std::vector<bool> spread;
      for (const std::unique_ptr<MG::Prolongation>& p : P)
      {
        Pptr.push_back(p->layout());
        massPtr.push_back(p->mass.get());
        owned.emplace_back(p->rowsOwned,p->colsOwned);
        spread.push_back(p->distributed);
      }
      Pptr.push_back(top->layout());
      massPtr.push_back(top->mass.get());
      owned.emplace_back(top->rowsOwned,top->colsOwned);
      spread.push_back(top->distributed);

      std::vector<const SystemMatrix*> Aptr;
      if (!galerkin)
        for (const std::unique_ptr<SystemMatrix>& A : levelOps[op.name])
          Aptr.push_back(A.get());

      // Collecting the lines of a mesh means walking its functions, joining
      // what the patch interfaces cut and cutting again where the partition
      // does, which is worth doing once per mesh rather than once per solve.
      // A level never changes after it has been added, and the mesh being
      // solved on only changes under a driver which refines it.
      std::vector<std::vector<IntVec>> subd;
      if (lineDir > 0) {
        std::vector<std::vector<IntVec>>& kept = levelLines[op.name];
        for (size_t i = kept.size(); i < sims.size(); i++)
          kept.push_back(this->lineSubdomains(*sims[i],op.block));

        std::vector<IntVec>& top = fineLines[op.name];
        if (top.empty() || this->refinesMesh())
          top = this->lineSubdomains(this->S1,op.block);

        subd = kept;
        subd.push_back(top);
      }

      // setMGHierarchy converts the transfer operators to PETSc format, so
      // the topmost one is not needed beyond this point. The level operators
      // are not copied, and stay owned by the level simulators.
      if (!pA->setMGHierarchy(op.block,Pptr,Aptr,subd,owned,massPtr,spread))
        return false;
    }

    return true;
  }

  //! Meshes used as coarse levels, coarsest first, read from geometry files
  std::vector<std::string> levelMeshes;

  bool galerkin = false; //!< Let PETSc form the coarse operators
  int  lineDir = 0;      //!< Direction mesh lines run along, 0 for no lines
  //! Topology set of the patches the lines are taken from, empty for all
  std::string lineSet;

  /*! The kept levels. A level is asked for its mesh only while the transfer
      operator onto the level above it is built and while its mesh lines are
      collected, both of which happen once, so there is a good deal of it
      which could be let go of afterwards: on a locally refined level the
      mesh is of the order of the operator itself.

      It is kept all the same, because an operator handed over by a level
      still points back at it. A PETScMatrix holds a reference to the
      ProcessAdm of the simulator which made it, and that destructor frees
      the communicator the operator was created on, which PETSc goes through
      whenever the preconditioner is set up. Letting go of a level therefore
      costs an invalid communicator on the next cycle. Making a matrix own
      its process administrator and solver parameters rather than refer to
      them is what this waits on. */
  std::vector<std::unique_ptr<MultigridProvider>> levels;
  std::vector<T1*> sims; //!< The kept levels, as simulators

  //! Operators assembled on each kept level, by operator name. They are
  //! handed over by the level simulators and outlive them.
  std::map<std::string,std::vector<std::unique_ptr<SystemMatrix>>> levelOps;
  //! Transfer operators between the kept levels, by operator name
  std::map<std::string,std::vector<std::unique_ptr<MG::Prolongation>>> prolong;
  //! Transfer operator onto the mesh being solved on, by operator name
  std::map<std::string,std::unique_ptr<MG::Prolongation>> topProlong;

  //! Mesh lines of each kept level, by operator name
  std::map<std::string,std::vector<std::vector<IntVec>>> levelLines;
  //! Mesh lines of the mesh being solved on, by operator name
  std::map<std::string,std::vector<IntVec>> fineLines;
};


//! \brief Stationary simulator driver using geometric multigrid.
template<class T1>
using SIMSolverStatMG = SIMSolverMGImpl<T1,SIMSolverStat<T1>>;

//! \brief Time stepping simulator driver using geometric multigrid.
//! \details The hierarchy is installed into each equation system the
//! simulator allocates, which is once for a mesh it keeps through the
//! simulation and once per step for one it does not.
template<class T1>
using SIMSolverMG = SIMSolverMGImpl<T1,SIMSolver<T1>>;

#endif
