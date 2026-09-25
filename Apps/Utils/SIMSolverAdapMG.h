// $Id$
//==============================================================================
//!
//! \file SIMSolverAdapMG.h
//!
//! \date Sep 23 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Stationary, adaptive SIM solver class template with geometric
//! multigrid built from the adaptive mesh sequence.
//!
//==============================================================================

#ifndef _SIM_SOLVER_ADAP_MG_H_
#define _SIM_SOLVER_ADAP_MG_H_

#include "MultigridProvider.h"
#include "SIMSolverAdap.h"

#include "ASMbase.h"
#include "IFEM.h"
#include "LogStream.h"
#include "MultigridTransfer.h"
#include "PETScMatrix.h"
#include "DomainDecomposition.h"
#include "ProcessAdm.h"
#include "SAM.h"
#include "SparseMatrix.h"
#include "TopologySet.h"
#include "LR/ASMLRSpline.h"
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
  \brief Adaptive simulator driver which installs a multigrid hierarchy.
  \details This is the AdaptiveSIM of SIMSolverAdap with one addition: before
  the equation system of a refinement step is solved, a hook is given the
  chance to install a multigrid hierarchy into it. The hook has to run at this
  point, since the system matrix does not exist until SIMbase::initSystem has
  been called from AdaptiveSIM::solveStep, and the preconditioner is built on
  the first solve after that.
*/

class AdaptiveMGSIM : public AdaptiveSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  AdaptiveMGSIM(SIMoutput& sim, bool sa = false) : AdaptiveSIM(sim,sa) {}

  //! \brief Sets the hook installing the multigrid hierarchy.
  void setMGInstaller(const std::function<bool()>& hook) { installMG = hook; }

protected:
  //! \brief Assembles and solves the linearized FE equation system.
  bool assembleAndSolveSystem() override
  {
    if (installMG && !installMG())
      return false;

    return this->AdaptiveSIM::assembleAndSolveSystem();
  }

private:
  std::function<bool()> installMG; //!< Installs the hierarchy, if given
};


/*!
  \brief Template class for adaptive simulator drivers using multigrid.

  \details The mesh sequence an adaptive simulation walks through is a
  multigrid hierarchy already: the spaces are nested, since refinement only
  inserts knot lines, and the meshes are graded towards whatever the error
  estimator points at. This driver keeps a chosen subset of those meshes
  alive, builds the transfer operators between them, and hands the result to
  PETSc as the levels of a PCMG preconditioner.

  Which meshes to keep is up to the input file:
  \code
    <multigrid levels="1 3"/>
    <multigrid stride="2"/>
  \endcode
  The listed levels are the \e coarse levels of the hierarchy; the mesh
  currently being solved on is always the finest, so \a levels="1 3" gives a
  three level hierarchy on the fifth refinement step. Leaving the tag out
  keeps every level, which is the finest hierarchy available but also the most
  expensive to store. Skipping levels is the interesting knob: it raises the
  coarsening ratio between neighbouring levels, which usually costs multigrid
  convergence rate and wants more smoothing steps in return.

  The coarse levels can also be meshes made outside the simulation, which is
  what a mesh generator writing one file per level gives:
  \code
    <geometry><patchfile>square-L3.g2</patchfile></geometry>
    <multigrid>
      <level>square-L1.g2</level>
      <level>square-L2.g2</level>
    </multigrid>
  \endcode
  These are geometry files, not input files: everything other than the mesh
  keeps coming from the one input file, which each level reads for itself
  after swapping its mesh for the one in the file. The geometry of the input
  file is the mesh being solved on, as it is in any other simulation, and the
  levels are listed coarsest first below it. They are built before the first
  solve and the simulator solving the problem never gives up its own mesh, so
  a hierarchy made this way costs the simulation nothing but the memory of
  the levels. Levels the adaptive loop keeps afterwards stack on top.

  The meshes have to be nested, each in the next, since the transfer operators
  are exact changes of basis and nothing else is meaningful. Generate them by
  inserting nested knot sets into one and the same geometry rather than by
  removing knots from the finest mesh, which perturbs the geometry unless the
  knots happen to be removable. MG::prolongation checks that the spaces really
  are nested and refuses to build an operator between meshes which are not.

  The operators a hierarchy is built for come from the simulator through the
  MultigridProvider interface, and are selected in the linear solver input
  with \a pc="gmg" on the matching block.
*/

template<class T1, template<class,class> class AdapImpl = SIMSolverAdapImpl>
class SIMSolverAdapMGImpl : public AdapImpl<T1,AdaptiveMGSIM>
{
  using Base = AdapImpl<T1,AdaptiveMGSIM>; //!< Base class alias

public:
  //! \brief The constructor forwards to the parent class constructor.
  //!
  //! \details \a AdapImpl is the adaptive driver to add the multigrid
  //! machinery to. It defaults to the plain SIMSolverAdapImpl, but an
  //! application with an adaptive driver of its own, such as the one Stokes
  //! uses to pick the norm to adapt on, passes that one instead.
  explicit SIMSolverAdapMGImpl(T1& s1) : Base(s1)
  {
    this->aSim.setMGInstaller([this]() { return this->installHierarchy(); });
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverAdapMGImpl() {}

  /*!
    \brief Solves the problem on a sequence of adaptively refined meshes.

    \details This is the loop of SIMSolverAdapImpl with the levels kept along
    the way. A level is kept \e after it has been solved on, since while it is
    being solved on it is the finest level of the hierarchy, not a coarse one.
  */
  int solveProblem(char* infile, const char* = nullptr) override
  {
    inputFile = infile;

    if (!this->aSim.initAdaptor())
      return 1;

    if (Base::exporter)
      Base::exporter->setFieldValue(this->exporterName, &this->S1,
                                    &this->aSim.getSolution(),
                                    &this->aSim.getProjections(),
                                    &this->aSim.getEnorm());

    if (!this->buildMeshLevels(infile))
      return 5;

    for (int iStep = 1; this->aSim.adaptMesh(iStep); iStep++)
      if (!this->aSim.solveStep(infile,iStep))
        return 1;
      else if (!this->aSim.writeGlv(infile,iStep))
        return 2;
      else if (Base::exporter && !Base::exporter->dumpTimeLevel(nullptr,true))
        return 3;
      else if (this->isCoarseLevel(iStep) && !this->keepLevel())
        return 4;

    return 0;
  }

protected:
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    if (strcasecmp(elem->Value(),"multigrid"))
      return this->Base::parse(elem);

    std::string levelList;
    if (utl::getAttribute(elem,"levels",levelList))
      utl::parseIntegers(coarseLevels,levelList.c_str());
    utl::getAttribute(elem,"stride",stride);
    utl::getAttribute(elem,"galerkin",galerkin);
    utl::getAttribute(elem,"lines",lineDir);
    utl::getAttribute(elem,"set",lineSet);

    for (const tinyxml2::XMLElement* child = elem->FirstChildElement();
         child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"level") && child->FirstChild())
        levelMeshes.push_back(child->FirstChild()->Value());

    IFEM::cout <<"\tGeometric multigrid: ";
    if (!levelMeshes.empty()) {
      IFEM::cout << levelMeshes.size() <<" mesh level(s)";
      for (const std::string& f : levelMeshes) IFEM::cout <<" "<< f;
      IFEM::cout <<", then ";
    }
    if (coarseLevels.empty())
      IFEM::cout <<"every "<< (stride > 1 ? std::to_string(stride)+". " : "")
                 <<"level";
    else {
      IFEM::cout <<"levels";
      for (int l : coarseLevels) IFEM::cout <<" "<< l;
    }
    if (galerkin)
      IFEM::cout <<", Galerkin coarse operators";
    IFEM::cout << std::endl;

    return true;
  }

  //! \brief Returns \e true if a refinement step is kept as a coarse level.
  bool isCoarseLevel(int iStep) const
  {
    if (coarseLevels.empty())
      return (iStep-1)%stride == 0;

    return std::find(coarseLevels.begin(),coarseLevels.end(),iStep)
           != coarseLevels.end();
  }

  /*!
    \brief Keeps the mesh of the current refinement step as a coarse level.

    \details The adaptive loop refines its patches in place, so a level has to
    be copied out before the next refinement destroys it. The copy is a
    simulator of the same type, reading the same input file so that it gets
    the same properties, boundary conditions and integrand, but given a
    private deep copy of the mesh as it stands instead of the one the input
    file describes. This is the same sequence AdaptiveSIM::solveStep goes
    through after a refinement.
  */
  bool keepLevel()
  {
    std::unique_ptr<MultigridProvider> level;
    T1* sim = this->makeLevel(inputFile,level);
    if (!sim)
      return false;

    const ASM::PatchVec& myModel = this->S1.getFEModel();
    if (sim->getFEModel().size() != myModel.size()) {
      std::cerr <<" *** SIMSolverAdapMG: The level has "
                << sim->getFEModel().size() <<" patches, the model has "
                << myModel.size() <<"."<< std::endl;
      return false;
    }

    for (size_t i = 0; i < myModel.size(); i++)
      if (!sim->getFEModel()[i]->copyMeshFrom(*myModel[i]))
        return false;

    // The topology was resolved while the input file was read, against the
    // mesh that file names, and the copy above has replaced that mesh. It
    // has to be resolved once more for the mesh now in place, or the patches
    // are left unjoined. Reading the file again does that and keeps the
    // model which is there, which is how a refinement is followed up.
    sim->clearProperties();
    if (!sim->read(inputFile))
      return false;

    if (!this->addLevel(sim,level,"Keeping"))
      return false;

    // The level holds a copy of the mesh being solved on, so it has to come
    // out the same size. It does not when the patches of a multi-patch model
    // fail to be joined along their interfaces after the copy, and a level
    // whose patches are loose is not the space it claims to be.
    const int nEq = sim->getSAM()->getNoEquations();
    const int nRef = this->S1.getSAM()->getNoEquations();
    if (nEq != nRef) {
      std::cerr <<" *** SIMSolverAdapMG: The level kept has "<< nEq
                <<" equations where the mesh it copied has "<< nRef
                <<". Its patches were not joined the same way."<< std::endl;
      return false;
    }

    return true;
  }

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
        std::cerr <<" *** SIMSolverAdapMG: The simulator did not create a"
                  <<" level of its own type."<< std::endl;
        return false;
      }

      sim->opt = this->S1.opt;
      if (!sim->readMesh(mesh) || !sim->read(infile))
        return false;

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

      const ASMLRSpline* lrPch = dynamic_cast<const ASMLRSpline*>(pch);
      std::vector<IntVec> lines;
      if (!lrPch || !lrPch->getLineDofs(lineDir,lines))
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

  //! \brief Creates a level simulator which has read the given input file.
  //! \param[in] file The input file the level reads
  //! \param[out] level The level, owning the simulator returned
  T1* makeLevel(const char* file, std::unique_ptr<MultigridProvider>& level)
  {
    level = this->S1.createMGLevel();
    T1* sim = dynamic_cast<T1*>(level.get());
    if (!sim) {
      std::cerr <<" *** SIMSolverAdapMG: The simulator did not create a level"
                <<" of its own type."<< std::endl;
      return nullptr;
    }

    sim->opt = this->S1.opt;
    if (!sim->read(file))
      return nullptr;

    return sim;
  }

  //! \brief Preprocesses a level, assembles its operators and keeps it.
  //! \param[in] sim The level simulator
  //! \param level The level, taken over by this driver
  //! \param[in] what What to call the level in the log
  bool addLevel(T1* sim, std::unique_ptr<MultigridProvider>& level,
                const char* what)
  {
    if (!sim->preprocess())
      return false;

    IFEM::cout <<"\t"<< what <<" level "<< 1+levels.size() <<" with "
               << sim->getSAM()->getNoEquations() <<" equations for multigrid"
               << std::endl;

    // Assemble the operators on this level while the mesh is current
    if (!galerkin)
      for (const MG::Operator& op : this->S1.getMGOperators()) {
        SystemMatrix* A = level->assembleMGOperator(op);
        if (!A) {
          std::cerr <<" *** SIMSolverAdapMG: Could not assemble the operator"
                    <<" \""<< op.name <<"\" on level "<< 1+levels.size()
                    <<"."<< std::endl;
          return false;
        }
        levelOps[op.name].push_back(A);
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
      std::vector<std::unique_ptr<SparseMatrix>>& P = prolong[op.name];
      std::vector<std::pair<int,int>>& L = layout[op.name];

      // The operators between the kept levels only have to be built once,
      // and so is the share of each one this process owns.
      for (size_t i = P.size(); i+1 < sims.size(); i++)
      {
        int nRow = 0, nCol = 0;
        if (!(P.emplace_back(MG::prolongation(*sims[i],*sims[i+1],op,
                                              MG::Transfer::CHANGE_OF_BASIS,
                                              &nRow,&nCol))).get())
          return false;
        L.emplace_back(nRow,nCol);
      }
      P.resize(sims.size()-1);
      L.resize(sims.size()-1);

      // The one onto the mesh being solved on changes with every refinement
      int nRow = 0, nCol = 0;
      std::unique_ptr<SparseMatrix> top =
        MG::prolongation(*sims.back(),this->S1,op,
                         MG::Transfer::CHANGE_OF_BASIS,&nRow,&nCol);
      if (!top)
        return false;

      std::vector<std::pair<int,int>> owned(L);
      owned.emplace_back(nRow,nCol);

      std::vector<const SparseMatrix*> Pptr;
      for (const std::unique_ptr<SparseMatrix>& p : P)
        Pptr.push_back(p.get());
      Pptr.push_back(top.get());

      std::vector<const SystemMatrix*> Aptr;
      if (!galerkin)
        for (const SystemMatrix* A : levelOps[op.name])
          Aptr.push_back(A);

      std::vector<std::vector<IntVec>> subd;
      if (lineDir > 0) {
        for (const T1* sim : sims)
          subd.push_back(this->lineSubdomains(*sim,op.block));
        subd.push_back(this->lineSubdomains(this->S1,op.block));
      }

      // setMGHierarchy converts the transfer operators to PETSc format, so
      // the topmost one is not needed beyond this point. The level operators
      // are not copied, and stay owned by the level simulators.
      if (!pA->setMGHierarchy(op.block,Pptr,Aptr,subd,owned))
        return false;
    }

    return true;
  }

  char* inputFile = nullptr; //!< Input file, re-read when keeping a level

  //! Meshes used as coarse levels, coarsest first, read from geometry files
  std::vector<std::string> levelMeshes;

  IntVec coarseLevels; //!< Refinement steps kept as coarse levels
  int    stride = 1;   //!< Keep every this many levels, if none are listed
  bool   galerkin = false; //!< Let PETSc form the coarse operators
  int    lineDir = 0;      //!< Direction mesh lines run along, 0 for no lines
  //! Topology set of the patches the lines are taken from, empty for all
  std::string lineSet;

  std::vector<std::unique_ptr<MultigridProvider>> levels; //!< The kept levels
  std::vector<T1*> sims; //!< The kept levels, as simulators

  //! Operators assembled on each kept level, by operator name.
  //! They are owned by the level simulators, which \ref levels keeps alive.
  std::map<std::string,std::vector<SystemMatrix*>> levelOps;
  //! Transfer operators between the kept levels, by operator name
  std::map<std::string,std::vector<std::unique_ptr<SparseMatrix>>> prolong;
  //! Rows and columns of each kept transfer operator this process owns
  std::map<std::string,std::vector<std::pair<int,int>>> layout;
};


//! Convenience alias template
template<class T1>
using SIMSolverAdapMG = SIMSolverAdapMGImpl<T1,SIMSolverAdapImpl>;

#endif
