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
//! multigrid.
//!
//==============================================================================

#ifndef _SIM_SOLVER_ADAP_MG_H_
#define _SIM_SOLVER_ADAP_MG_H_

#include "SIMSolverAdap.h"
#include "SIMSolverMG.h"


/*!
  \brief Template class for adaptive simulator drivers using multigrid.

  \details The mesh sequence an adaptive simulation walks through is a
  multigrid hierarchy already: the spaces are nested, since refinement only
  inserts knot lines, and the meshes are graded towards whatever the error
  estimator points at. This driver is SIMSolverMGImpl with those meshes added
  to the hierarchy as they are walked through, a chosen subset of them kept
  alive rather than refined away.

  Which meshes to keep is up to the input file:
  \code
    <multigrid levels="1 3"/>
    <multigrid stride="2"/>
    <multigrid levels="fixed"/>
  \endcode
  The listed levels are the \e coarse levels of the hierarchy; the mesh
  currently being solved on is always the finest, so \a levels="1 3" gives a
  three level hierarchy on the fifth refinement step. Leaving the tag out
  keeps every level, which is the finest hierarchy available but also the most
  expensive to store. Skipping levels is the interesting knob: it raises the
  coarsening ratio between neighbouring levels, which usually costs multigrid
  convergence rate and wants more smoothing steps in return.

  \a levels="fixed" keeps none of them, which leaves the hierarchy as the
  levels read from geometry files built it, however far the mesh being solved
  on is refined past the finest of those.

  Levels the adaptive loop keeps stack on top of the ones read from files,
  which are described in SIMSolverMGImpl along with everything else the
  \a multigrid tag takes.
*/

template<class T1, template<class,class> class AdapImpl = SIMSolverAdapImpl>
class SIMSolverAdapMGImpl : public SIMSolverMGImpl<T1,AdapImpl<T1,AdaptiveSIM>>
{
  //! Base class alias
  using Base = SIMSolverMGImpl<T1,AdapImpl<T1,AdaptiveSIM>>;

public:
  //! \brief The constructor forwards to the parent class constructor.
  //!
  //! \details \a AdapImpl is the adaptive driver the multigrid machinery is
  //! added to. It defaults to the plain SIMSolverAdapImpl, but an application
  //! with an adaptive driver of its own, such as the one Stokes uses to pick
  //! the norm to adapt on, passes that one instead.
  explicit SIMSolverAdapMGImpl(T1& s1) : Base(s1) {}

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
  //! \brief Parses the refinement steps whose meshes are kept as levels.
  bool parseMore(const tinyxml2::XMLElement* elem) override
  {
    std::string levelList;
    if (utl::getAttribute(elem,"levels",levelList))
    {
      if (levelList == "fixed")
        fixedLevels = true;
      else
        utl::parseIntegers(coarseLevels,levelList.c_str());
    }
    utl::getAttribute(elem,"stride",stride);

    return true;
  }

  //! \brief Returns which of the refinement steps the hierarchy keeps.
  std::string moreLevels() const override
  {
    if (fixedLevels)
      return ""; // the hierarchy is what it was built as

    if (coarseLevels.empty())
      return "every " + (stride > 1 ? std::to_string(stride)+". "
                                    : std::string()) + "level";

    std::string what("levels");
    for (int l : coarseLevels)
      what += " " + std::to_string(l);

    return what;
  }

  //! \brief Returns that the mesh being solved on is refined as it goes.
  bool refinesMesh() const override { return true; }

  //! \brief Returns \e true if a refinement step is kept as a coarse level.
  bool isCoarseLevel(int iStep) const
  {
    if (fixedLevels)
      return false; // the hierarchy is what it was built as

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
    {
      typename Base::Quiet quiet;
      sim->clearProperties();
      if (!sim->read(inputFile))
        return false;
    }

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

    // A level kept mid-solve is followed by the norms of the step it was
    // kept after, which it would otherwise run straight into.
    IFEM::cout << std::endl;

    return true;
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

    typename Base::Quiet quiet;
    if (!sim->read(file))
      return nullptr;

    return sim;
  }

  char* inputFile = nullptr; //!< Input file, re-read when keeping a level

  IntVec coarseLevels; //!< Refinement steps kept as coarse levels
  int    stride = 1;   //!< Keep every this many levels, if none are listed
  //! Whether the hierarchy stays as it was built, keeping no level of the
  //! refinements. A hierarchy read from geometry files is then the whole of
  //! it, however far the mesh being solved on is refined past the finest of
  //! them.
  bool   fixedLevels = false;
};


//! Convenience alias template
template<class T1>
using SIMSolverAdapMG = SIMSolverAdapMGImpl<T1,SIMSolverAdapImpl>;

#endif
