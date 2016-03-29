// $Id$
//==============================================================================
//!
//! \file DomainDecomposition.h
//!
//! \date Feb 23 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Domain decomposition related partitioning for structured models.
//!
//==============================================================================

#ifndef _DOMAIN_DECOMPOSITION_H
#define _DOMAIN_DECOMPOSITION_H

#include <map>
#include <vector>
#include <cstddef>


class LinSolParams;
class ProcessAdm;
class SAMpatch;
class SIMbase;


/*!
  \brief Class containing domain decomposition related partitioning.
*/

class DomainDecomposition
{
public:
  //! \brief Struct defining a domain interface.
  struct Interface {
    int master; //!< Master patch (global number).
    int slave;  //!< Slave patch (global number).
    int midx;   //!< Index of boundary on master:
    int sidx;   //!< Index of boundary on slave.
    int orient; //!< Orientation.
    int dim;    //!< Dimension of boundary.
  };

  std::vector<Interface> ghostConnections; //!< Connections to other processes.

  //! \brief Returns number of spatial dimensions
  int getNoSpaceDim() const { return nsd; }

  //! \brief Setup domain decomposition
  bool setup(const ProcessAdm& adm, const SIMbase& sim, const LinSolParams* spar);

  //! \brief Get first equation owned by this process.
  int getMinEq() const { return minEq; };

  //! \brief Get last equation owned by this process.
  int getMaxEq() const { return maxEq; };

  //! \brief Get first node owned by this process.
  int getMinNode() const { return minNode; }

  //! \brief Get last node owned by this process.
  int getMaxNode() const { return maxNode; }

  //! \brief Get first DOF owned by this process.
  int getMinDOF() const { return minDof; }

  //! \brief Get last DOF owned by this process.
  int getMaxDOF() const { return maxDof; }

  //! \brief Calculates a partitioning with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] nel2 Number of knot-spans in second parameter direction.
  //! \param[in] nel3 Number of knot-spans in third parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] g2 Number of subdomains in second parameter direction.
  //! \param[in] g3 Number of subdomains in third parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  //! \details nel values determine the dimensionality.
  void calcAppropriateGroups(size_t nel1, size_t nel2, size_t nel3,
                             size_t g1, size_t g2, size_t g3, size_t overlap);

  //! \brief Calculates a 1D partitioning with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  void calcGroups(size_t nel1, size_t g1, size_t overlap);

  //! \brief Calculates a 2D partitioning with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] nel2 Number of knot-spans in second parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] g2 Number of subdomains in second parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  void calcGroups(size_t nel1, size_t nel2,
                  size_t g1, size_t g2, size_t overlap);

  //! \brief Calculates a 3D partitioning with a given overlap.
  //! \param[in] nel1 Number of knot-spans in first parameter direction.
  //! \param[in] nel2 Number of knot-spans in second parameter direction.
  //! \param[in] nel2 Number of knot-spans in third parameter direction.
  //! \param[in] g1 Number of subdomains in first parameter direction.
  //! \param[in] g2 Number of subdomains in second parameter direction.
  //! \param[in] g3 Number of subdomains in third parameter direction.
  //! \param[in] overlap Overlap of subdomains.
  void calcGroups(size_t nel1, size_t nel2, size_t nel3,
                  size_t g1, size_t g2, size_t g3, size_t overlap);

  //! \brief Indexing operator.
  const std::vector<size_t>& operator[](int i) const { return subdomains[i]; }

  //! \brief Returns number of subdomains.
  size_t getNoSubdomains() const { return subdomains.size(); }

  //! \brief Set owner for a patch.
  void setPatchOwner(size_t p, size_t owner) { patchOwner[p-1] = owner; }

  //! \brief Get process owning patch p.
  int getPatchOwner(size_t p) { return patchOwner[p-1]; }

  //! \brief Get global equation number
  int getGlobalEq(int lEq) const;

  //! \brief Obtain local-to-global eq mapping.
  const std::vector<int>& getMLGEQ() const { return MLGEQ; }

  //! \brief Obtain local-to-global node mapping.
  const std::vector<int>& getMLGN() const { return MLGN; }

  //! \brief Returns associated SAM
  const SAMpatch* getSAM() const { return sam; }

private:
  std::vector<std::vector<size_t>> subdomains; //!< Number of local subdomains

  //! \brief Calculate the global node numbers for given finite element model.
  bool calcGlobalNodeNumbers(const ProcessAdm& adm, const SIMbase& sim);

  //! \brief Calculate the global equation numbers for given finite element model.
  bool calcGlobalEqNumbers(const ProcessAdm& adm, const SIMbase& sim);

  bool setupBlocks(const LinSolParams& solParams, const ProcessAdm& adm,
                   const SIMbase& sim);

  std::map<int,int> patchOwner; //!< Process that owns a particular patch

  std::vector<int> MLGN; //!< Process-local-to-global node numbers
  std::vector<int> MLGEQ; //!< Process-local-to-global equation numbers
  int minEq; //!< First equation we own
  int maxEq; //!< Last equation we own
  int minDof; //!< First DOF we own
  int maxDof; //!< Last DOF we own
  int minNode; //!< First node we own
  int maxNode; //!< Last node we own
  int nsd; //!< Number of spatial dimensions

  //! \brief Holds information about a block matrix
  struct MatrixBlock {
    int minEq; //!< First equation we own in block
    int maxEq; //!< Last equation we own in block
    std::vector<int> eqs; //!< Equations numbers belonging to matrix block on this process
    std::vector<int> MLGEQ; //!< Global equation numbers in block
  };

  std::vector<MatrixBlock> blocks; //!< Equations belonging to matrix blocks

  const SAMpatch* sam; //!< The SAM the DD is setup across
};

#endif
