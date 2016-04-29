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
#include <set>
#include <vector>
#include <cstddef>

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

  //! \brief Functor to order ghost connections.
  class SlaveOrder {
    public:
      //! \brief Compare Interface by master index, if equal by slave index.
      bool operator()(const Interface& A, const Interface& B) const
      {
        return A.master != B.master ? A.master < B.master : A.slave < B.slave;
      }
  };

  std::set<Interface, SlaveOrder> ghostConnections; //!< Connections to other processes.

  DomainDecomposition() : blocks(1) {}

  //! \brief Get number of spatial dimensions
  int getNoSpaceDim() const { return nsd; }

  //! \brief Setup domain decomposition
  bool setup(const ProcessAdm& adm, const SIMbase& sim);

  //! \brief Get first equation owned by this process.
  int getMinEq(size_t idx = 0) const { return blocks[idx].minEq; };
  //! \brief Get last equation owned by this process.
  int getMaxEq(size_t idx = 0) const { return blocks[idx].maxEq; };
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
  void setPatchOwner(size_t p, size_t owner) { patchOwner[p] = owner; }

  //! \brief Get process owning patch p.
  int getPatchOwner(size_t p) const;

  //! \brief Get global equation number
  //! \param lEq Local equation number
  //! \param idx Block to get index for
  int getGlobalEq(int lEq, size_t idx=0) const;

  //! \brief Obtain local-to-global equation mapping.
  const std::vector<int>& getMLGEQ(size_t idx = 0) const { return blocks[idx].MLGEQ; }

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

  std::map<int,int> patchOwner; //!< Process that owns a particular patch

  //! \brief Struct with information per matrix block.
  struct BlockInfo {
    std::vector<int> MLGEQ; //!< Process-local-to-global equation numbers for block.
    int minEq; //!< First equation we own in block.
    int maxEq; //!< Last equation we own in block.
    std::set<int> localEqs; //!< Local equations belonging to the block.
  };

  std::vector<int> MLGN; //!< Process-local-to-global node numbers
  std::vector<BlockInfo> blocks; //!< Equation mappings for all matrix blocks. First entry is for the total matrix.
  int minDof; //!< First DOF we own
  int maxDof; //!< Last DOF we own
  int minNode; //!< First node we own
  int maxNode; //!< Last node we own
  int nsd; //!< Number of spatial dimensions

  const SAMpatch* sam; //!< The SAM the DD is setup across
};

#endif
