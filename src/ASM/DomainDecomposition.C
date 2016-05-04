// $Id$
//==============================================================================
//!
//! \file DomainDecomposition.C
//!
//! \date Feb 23 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Domain decomposition partitioning for structured models.
//!
//==============================================================================

#include "DomainDecomposition.h"
#include "ASMbase.h"
#include "ASMstruct.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#include "SAMpatch.h"
#include "SIMbase.h"
#include "IFEM.h"
#include "Utilities.h"
#include <cassert>
#include <numeric>
#include <iostream>


//! \brief Iterator for matching nodes on edges/faces with a given orientation and index.
class NodeIterator {
public:
  //! \brief Constructor.
  //! \param pch Slave patch
  //! \param orient Orientation of boundary on slave
  //! \param lIdx Index of boundary on slave
  //! \param basis Basis to use for boundary on slave
  NodeIterator(const ASMstruct* pch, int orient, int lIdx, int basis)
  {
    int n1, n2, n3;
    pch->getSize(n1,n2,n3,basis);
    int nsd = pch->getNoSpaceDim();
    if (nsd == 3) {
      int dim1, dim2;
      if (lIdx == 1 || lIdx == 2)
        dim1 = n2, dim2 = n3;
      else if (lIdx == 3 || lIdx == 4)
        dim1 = n1, dim2 = n3;
      else //if (lIdx == 5 || lIdx == 6)
        dim1 = n1, dim2 = n2;

      nodes.resize(dim1*dim2);
      if (orient == 0)
        std::iota(nodes.begin(), nodes.end(), 0);
      else if (orient == 1) {
        auto it = nodes.begin();
        for (int n = dim2-1; n >= 0; --n, it += dim1)
          std::iota(it, it+dim1, n*dim1);
      }
      else if (orient == 2) {
        for (int n = 0; n < dim2; ++n) {
          int idx = (n+1)*dim1;
           std::generate(nodes.begin()+n*dim1, nodes.begin()+(n+1)*dim1,
                       [&idx] { return --idx; });
        }
      }
      else if (orient == 3)
        std::iota(nodes.rbegin(), nodes.rend(), 0);
      else if (orient == 4) {
        for (int n = 0; n < dim2; ++n) {
          int idx = n-dim1;
          std::generate(nodes.begin()+n*dim1, nodes.begin()+(n+1)*dim1,
                        [dim1,&idx] { return idx += dim1; });
        }
      } else if (orient == 5) {
        for (int n = 0; n < dim2; ++n) {
          int idx = -n-1;
          std::generate(nodes.begin()+n*dim1, nodes.begin()+(n+1)*dim1,
                         [dim1,&idx] { return idx += dim1; });
        }
      } else if (orient == 6) {
        for (int n = 0; n < dim2; ++n) {
          int idx = (dim2-1)*dim1+n + dim1;
          std::generate(nodes.begin()+n*dim1, nodes.begin()+(n+1)*dim1,
                        [dim1,&idx] { return idx -= dim1; });
        }
      } else if (orient == 7) {
        for (int n = 0; n < dim2; ++n) {
          int idx = dim2*dim1-n-1 + dim1;
          std::generate(nodes.begin()+n*dim1, nodes.begin()+(n+1)*dim1,
                        [dim1,&idx] { return idx -= dim1; });
        }
      }
    } else {
      if (lIdx == 1 || lIdx == 2)
        nodes.resize(n2);
      else
        nodes.resize(n1);

      if (orient == 0)
        std::iota(nodes.begin(), nodes.end(), 0);
      else
        std::iota(nodes.rbegin(), nodes.rend(), 0);
    }
  }

  //! \brief Obtain start of node numbers.
  std::vector<int>::const_iterator begin() { return nodes.begin(); }
  //! \brief Obtain end of node numbers.
  std::vector<int>::const_iterator end() { return nodes.end(); }
  //! \brief Obtain number of nodes
  size_t size() const { return nodes.size(); }
protected:
  //! \brief Node number on boundary.
  std::vector<int> nodes;
};


void DomainDecomposition::calcAppropriateGroups(size_t nel1, size_t nel2, size_t nel3,
                                                size_t g1, size_t g2, size_t g3, size_t overlap)
{
  if (nel3 > 0)
    calcGroups(nel1, nel2, nel3, g1, g2, g3, overlap);
  else if (nel2 > 0)
    calcGroups(nel1, nel2, g1, g2, overlap);
  else
    calcGroups(nel1, g1, overlap);
}


void DomainDecomposition::calcGroups (size_t nel1, size_t g1, size_t overlap)
{
  if (g1 == 1)
  {
    subdomains.resize(1);
    subdomains[0].resize(nel1);
    std::iota(subdomains[0].begin(),subdomains[0].end(),0);
  }
  else
  {
    size_t nel1_sub = floor(double(nel1)/g1 + 1);
    subdomains.resize(g1);
    size_t ofs1 = 0;
    for (size_t gu = 0; gu < g1; ++gu) {
      if (gu == g1 - 1)
        nel1_sub = nel1 - ofs1;
      subdomains[gu].reserve(nel1_sub);
      for (size_t i = 0; i < nel1_sub && ofs1 + i < nel1; ++i)
        subdomains[gu].push_back(ofs1 + i);
      ofs1 += nel1_sub - overlap;
    }
  }
}


void DomainDecomposition::calcGroups (size_t nel1, size_t nel2,
                                      size_t g1, size_t g2, size_t overlap)
{
  if (g1*g2 == 1)
  {
    subdomains.resize(1);
    subdomains[0].resize(nel1*nel2);
    std::iota(subdomains[0].begin(),subdomains[0].end(),0);
  }
  else
  {
    size_t nel2_sub = floor(double(nel2)/g2 + 1);
    subdomains.resize(g1*g2);
    size_t g = 0;
    size_t ofs2 = 0;
    for (size_t gv = 0; gv < g2; ++gv) {
      if (gv == g2 - 1)
        nel2_sub = nel2 - ofs2;
      size_t nel1_sub = floor(double(nel1)/g1 + 1);
      size_t ofs1 = 0;
      for (size_t gu = 0; gu < g1; ++gu, ++g) {
        if (gu == g1 - 1)
          nel1_sub = nel1-ofs1;

        subdomains[g].reserve(nel1_sub*nel2_sub);
        for (size_t j = 0; j < nel2_sub && ofs2 + j < nel2; ++j)
          for (size_t i = 0; i < nel1_sub && ofs1 + i < nel1; ++i)
            subdomains[g].push_back((ofs2+j)*nel1 + ofs1 + i);
        ofs1 += nel1_sub - overlap;
      }
      ofs2 += nel2_sub - overlap;
    }
  }
}


void DomainDecomposition::calcGroups (size_t nel1, size_t nel2, size_t nel3,
                                      size_t g1, size_t g2, size_t g3, size_t overlap)
{
  if (g1*g2*g3 == 1)
  {
    subdomains.resize(1);
    subdomains[0].resize(nel1*nel2*nel3);
    std::iota(subdomains[0].begin(),subdomains[0].end(),0);
  }
  else
  {
    subdomains.resize(g1*g2*g3);
    size_t g = 0;
    size_t ofs3 = 0;
    size_t nel3_sub = floor(double(nel3)/g3+1);
    for (size_t gw = 0; gw < g3; ++gw) {
      if (gw == g3 - 1)
        nel3_sub = nel3-ofs3;
      size_t ofs2 = 0;
      size_t nel2_sub = floor(double(nel2)/g2+1);
      for (size_t gv = 0; gv < g2; ++gv) {
        if (gv == g2 - 1)
          nel2_sub = nel2-ofs2;
        size_t ofs1 = 0;
        size_t nel1_sub = floor(double(nel1)/g1+1);
        for (size_t gu = 0; gu < g1; ++gu, ++g) {
          if (gu == g1 - 1)
            nel1_sub = nel1-ofs1;
          subdomains[g].reserve(nel1_sub*nel2_sub*nel3_sub);
          for (size_t k = 0; k < nel3_sub && ofs3 + k < nel3; ++k)
            for (size_t j = 0; j < nel2_sub && ofs2 + j < nel2; ++j)
              for (size_t i = 0; i < nel1_sub && ofs1 + i < nel1; ++i)
                subdomains[g].push_back((ofs3 + k)*nel1*nel2 + (ofs2+j)*nel1 + ofs1 + i);
          ofs1 += nel1_sub - overlap;
        }
        ofs2 += nel2_sub - overlap;
      }
      ofs3 += nel3_sub - overlap;
    }
  }
}


bool DomainDecomposition::calcGlobalNodeNumbers(const ProcessAdm& adm,
                                                const SIMbase& sim)
{
  minDof = minNode = 1;
  maxDof = sim.getSAM()->getNoDOFs();
  maxNode = sim.getSAM()->getNoNodes();
#ifdef HAVE_MPI
  if (adm.getNoProcs() == 1)
    return true;

  minDof = 1;
  maxDof = sim.getSAM()->getNoDOFs();
  minNode = 1;
  maxNode = sim.getSAM()->getNoNodes();
  if (adm.getProcId() > 0) {
    adm.receive(minNode, adm.getProcId()-1);
    adm.receive(minDof, adm.getProcId()-1);
    maxDof  = minDof++;
    maxNode = minNode++;
  }

  MLGN.resize(sim.getSAM()->getNoNodes());
  std::iota(MLGN.begin(), MLGN.end(), minNode);

  std::map<int,int> old2new;
  for (const auto& it : ghostConnections) {
    int sidx = sim.getLocalPatchIndex(it.slave);
    if (sidx < 1)
      continue;

    IntVec lNodes;
    sim.getPatch(sidx)->getBoundaryNodes(it.sidx, lNodes);
    int nRecv;
    adm.receive(nRecv, getPatchOwner(it.master));
    if (nRecv =! lNodes.size()) {
      std::cerr <<"\n *** DomainDecomposition::calcGlobalNodeNumbers(): Topology error, boundary size "
        << nRecv << ", expected " << lNodes.size() << std::endl;
      return false;
    }
    IntVec glbNodes(lNodes.size());
    adm.receive(glbNodes, getPatchOwner(it.master));
    size_t ofs = 0;
    for (size_t b = 1; b <= sim.getPatch(sidx)->getNoBasis(); ++b) {
      NodeIterator iter(dynamic_cast<const ASMstruct*>(sim.getPatch(sidx)),
                        it.orient, it.sidx, b);
      auto it_n = iter.begin();
      for (size_t i = 0; i < iter.size(); ++i, ++it_n) {
        int node = MLGN[lNodes[i+ofs]-1];
        old2new[node] = glbNodes[*it_n + ofs];
      }
      ofs += iter.size();
    }
  }

  // remap ghost nodes
  for (auto& it : MLGN)
    utl::renumber(it, old2new, false);

  // remap rest of our nodes
  for (int i = 0; i < sim.getSAM()->getNoNodes() && adm.getProcId() != 0; ++i)
    if (old2new.find(i + minNode) == old2new.end()) {
      std::map<int,int> old2new2;
      old2new2[i + minNode] = ++maxNode;
      auto dof = sim.getSAM()->getNodeDOFs(i);
      maxDof += dof.second-dof.first+1;
      for (auto& it : MLGN)
        utl::renumber(it, old2new2, false);
    }

  if (adm.getProcId() < adm.getNoProcs()-1) {
    adm.send(maxNode, adm.getProcId()+1);
    adm.send(maxDof, adm.getProcId()+1);
  }

  for (const auto& it : ghostConnections) {
    int midx = sim.getLocalPatchIndex(it.master);
    if (midx < 1)
      continue;

    std::vector<int> glbNodes;
    sim.getPatch(midx)->getBoundaryNodes(it.midx, glbNodes);
    for (size_t i = 0; i < glbNodes.size(); ++i)
      glbNodes[i] = MLGN[glbNodes[i]-1];

    adm.send(int(glbNodes.size()), getPatchOwner(it.slave));
    adm.send(glbNodes, getPatchOwner(it.slave));
  }
#endif

  return true;
}


std::vector<int> DomainDecomposition::setupEquationNumbers(const SIMbase& sim,
                                                           int pidx, int lidx)
{
  std::vector<IntVec> lNodes(sim.getPatch(pidx)->getNoBasis());
  std::vector<int> result;

  for (size_t block = 0; block < blocks.size(); ++block) {
    std::set<int> bases;
    if (block == 0) {
      for (size_t i = 1; i <= sim.getPatch(pidx)->getNoBasis(); ++i)
        bases.insert(i);
    } else
      bases = utl::getDigits(sim.getSolParams()->getBlock(block-1).basis);

    for (const int& basis : bases) {
      if (lNodes[basis-1].empty())
        sim.getPatch(pidx)->getBoundaryNodes(lidx, lNodes[basis-1], basis);

      std::set<int> components;
      if (block == 0 || sim.getSolParams()->getBlock(block-1).comps == 0) {
        for (size_t i = 1; i <= sim.getPatch(pidx)->getNoFields(basis); ++i)
          components.insert(i);
      } else
        components = utl::getDigits(sim.getSolParams()->getBlock(block-1).comps);

      for (size_t i = 0; i < lNodes[basis-1].size(); ++i) {
        int node = lNodes[basis-1][i];
        for (const int& dof : components) {
          int eq = sim.getSAM()->getEquation(node, dof);
          if (eq > 0) {
            if (block == 0)
              result.push_back(blocks[block].MLGEQ[eq-1]);
            else {
              auto it = blocks[block].localEqs.find(eq);
              result.push_back(blocks[block].MLGEQ[std::distance(blocks[block].localEqs.begin(), it)]);
            }
          } else
            result.push_back(0);
        }
      }
    }
  }

  return result;
}


bool DomainDecomposition::calcGlobalEqNumbers(const ProcessAdm& adm,
                                              const SIMbase& sim)
{
  // defaults in serial and non-parallel runs with MPI enabled builds.
  blocks[0].minEq = 1;
  blocks[0].maxEq = sim.getSAM()->getNoEquations();
  for (size_t i = 1; i < blocks.size(); ++i) {
    blocks[i].minEq = 1;
    blocks[i].maxEq = blocks[i].localEqs.size();
  }

#ifdef HAVE_MPI
  if (adm.getNoProcs() == 1)
    return true;

  std::vector<int> nEqs(blocks.size());
  if (adm.getProcId() > 0)
    adm.receive(nEqs, adm.getProcId()-1);

  for (size_t block = 0; block < blocks.size(); ++block) {
    size_t size = block == 0 ? sim.getSAM()->getNoEquations() :
                               blocks[block].localEqs.size();
    blocks[block].MLGEQ.resize(size);
    std::iota(blocks[block].MLGEQ.begin(),
              blocks[block].MLGEQ.end(), nEqs[block]+1);
    if (adm.getProcId() > 0) {
      blocks[block].minEq = nEqs[block]+1;
      blocks[block].maxEq = nEqs[block];
    } else {
      blocks[block].minEq = 1;
      blocks[block].maxEq = size;
    }
  }

  std::vector<std::map<int,int>> old2new(blocks.size());
  for (const auto& it : ghostConnections) {
    int sidx = sim.getLocalPatchIndex(it.slave);
    if (sidx < 1)
      continue;

    IntVec locEqs = setupEquationNumbers(sim, sidx, it.sidx);

    int nRecv;
    adm.receive(nRecv, getPatchOwner(it.master));
    if (nRecv =! locEqs.size()) {
      std::cerr <<"\n *** DomainDecomposition::calcGlobalNodeNumbers(): Topology error, number of equations "
        << nRecv << ", expected " << locEqs.size() << std::endl;
      return false;
    }

    IntVec glbEqs(locEqs.size());
    adm.receive(glbEqs, getPatchOwner(it.master));

    size_t ofs = 0;
    for (size_t block = 0; block < blocks.size(); ++block) {
      std::set<int> bases;
      if (block == 0) {
        for (size_t i = 1; i <= sim.getPatch(sidx)->getNoBasis(); ++i)
          bases.insert(i);
      } else
        bases = utl::getDigits(sim.getSolParams()->getBlock(block-1).basis);

      for (const int& basis : bases) {
        std::set<int> components;
        if (block == 0 || sim.getSolParams()->getBlock(block-1).comps == 0) {
          for (size_t i = 1; i <= sim.getPatch(sidx)->getNoFields(basis); ++i)
            components.insert(i);
        } else
          components = utl::getDigits(sim.getSolParams()->getBlock(block-1).comps);

        NodeIterator iter(dynamic_cast<const ASMstruct*>(sim.getPatch(sidx)),
                          it.orient, it.sidx, basis);
        auto it_n = iter.begin();
        int nodeDofs = components.size();
        for (size_t i = 0; i < iter.size(); ++i, ++it_n) {
          for (size_t comp = 0; comp < components.size(); ++comp) {
            int leq = locEqs[ofs+i*nodeDofs+comp];
            int geq = glbEqs[ofs+(*it_n)*nodeDofs+comp];
            if (leq < 1 && geq > 0) {
              std::cerr <<"\n *** DomainDecomposition::calcGlobalNodeNumbers(): "
                        << "Topology error on process " << adm.getProcId()
                        << ", dof constraint mismatch "
                        << "local: " << leq << " global " << geq << std::endl;
              return false;
            }
            if (leq < 1)
              continue;

            old2new[block][leq] = geq;
          }
        }
        ofs += iter.size()*nodeDofs;
      }
    }
  }

  for (size_t block = 0; block < blocks.size() && adm.getProcId() > 0; ++block) {
    // remap ghost equations
    for (auto& it : blocks[block].MLGEQ)
      utl::renumber(it, old2new[block], false);

    // remap the rest of our equations
    size_t size = block == 0 ? sim.getSAM()->getNoEquations() :
                               blocks[block].localEqs.size();
    for (size_t i = 1; i <= size; ++i) {
      if (old2new[block].find(i + nEqs[block]) == old2new[block].end()) {
        std::map<int,int> old2new2;
        old2new2[i + nEqs[block]] = ++blocks[block].maxEq;
        for (auto& it : blocks[block].MLGEQ)
          utl::renumber(it, old2new2, false);
      }
    }
  }

  if (adm.getProcId() < adm.getNoProcs()-1) {
    std::vector<int> maxEqs;
    for (auto& it : blocks)
      maxEqs.push_back(it.maxEq);

    adm.send(maxEqs, adm.getProcId()+1);
  }

  for (const auto& it : ghostConnections) {
    int midx = sim.getLocalPatchIndex(it.master);
    if (midx < 1)
      continue;

    IntVec glbEqs = setupEquationNumbers(sim, midx, it.midx);
    adm.send(int(glbEqs.size()), getPatchOwner(it.slave));
    if (glbEqs.size())
      adm.send(glbEqs, getPatchOwner(it.slave));
  }
#endif

  return true;
}


int DomainDecomposition::getGlobalEq(int lEq, size_t idx) const
{
  if (lEq < 1 || idx >= blocks.size())
    return 0;

  if (blocks[0].MLGEQ.empty()) { // serial
    if (lEq > blocks[idx].maxEq)
      return 0;

    return lEq;
  }

  if (!blocks[idx].MLGEQ.empty() && lEq > (int)blocks[idx].MLGEQ.size())
    return 0;

  return blocks[idx].MLGEQ[lEq-1];
}


bool DomainDecomposition::setup(const ProcessAdm& adm, const SIMbase& sim)
{
//  nsd = sim.getNoSpaceDim();
  sam = dynamic_cast<const SAMpatch*>(sim.getSAM());

  // Establish global node numbers
  if (!calcGlobalNodeNumbers(adm, sim))
    return false;

  // Establish local equation mappings for each block.
  if (sim.getSolParams() && sim.getSolParams()->getNoBlocks() > 1) {
    const LinSolParams& solParams = *sim.getSolParams();
    blocks.resize(solParams.getNoBlocks()+1);
    
    // Find local equations for each block
    for (size_t i = 0; i < solParams.getNoBlocks(); ++i) {
      // grab DOFs of given type(s)
      int basis = solParams.getBlock(i).basis;
      char dofType = basis  == 1 ? 'D' : 'P'+basis-2;
      if (solParams.getBlock(i).comps != 0) {
        std::set<int> comps = utl::getDigits(solParams.getBlock(i).comps);
        for (auto& c : comps) {
          std::set<int> tmp = sam->getEquations(dofType, c);
          blocks[i+1].localEqs.insert(tmp.begin(), tmp.end());
        }
      } else {
        std::set<int> bases = utl::getDigits(basis);
        for (auto& b : bases) {
          int cb = b;
          dofType = cb == 1 ? 'D' : 'P'+cb-2;
          std::set<int> tmp = adm.dd.getSAM()->getEquations(dofType);
          blocks[i+1].localEqs.insert(tmp.begin(), tmp.end());
          // hack. stick multipliers in the second block. Correct thing to do for average pressure constraint in Stokes.
          if (i == 1) {
            tmp = sam->getEquations('L', 1);
            blocks[i+1].localEqs.insert(tmp.begin(), tmp.end());
          }
        }
      }
    }
  }

  // Establish global equation numbers for all blocks.
  if (!calcGlobalEqNumbers(adm, sim))
    return false;

  return true;
}


int DomainDecomposition::getPatchOwner(size_t p) const
{
  auto it = patchOwner.find(p);
  if (it == patchOwner.end())
    return -1;

  return it->second;
}
