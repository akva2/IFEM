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
#include "ProcessAdm.h"
#include "SAM.h"
#include "SIMbase.h"
#include "Utilities.h"
#include <numeric>
#include <iostream>


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
#ifdef PARALLEL_PETSC
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
  for (int i = 0; i < sim.getSAM()->getNoNodes();++i)
    MLGN[i] = i + minNode;

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
    for (size_t i = 0; i < glbNodes.size(); ++i) {
      int node = MLGN[lNodes[it.orient==1?glbNodes.size()-1-i:i]-1];
      old2new[node] = glbNodes[i];
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
    adm.send(maxDof, adm.getProcId()+1);
    adm.send(maxNode, adm.getProcId()+1);
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

//  for (int p = 0; p < adm.getNoProcs(); ++p) {
//    std::flush(std::cout);
//    MPI_Barrier(PETSC_COMM_WORLD);
//    if (adm.getProcId() == p) {
//      std::cout << "Min/max node " << minNode << " " << maxNode << " dof " << minDof << " " << maxDof << std::endl;
//      for (const auto& it : MLGN)
//        std::cout << it << " ";
//      std::cout << std::endl;
//    }
//  }
#else
  minDof = minNode = 1;
  maxDof = sim.getSAM()->getNoDOFs();
  maxNode = sim.getSAM()->getNoNodes();
#endif

  return true;
}


bool DomainDecomposition::calcGlobalEqNumbers(const ProcessAdm& adm,
                                              const SIMbase& sim)
{
#ifdef PARALLEL_PETSC
  int nEq = 0;
  if (adm.getProcId() > 0)
    adm.receive(nEq, adm.getProcId()-1);
  else
    adm.send(sim.getSAM()->getNoEquations(), 1);

  MLGEQ.resize(sim.getSAM()->getNoEquations());
  for (int i = 0; i < sim.getSAM()->getNoEquations();++i)
    MLGEQ[i] = i + 1 + nEq;

  minEq = nEq+1;
  maxEq = nEq;

  std::map<int,int> old2new;
  for (const auto& it : ghostConnections) {
    int sidx = sim.getLocalPatchIndex(it.slave);
    if (sidx < 1)
      continue;

    IntVec lNodes;
    sim.getPatch(sidx)->getBoundaryNodes(it.sidx, lNodes);
    std::vector<int> locEqs;
    for (size_t i = 0; i < lNodes.size(); ++i) {
      int node = lNodes[i];
      std::pair<int,int> dofs = sim.getSAM()->getNodeDOFs(node);
      for (int dof = dofs.first; dof <= dofs.second; ++dof) {
        int eq = sim.getSAM()->getEquation(node, dof-dofs.first+1);
        if (eq > 0)
          locEqs.push_back(MLGEQ[eq-1]);
        else
          locEqs.push_back(0);
      }
    }
    int nRecv;
    adm.receive(nRecv, getPatchOwner(it.master));
    if (nRecv =! locEqs.size()) {
      std::cerr <<"\n *** DomainDecomposition::calcGlobalNodeNumbers(): Topology error, number of equations "
        << nRecv << ", expected " << locEqs.size() << std::endl;
      return false;
    }

    IntVec glbEqs(locEqs.size());
    adm.receive(glbEqs, getPatchOwner(it.master));

    // check that dof types match
    for (size_t i = 0; i < locEqs.size(); ++i)
      if (locEqs[i] < 1 && glbEqs[i] > 0) {
        std::cerr <<"\n *** DomainDecomposition::calcGlobalNodeNumbers(): Topology error, dof constraint mismatch" << std::endl;
        return false;
      }

    for (size_t i = 0; i < glbEqs.size(); ++i) {
      if (locEqs[i] < 1)
        continue;
      int leq = locEqs[it.orient==1?glbEqs.size()-1-i:i];
      old2new[leq] = glbEqs[i];
    }
  }

  // remap ghost equations
  for (auto& it : MLGEQ)
    utl::renumber(it, old2new, false);

  // remap the rest of our equations
  for (int i = 1; i <= sim.getSAM()->getNoEquations(); ++i)
    if (old2new.find(i + nEq) == old2new.end()) {
      std::map<int,int> old2new2;
      old2new2[i + nEq] = ++maxEq;
      for (auto& it : MLGEQ)
        utl::renumber(it, old2new2, false);
    }

  if (adm.getProcId() < adm.getNoProcs()-1 && adm.getProcId() != 0)
    adm.send(maxEq, adm.getProcId()+1);

  for (const auto& it : ghostConnections) {
    int midx = sim.getLocalPatchIndex(it.master);
    if (midx < 1)
      continue;

    std::vector<int> glbEqs;
    std::vector<int> glbNodes;
    sim.getPatch(midx)->getBoundaryNodes(it.midx, glbNodes);
    for (size_t i = 0; i < glbNodes.size(); ++i) {
      int node = glbNodes[i];
      std::pair<int,int> dofs = sim.getSAM()->getNodeDOFs(node);
      for (int dof = dofs.first; dof <= dofs.second; ++dof) {
        int eq = sim.getSAM()->getEquation(node, dof-dofs.first+1);
        if (eq > 0)
          glbEqs.push_back(MLGEQ[eq-1]);
        else
          glbEqs.push_back(0);
      }
    }
    adm.send(int(glbEqs.size()), getPatchOwner(it.slave));
    if (glbEqs.size())
      adm.send(glbEqs, getPatchOwner(it.slave));
  }

//  for (int p = 0; p < adm.getNoProcs(); ++p) {
//    std::flush(std::cout);
//    MPI_Barrier(PETSC_COMM_WORLD);
//    if (adm.getProcId() == p) {
//      std::cout << "proc " << p <<": ";
//      for (const auto& it : MLGEQ)
//        std::cout << it << " ";
//      std::cout << std::endl;
//    }
//  }
//  exit(1);

#else
  minEq = 1;
  maxEq = sim.getSAM()->getNoEquations();
#endif
  return true;
}


int DomainDecomposition::getGlobalEq(int lEq) const
{
  if (lEq < 1 || (!MLGEQ.empty() && lEq > (int)MLGEQ.size()))
    return 0;

  if (MLGEQ.empty()) {
    if (lEq > maxEq)
      return 0;

    return lEq;
  }

  return MLGEQ[lEq-1];
}


bool DomainDecomposition::setup(const ProcessAdm& adm, const SIMbase& sim)
{
  return calcGlobalNodeNumbers(adm, sim) && calcGlobalEqNumbers(adm, sim);
}


int DomainDecomposition::getPatchOwner(size_t p) const
{
  auto it = patchOwner.find(p);
  if (it == patchOwner.end())
    return -1;

  return it->second;
}
