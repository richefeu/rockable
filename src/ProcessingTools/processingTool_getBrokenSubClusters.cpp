//  Copyright or © or Copr. Rockable
//
//  vincent.richefeu@3sr-grenoble.fr
//
//  This software is a computer program whose purpose is
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces.
//  It is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include "processingTool_getBrokenSubClusters.hpp"

#include <map>

/**
    @brief      Get the set of sub-parts (broken clusters).
    @attention  The clusterId in subParts is NOT the original cluster number (when it was not broken).
                Let say clusterId should be called partId in 'clusterParticles'

    @param[out]   subParts   A vector of 'clusterParticles'

    @remark This method has been carfully checked with Marta. It seems to work correctly
*/
void getBrokenSubClusters(Rockable* box, std::vector<clusterParticles>& subParts) {
  subParts.clear();  // clear clusters in case it's not empty

  // Union-find: two particles of the same cluster linked by an interface (still
  // unbroken) belong to the same sub-cluster, and so on through the interfaces
  size_t n = box->Particles.size();
  std::vector<size_t> parent(n);
  for (size_t i = 0; i < n; i++) parent[i] = i;
  auto find = [&parent](size_t i) {
    while (parent[i] != i) {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  };

  for (size_t i = 0; i < box->Interfaces.size(); i++) {
    for (auto it = box->Interfaces[i].begin(); it != box->Interfaces[i].end(); ++it) {
      size_t I = it->i;
      size_t J = it->j;
      if (box->Particles[I].cluster == box->Particles[J].cluster) {
        size_t a = find(I);
        size_t b = find(J);
        if (a != b) parent[a] = b;
      }
    }
  }

  // One sub-cluster per root, for the particles that are not driven (a particle
  // without any interface is a sub-cluster on its own)
  std::map<size_t, size_t> subOfRoot;
  for (size_t p = box->nDriven; p < n; p++) {
    size_t r = find(p);
    auto found = subOfRoot.find(r);
    if (found == subOfRoot.end()) {
      clusterParticles C;
      C.clusterId = subParts.size();
      subOfRoot[r] = subParts.size();
      subParts.push_back(C);
      found = subOfRoot.find(r);
    }
    subParts[found->second].particleId.push_back(p);
  }
}
