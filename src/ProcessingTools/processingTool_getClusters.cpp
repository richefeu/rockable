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

#include "processingTool_getClusters.hpp"

/**
   @brief Get the set of clusters according to the cluster-ID that each particle has.
          Each cluster holds a vector of particle numbers (index in the vector @c Rockable::Particles)

   @param[out]  clusters  A vector of 'clusterParticles' that hold the clusterId and the list of involved particleId
*/
void getClusters(Rockable *box, std::vector<clusterParticles>& clusters) {
  clusters.clear();  // clear clusters if not empty
  std::set<clusterParticles> clusterSet;
  clusterParticles C;
  for (size_t i = box->nDriven; i < box->Particles.size(); i++) {
    C.clusterId = box->Particles[i].cluster;
    auto itClust = clusterSet.find(C);
    if (itClust == clusterSet.end()) {  // if not found
      auto p = clusterSet.insert(C);
      itClust = p.first;
    }

    // An element in the set clusterSet cannot be modified (because a set is
    // ordered) so we get a pointer to it, and change only the parameters that
    // will not affect the order
    clusterParticles* cp = const_cast<clusterParticles*>(std::addressof(*itClust));
    cp->particleId.push_back(i);
  }

  clusters.assign(clusterSet.begin(), clusterSet.end());

#if 0
  for (size_t c = 0 ; c < clusters.size() ; c++) {
    std::cout << std::endl;
    __SHOW(c);
    __SHOW( clusters[c].particleId.size() );
    for (size_t i = 0 ; i < clusters[c].particleId.size() ; i++) {
      __SHOW( clusters[c].particleId[i] );
    }
  }
#endif
}
