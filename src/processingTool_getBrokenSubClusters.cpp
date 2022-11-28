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

/**
    @brief      Get the set of sub-parts (broken clusters).
    @attention  The clusterId in subParts is NOT the original cluster number (when it was not broken).
                Let say clusterId should be called partId in 'clusterParticles'

    @param[out]   subParts   A vector of 'clusterParticles'

    @remark This method has been carfully checked with Marta. It seems to work correctly
*/
void getBrokenSubClusters(Rockable *box, std::vector<clusterParticles>& subParts) {
  subParts.clear();  // clear clusters in case it's not empty

  std::vector<std::set<int>> subs;
  std::set<int> addedParticles;
  for (size_t i = 0; i < box->Interfaces.size(); i++) {
    for (auto it = box->Interfaces[i].begin(); it != box->Interfaces[i].end(); ++it) {
      size_t I = it->i;
      size_t J = it->j;
      if (box->Particles[I].cluster == box->Particles[J].cluster) {

        int afound = -1;
        if (!subs.empty()) {
          for (size_t a = subs.size(); a-- > 0;) {
            auto itI = subs[a].find(I);
            if (itI != subs[a].end()) {
              afound = a;
              break;
            }
            auto itJ = subs[a].find(J);
            if (itJ != subs[a].end()) {
              afound = a;
              break;
            }
          }
        }

        if (afound >= 0) {
          subs[afound].insert(I);
          subs[afound].insert(J);
          addedParticles.insert(I);
          addedParticles.insert(J);
        } else {
          std::set<int> s;
          s.insert(I);
          s.insert(J);
          subs.push_back(s);
          addedParticles.insert(I);
          addedParticles.insert(J);
        }

      }  // end 'if same cluster'
    }    // end 'for it'
  }      // end 'for i'

  // Merge sets that have common particle.
  // Remark: the merge of more than 2 sets is not possible because an interface can onle involve 2 particles
  std::vector<std::set<int>> subsFinal;
  for (size_t a = 0; a < subs.size(); a++) {
    bool hasBeenMerged = false;
    for (size_t b = a + 1; b < subs.size(); b++) {
      bool hasCommon = false;
      for (auto it = subs[a].begin(); it != subs[a].end(); ++it) {
        if (subs[b].find(*it) != subs[b].end()) {  // found
          hasCommon = true;
          break;
        }
      }
      if (hasCommon == true) {
        std::set<int> s;
        for (auto ia = subs[a].begin(); ia != subs[a].end(); ++ia) s.insert(*ia);
        for (auto ib = subs[b].begin(); ib != subs[b].end(); ++ib) s.insert(*ib);
        subsFinal.push_back(s);
        subs.erase(subs.begin() + b);
        hasBeenMerged = true;
        break;
      }
    }

    if (hasBeenMerged == false) {
      // add a to subsFinal
      std::set<int> s;
      for (auto ia = subs[a].begin(); ia != subs[a].end(); ++ia) s.insert(*ia);
      subsFinal.push_back(s);
    }
  }

  // Add parts composed of single particle
  for (int n = box->nDriven; n < (int)box->Particles.size(); n++) {
    if (addedParticles.find(n) == addedParticles.end()) {
      std::set<int> s;
      s.insert(n);
      subsFinal.push_back(s);
    }
  }

  for (size_t a = 0; a < subsFinal.size(); a++) {
    clusterParticles C;
    C.clusterId = a;
    C.particleId.assign(subsFinal[a].begin(), subsFinal[a].end());
    subParts.push_back(C);
  }

#if 0
  std::cout << " --------- " << std::endl;
  for (size_t c = 0 ; c < subParts.size() ; c++) {
    std::cout << std::endl;
    __SHOW(c);
    __SHOW( subParts[c].particleId.size() );
    for (size_t i = 0 ; i < subParts[c].particleId.size() ; i++) {
      __SHOW( subParts[c].particleId[i] );
    }
  }
#endif
}
