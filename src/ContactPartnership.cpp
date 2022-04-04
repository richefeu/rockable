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

#include <iostream>
#include <utility>

#include "message.hpp"

#include "ContactPartnership.hpp"
#include "Rockable.hpp"

ContactPartnership::ContactPartnership() : name("None"), update(nullptr), getWeight(nullptr) {}
ContactPartnership::~ContactPartnership() {}

void ContactPartnership::setModel(std::string& modelName) {
  auto log = spdlog::get("console");
  if (modelName == "NumberWeight") {
    name = modelName;
    log->info("ContactPartnership has been set to 'NumberWeight'");

    update = [this](Rockable& box) -> void {
      weightMap.clear();

      // FIXME: The following cannot be easily parallelized because of
      //        concurrent access to 'weightMap'

      //#pragma omp parallel for default(shared)
      for (size_t k = 0; k < box.Interactions.size(); ++k) {
        for (auto it = box.Interactions[k].begin(); it != box.Interactions[k].end(); ++it) {
          Interaction* I = const_cast<Interaction*>(std::addressof(*it));
          if (I->stick != nullptr || I->dn >= 0.0) continue;
          size_t i = I->i;
          size_t j = I->j;
          if (j < i) std::swap(i, j);
          std::pair<size_t, size_t> p(i, j);

          auto ref = weightMap.find(p);
          if (ref != weightMap.end()) {
            //#pragma omp atomic
            ref->second += 1.0;
          } else {
            //#pragma omp critical(addInMap)
            weightMap[p] = 1.0;
          }
        }
      }

      // todo: replace by inverse so that the getWeight is faster
    };

    getWeight = [this](Interaction& I) -> double {
      size_t i = I.i;
      size_t j = I.j;
      if (j < i) std::swap(i, j);
      std::pair<size_t, size_t> p(i, j);
      auto ref = weightMap.find(p);
      if (ref != weightMap.end()) {
        return (1.0 / ref->second);
      }
      return 1.0;
    };

  } else if (modelName == "OverlapWeight") {
    name = modelName;
    log->info("ContactPartnership has been set to 'OverlapWeight'");

    update = [this](Rockable& box) -> void {
      weightMap.clear();

      //#pragma omp parallel for default(shared)
      for (size_t k = 0; k < box.Interactions.size(); ++k) {
        for (auto it = box.Interactions[k].begin(); it != box.Interactions[k].end(); ++it) {
          Interaction* I = const_cast<Interaction*>(std::addressof(*it));
          if (I->stick != nullptr || I->dn >= 0.0) continue;
          size_t i = I->i;
          size_t j = I->j;
          if (j < i) std::swap(i, j);
          std::pair<size_t, size_t> p(i, j);

          auto ref = weightMap.find(p);
          if (ref != weightMap.end()) {
            //#pragma omp atomic
            ref->second += I->dn;
          } else {
            //#pragma omp critical(addInMap)
            weightMap[p] = I->dn;
          }
        }
      }

      // todo: replace by inverse so that the getWeight is faster
    };

    getWeight = [this](Interaction& I) -> double {
      size_t i = I.i;
      size_t j = I.j;
      if (j < i) std::swap(i, j);
      std::pair<size_t, size_t> p(i, j);
      auto ref = weightMap.find(p);
      if (ref != weightMap.end()) {
        return (I.dn / ref->second);
      }
      return 1.0;
    };

  } else if (modelName == "SurfaceWeight") {
    name = modelName;
    log->info("ContactPartnership has been set to 'SurfaceWeight'");
    ///////////// THAT ONE IS NOT YET IMPLEMENTED ////////////
    ///////////// and maybe it will never be!     ////////////
  } else if (modelName == "None") {
    name = modelName;
    log->info("No Partnership has been set");
    update = nullptr;
    getWeight = nullptr;
  } else {
    log->info("ContactPartnership '{}' is unknown. So no Partnership has been set", modelName);
    modelName = "Unknown";
    update = nullptr;
    getWeight = nullptr;
  }
}
