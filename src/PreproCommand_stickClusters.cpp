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

#include "factory.hpp"

#include "PreproCommand_stickClusters.hpp"
#include "Rockable.hpp"

static Registrar<PreproCommand, StickClusters> registrar("stickClusters");

StickClusters::StickClusters() { }

void StickClusters::addCommand() {  
  box->parser.kwMap["stickClusters"] = [this](std::istream& conf) {
    conf >> this->epsilonDist;
    exec();
  };
}

/**
   @brief Create sticked interface between two sub-elements of two distinct particles
          that do not belong to the same cluster (different cluster-ID).
          The idea is to update the neighbor list and then replace the contacts
          between different clusters by a SINGLE sticked link.

   @param[in]  epsilonDist  A small distance below which the 'glue' is not added.
*/
void StickClusters::exec() {
  // In case the neighbor list has not been yet updated
  box->UpdateNL();

  // Associate pair(i, j) with vector of interactions
  std::map<std::pair<size_t, size_t>, std::vector<Interaction*>> ctc_packets;

  // Find Packets of contacts
  for (size_t k = 0; k < box->Interactions.size(); ++k) {
    for (auto it = box->Interactions[k].begin(); it != box->Interactions[k].end(); ++it) {
      size_t i = it->i;
      size_t j = it->j;

      if (box->glue_with_walls == false) {
        if (i < box->nDriven || j < box->nDriven) continue;  // no sticked link involving a driven body ('walls')
      }

      // Stick is only possible in-between different clusters
      if (box->Particles[i].cluster == box->Particles[j].cluster) continue;

      // update the interaction
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      bool valid = Interaction::UpdateDispatcher[it->type](*I, box->Particles[i], box->Particles[j]);

      if (valid && it->dn < epsilonDist) {
        if (j < i) std::swap(i, j);  // NORMALLY, NOT NECESSARY BECAUSE THE NEIGHBOR LIST IS CONSTRUCTED SO THAT i < j
        std::pair<size_t, size_t> duo(i, j);
        auto itf = ctc_packets.find(duo);
        if (itf == ctc_packets.end()) {  // if not found
          std::vector<Interaction*> v;
          auto p = ctc_packets.insert(std::pair<std::pair<size_t, size_t>, std::vector<Interaction*>>(duo, v));
          itf = p.first;
        }
        itf->second.push_back(I);
      }
    }
  }

  if (ctc_packets.empty()) {
    box->console->info("@PreproCommand::stickClusters, No possible glued points, ctc_packets is empty");
    return;
  }

  // Replace contact packets by single sticked link:
  for (auto p : ctc_packets) {

    // -- barycenter of each packet
    vec3r c;
    size_t nb = p.second.size();
    for (size_t ip = 0; ip < nb; ip++) {
      c += (p.second)[ip]->pos;
    }
    c /= nb;  // remark: we know that nb >= 1

    // -- select the contact point that is the closest (for each packet)
    double d2min = norm2((p.second)[0]->pos - c);
    size_t imin = 0;
    for (size_t ip = 1; ip < nb; ip++) {
      double d2 = norm2((p.second)[ip]->pos - c);
      if (d2 < d2min) {
        d2min = d2;
        imin = ip;
      }
    }

    //  -- Create the interface and plug it to the selected interaction (it will then be 'glued')
    size_t i = (p.second)[imin]->i;
    size_t j = (p.second)[imin]->j;
    if (j < i) std::swap(i, j);
    BreakableInterface BI_toInsert(i, j);
    BI_toInsert.isInner = 0;
    BI_toInsert.dn0 = (p.second)[imin]->dn;
    (p.second)[imin]->mom.reset();
    (p.second)[imin]->ft.reset();
    (p.second)[imin]->fn = 0.0;
    BI_toInsert.kn = box->dataTable.get(box->idKnOuterBond, box->Particles[i].group, box->Particles[j].group);
    BI_toInsert.kt = box->dataTable.get(box->idKtOuterBond, box->Particles[i].group, box->Particles[j].group);
    BI_toInsert.kr = box->dataTable.get(box->idKrOuterBond, box->Particles[i].group, box->Particles[j].group);
    BI_toInsert.fn0 = box->dataTable.get(box->idFn0OuterBond, box->Particles[i].group, box->Particles[j].group);
    BI_toInsert.ft0 = box->dataTable.get(box->idFt0OuterBond, box->Particles[i].group, box->Particles[j].group);
    BI_toInsert.mom0 = box->dataTable.get(box->idMom0OuterBond, box->Particles[i].group, box->Particles[j].group);
    BI_toInsert.power = box->dataTable.get(box->idPowOuterBond, box->Particles[i].group, box->Particles[j].group);

    std::pair<std::set<BreakableInterface>::iterator, bool> ret;
    ret = box->Interfaces[i].insert(BI_toInsert);
    BreakableInterface* BI = const_cast<BreakableInterface*>(std::addressof(*(ret.first)));

    Interaction* Iptr = const_cast<Interaction*>(std::addressof(*((p.second)[imin])));

    Iptr->stick = BI;
    BI->concernedBonds.push_back(Iptr);
  }

  // -- remove 'not-glued' contacts between clusters
  for (size_t k = 0; k < box->Interactions.size(); ++k) {
    for (auto it = box->Interactions[k].begin(); it != box->Interactions[k].end();) {
      size_t i = it->i;
      size_t j = it->j;
      if (i < box->nDriven || j < box->nDriven) {
        ++it;
        continue;
      }  // contact with walls are not erased
      if (box->Particles[i].cluster == box->Particles[j].cluster) {
        ++it;
        continue;
      }  // inner bonds are not erased
      if (it->stick != nullptr) {
        ++it;
        continue;
      }  // bonded links are not erased

      // erase the interaction (when not sticked)
      it = box->Interactions[k].erase(it);  // 'it' points now to the next iterator
    }
  }

  // activeInteractions needs to be re-set because 'saveConf' uses it
  box->activeInteractions.clear();
  for (size_t k = 0; k < box->Interactions.size(); ++k) {
    for (auto it = box->Interactions[k].begin(); it != box->Interactions[k].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      if (it->dn < 0.0 || it->stick != nullptr) {
        box->activeInteractions.push_back(I);
      }
    }
  }
}
