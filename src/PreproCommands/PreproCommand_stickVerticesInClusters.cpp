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

#include "PreproCommand_stickVerticesInClusters.hpp"
#include "Core/Rockable.hpp"

StickVerticesInClusters::StickVerticesInClusters() {}

void StickVerticesInClusters::addCommand() {
  box->parser.kwMap["stickVerticesInClusters"] = [this](std::istream& conf) {
    conf >> this->epsilonDist;
    exec();
  };
}

/**
   @brief Create sticked contacts between vertices (spheres) of particles
          that belong to the same cluster (same cluster-ID)

   @param[in]  epsilonDist  A small distance above which the 'glue' is not added
*/
void StickVerticesInClusters::exec() {
  for (size_t i = 0; i < box->Particles.size(); i++) {
    box->Particles[i].updateObb();
  }

  for (size_t i = box->nDriven; i < box->Particles.size(); i++) {

    OBB obbi = box->Particles[i].obb;
    obbi.enlarge(0.5 * box->DVerlet);

    for (size_t j = i + 1; j < box->Particles.size(); j++) {

      if (box->Particles[i].cluster != box->Particles[j].cluster) continue;

      OBB obbj = box->Particles[j].obb;
      obbj.enlarge(0.5 * box->DVerlet);

      // Check intersection
      if (obbi.intersect(obbj)) {

        for (size_t iv = 0; iv < box->Particles[i].shape->vertex.size(); iv++) {
          vec3r Vi = box->Particles[i].GlobVertex(iv);
          double Ri = box->Particles[i].MinskowskiRadius();
          for (size_t jv = 0; jv < box->Particles[j].shape->vertex.size(); jv++) {
            vec3r Vj = box->Particles[j].GlobVertex(jv);
            double Rj = box->Particles[j].MinskowskiRadius();

            double distSqr = norm2(Vj - Vi);
            double dMaxSqr = epsilonDist + (Ri + Rj);
            dMaxSqr *= dMaxSqr;

            if (distSqr < dMaxSqr) {
              // It is necessarily two free bodies (spheres at vertices) that interact
              double meff =
                  (box->Particles[i].mass * box->Particles[j].mass) / (box->Particles[i].mass + box->Particles[j].mass);
              double en2 = box->dataTable.get(box->idEn2InnerBond, box->Particles[i].group, box->Particles[j].group);
              double kn = box->dataTable.get(box->idKnInnerBond, box->Particles[i].group, box->Particles[j].group);
              double Damp = 0.0;
              if (en2 > 0.0 && en2 < 1.0) {
                double logen = 0.5 * log(en2);
                double dampRate = -logen / sqrt(logen * logen + Mth::piSqr);
                Damp = dampRate * 2.0 * sqrt(kn * meff);
              } else if (en2 <= 0.0)
                Damp = 2.0 * sqrt(kn * meff);
              else
                Damp = 0.0;

              BreakableInterface BI_toInsert(i, j);
              BI_toInsert.isInner = 1;
              BI_toInsert.kn = kn;
              BI_toInsert.kt = box->dataTable.get(box->idKtInnerBond, box->Particles[i].group, box->Particles[j].group);
              BI_toInsert.kr = 0.0;
              box->dataTable.set(box->idKrInnerBond, box->Particles[i].group, box->Particles[j].group, 0.0);
              BI_toInsert.fn0 =
                  box->dataTable.get(box->idFn0InnerBond, box->Particles[i].group, box->Particles[j].group);
              BI_toInsert.ft0 =
                  box->dataTable.get(box->idFt0InnerBond, box->Particles[i].group, box->Particles[j].group);
              BI_toInsert.mom0 = 1.0e10;
              box->dataTable.set(box->idMom0InnerBond, box->Particles[i].group, box->Particles[j].group, 1.0e10);
              BI_toInsert.power =
                  box->dataTable.get(box->idPowInnerBond, box->Particles[i].group, box->Particles[j].group);
              std::pair<std::set<BreakableInterface>::iterator, bool> ret;
              ret = box->Interfaces[i].insert(BI_toInsert);
              BreakableInterface* BI = const_cast<BreakableInterface*>(std::addressof(*(ret.first)));
              Interaction I(i, j, vvType, iv, jv, Damp, BI);
              Interaction::UpdateDispatcher[vvType](I, box->Particles[i], box->Particles[j]);
              Interaction* Iptr = const_cast<Interaction*>(std::addressof(*((box->Interactions[i].insert(I)).first)));
              BI->concernedBonds.push_back(Iptr);
              box->activeInteractions.push_back(Iptr);
            }
          }
        }

      }  // if obb intersect
    }    // j
  }      // i

  std::sort(box->activeInteractions.begin(), box->activeInteractions.end(), std::less<Interaction*>());
}
