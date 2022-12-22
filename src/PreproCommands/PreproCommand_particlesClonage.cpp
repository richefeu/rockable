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

#include "PreproCommand_particlesClonage.hpp"
#include "Core/Rockable.hpp"

particlesClonage::particlesClonage() {}

void particlesClonage::addCommand() {
  box->parser.kwMap["particlesClonage"] = [this](std::istream& conf) {
    conf >> this->ifirst >> this->ilast >> this->translation;
    exec();
  };
}

/**
   Usage in input conf-file:
   particlesClonage idFirst idLast dX dY dZ
*/
void particlesClonage::exec() {
  if (box->Particles.empty()) return;

  int numClusterMax = box->Particles[0].cluster;
  for (size_t i = 1; i < box->Particles.size(); i++) {
    if (box->Particles[i].cluster > numClusterMax) numClusterMax = box->Particles[i].cluster;
  }

  Particle P;
  for (size_t i = ifirst; i <= ilast; i++) {
    P.group = box->Particles[i].group;
    P.cluster = box->Particles[i].cluster - ifirst + numClusterMax + 1;
    P.homothety = box->Particles[i].homothety;
    P.pos = box->Particles[i].pos + translation;
    P.vel = box->Particles[i].vel;
    P.acc = box->Particles[i].acc;
    P.Q = box->Particles[i].Q;
    P.vrot = box->Particles[i].vrot;
    P.arot = box->Particles[i].arot;
    P.shape = box->Particles[i].shape;
    P.homothety = box->Particles[i].homothety;
    P.mass = box->Particles[i].mass;
    P.inertia = box->Particles[i].inertia;
    P.obb = box->Particles[i].obb;
    P.force = box->Particles[i].force;
    P.moment = box->Particles[i].moment;
    box->Particles.push_back(P);
  }

  if (box->Interactions.size() != box->Particles.size()) box->Interactions.resize(box->Particles.size());
  if (box->Interfaces.size() != box->Particles.size()) box->Interfaces.resize(box->Particles.size());
}
