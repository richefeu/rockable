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

#include "PreproCommand_homothetyRange.hpp"
#include "Rockable.hpp"

//static Registrar<PreproCommand, homothetyRange> registrar("homothetyRange");

homothetyRange::homothetyRange() {}

void homothetyRange::addCommand() {
  box->parser.kwMap["homothetyRange"] = [this](std::istream& conf) {
    int timeSeededInt;
    conf >> this->ifirst >> this->ilast >> this->hmin >> this->hmax >> timeSeededInt;
    this->timeSeeded = (bool)timeSeededInt;
    exec();
  };
}

/**
   @brief Set Uniform distribution of homothety of a sub-set of particles
   @param[in]  idFirst     ID-Number of the first particle
   @param[in]  idLast      ID-Number of the last particle
   @param[in]  hmin        Minimum of the homothety range
   @param[in]  hmax        Maximum of the homothety range
   @param[in]  timeSeeded  if true, the random generator is seeded with current time

   Usage in input conf-file:
   homothetyRange idFirst idLast hmin hmax timeSeeded(0/1)
*/
void homothetyRange::exec() {
  std::default_random_engine generator;
  if (timeSeeded == true) {
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }
  std::uniform_real_distribution<> distribution(hmin, hmax);

  for (size_t i = ifirst; i <= ilast; i++) {
    double h = distribution(generator);
    box->Particles[i].homothety = h;
    box->Particles[i].mass =
        (h * h * h * box->Particles[i].shape->volume) * box->properties.get(box->idDensity, box->Particles[i].group);
    box->Particles[i].inertia = (h * h * box->Particles[i].shape->inertia_mass) * box->Particles[i].mass;
  }
}
