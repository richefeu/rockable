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

#include <limits>

#include "factory.hpp"

#include "Core/Rockable.hpp"
#include "DataExtractor_TrackStress.hpp"

// Component xx of the stress (body frame) of the soft particles, from the particle istart to the last one
// (written by Mukesh Singh Bisht). The rigid particles are recorded with a zero stress.
TrackStress::TrackStress() : istart(0) {}

void TrackStress::read(std::istream& is) {
  is >> istart >> filename >> nrec;
  if (box->isInteractive() == false) recordFile.open(filename.c_str());
  nstep = std::numeric_limits<int>::max();

#ifndef ROCKABLE_ENABLE_SOFT_PARTICLES
  Logger::warn("TrackStress: Rockable was compiled without ROCKABLE_ENABLE_SOFT_PARTICLES, nothing will be recorded");
#endif

  // Documentation of the output file
  docString << "istart = " << istart << ", nrec = " << nrec;
  columnDoc.clear();
  columnDoc.push_back("Time");
  for (size_t i = istart; i < box->Particles.size(); ++i) {
    columnDoc.push_back("sigma_xx of particle " + std::to_string(i));
  }
}

void TrackStress::init() {}

void TrackStress::exec() {}

void TrackStress::record() {
#ifdef ROCKABLE_ENABLE_SOFT_PARTICLES
  recordFile << box->t;
  for (size_t i = istart; i < box->Particles.size(); ++i) {
    const auto& SPP = box->Particles[i].SPP;
    recordFile << ' ' << ((SPP != nullptr) ? SPP->stress.xx : 0.0);
  }
  recordFile << std::endl << std::flush;
#endif
}

void TrackStress::end() {}
