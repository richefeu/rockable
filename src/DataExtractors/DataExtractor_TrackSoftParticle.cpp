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
#include "DataExtractor_TrackSoftParticle.hpp"

// State of a soft particle: deformation gradient, stress and strain (body frame)
TrackSoftParticle::TrackSoftParticle() : ibody(0) {}

void TrackSoftParticle::read(std::istream& is) {
  is >> ibody >> filename >> nrec;
  if (box->isInteractive() == false) recordFile.open(filename.c_str());
  nstep = std::numeric_limits<int>::max();

#ifndef ROCKABLE_ENABLE_SOFT_PARTICLES
  Logger::warn(
      "TrackSoftParticle: Rockable was compiled without ROCKABLE_ENABLE_SOFT_PARTICLES, nothing will be recorded");
#endif

  // Documentation of the output file
  docString << "ibody = " << ibody << ", nrec = " << nrec;
  columnDoc.clear();
  columnDoc.push_back("Time");
  columnDoc.push_back("det(F) (volume ratio)");
  columnDoc.push_back("trace of the stress");
  columnDoc.push_back("trace of the strain (volumetric strain)");
  const char* comp[] = {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};
  for (const char* c : comp) columnDoc.push_back(std::string("F.") + c);
  for (const char* c : comp) columnDoc.push_back(std::string("stress.") + c);
  for (const char* c : comp) columnDoc.push_back(std::string("strain.") + c);
}

void TrackSoftParticle::init() {
#ifdef ROCKABLE_ENABLE_SOFT_PARTICLES
  if (ibody >= box->Particles.size() || box->Particles[ibody].SPP == nullptr) {
    Logger::warn("TrackSoftParticle: the particle {} is not a soft particle, nothing will be recorded", ibody);
  }
#endif
}

void TrackSoftParticle::exec() {}

void TrackSoftParticle::record() {
#ifdef ROCKABLE_ENABLE_SOFT_PARTICLES
  if (ibody >= box->Particles.size() || box->Particles[ibody].SPP == nullptr) {
    return;
  }
  const auto& SPP = box->Particles[ibody].SPP;
  recordFile << box->t << ' ' << SPP->F.det() << ' ' << SPP->stress.trace() << ' ' << SPP->strain.trace() << ' '
             << SPP->F << ' ' << SPP->stress << ' ' << SPP->strain << std::endl
             << std::flush;
#endif
}

void TrackSoftParticle::end() {}
