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

#include "Core/Rockable.hpp"
#include "PreproCommand_randomlyOrientedVelocitiesClusters.hpp"
#include "ProcessingTools/processingTool_getClusters.hpp"

randomlyOrientedVelocitiesClusters::randomlyOrientedVelocitiesClusters() { }

void randomlyOrientedVelocitiesClusters::addCommand() {  
  box->parser.kwMap["randomlyOrientedVelocitiesClusters"] = [this](std::istream& conf) {
    conf >> this->velocityMagnitude >> this->opt;
    exec();
  };
}

/**
   @brief Set random orientation to the clusters
   @param[in]  velocityMagnitude  Magnitude of all velocities (only orientations change)
   @param[in]  opt                An option. if opt = 1 then all the velocity vectors
                                  will be oriented towards negative y (downward)
*/
void randomlyOrientedVelocitiesClusters::exec() {
  std::vector<clusterParticles> clusters;
  getClusters(box, clusters);

  quat q;
  q.randomize(true);
  vec3r u(velocityMagnitude, 0.0, 0.0);
  for (size_t c = 0; c < clusters.size(); c++) {
    q.randomize();
    vec3r v = q * u;
    if (opt == 1) v.y = -fabs(v.y);
    for (size_t i = 0; i < clusters[c].particleId.size(); i++) {
      box->Particles[clusters[c].particleId[i]].vel = v;
    }
  }
}
