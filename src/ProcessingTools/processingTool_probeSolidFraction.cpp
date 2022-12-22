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

#include "processingTool_probeSolidFraction.hpp"

/**
 * @brief xxx
 *
 * @param aabb    The probe as an Axis Aligned Bounding Box
 * @param MCnstep Number of Monte Carlo steps
 * @return double Estimate of the solid fraction inside the probe
 */
double probeSolidFraction(Rockable *box, AABB& aabb, size_t MCnstep) {
  if (MCnstep == 0) return -1.0;

  // select the concerned particles
  OBB zone;
  zone.center = 0.5 * (aabb.min + aabb.max);
  zone.extent.set(0.5 * (aabb.max.x - aabb.min.x), 0.5 * (aabb.max.y - aabb.min.y), 0.5 * (aabb.max.z - aabb.min.z));
  std::vector<size_t> pid;
  for (size_t i = 0; i < box->Particles.size(); ++i) {
    box->Particles[i].updateObb();
    if (zone.intersect(box->Particles[i].obb)) {
      pid.push_back(i);
    }
  }

  vec3r pt3;
  std::vector<double> vv(3);
  Mth::sobolSequence(-3, vv);  // Initialize the Sobol sequence
  size_t count = 0;
  for (size_t imc = 0; imc < MCnstep; ++imc) {
    Mth::sobolSequence(3, vv);
    pt3.set(aabb.min.x + vv[0] * (aabb.max.x - aabb.min.x), aabb.min.y + vv[1] * (aabb.max.y - aabb.min.y),
            aabb.min.z + vv[2] * (aabb.max.z - aabb.min.z));

    bool inSolid = false;
    for (size_t ii = 0; ii < pid.size(); ii++) {
      size_t i = pid[ii];
      vec3r ptTest = pt3 - box->Particles[i].pos;
      quat Qinv = box->Particles[i].Q.get_conjugated();
      ptTest = Qinv * ptTest;
      ptTest /= box->Particles[i].homothety;

      if (box->Particles[i].shape->inside(ptTest)) {
        inSolid = true;
        break;
      }
    }
    if (inSolid) count++;
  }

  return ((double)count / (double)MCnstep);
}