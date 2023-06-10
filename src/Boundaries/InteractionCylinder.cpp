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

#include "InteractionCylinder.hpp"

InteractionCylinder::InteractionCylinder() : InteractionBoundary() {}
InteractionCylinder::InteractionCylinder(const InteractionCylinder& I) : InteractionBoundary(I) {}
InteractionCylinder::InteractionCylinder(size_t I, size_t Isub, double Damp) : InteractionBoundary(I,Isub,Damp) {}
//InteractionBoundary& InteractionCylinder::operator=(const InteractionCylinder& other) {
//  return InteractionBoundary::operator=(other);
//}

bool InteractionCylinder::update(Boundary* Bi, Particle& Pj) {
  Cylinder* Pi = static_cast<Cylinder*>(Bi);

  vec3r posi = Pi->pos;
  vec3r posj = Pj.GlobVertex(isub);
  vec3r pos_ij = posi - posj;  // from j to i
  double Ri = Pi->rCyl;
  double Rj = Pj.MinskowskiRadius();

  prev_n = n;
  vec3r N = pos_ij - dot(pos_ij, Pi->nCyl) * Pi->nCyl;
  double l = norm(N);
  if (l == 0.0) return false;
  N /= l;
  n = dot(pos_ij, N) * N - Ri * N;
  l = n.normalize();
  prev_dn = dn;
  dn = l - Rj;
  pos = posj - n * (Rj + 0.5 * dn);

  if (dn > 0.0) return false;
  vec3r PPi = pos - Pi->pos;  // PPi.round(periodicity);
  vec3r PPj = pos - Pj.pos;   // PPj.round(periodicity);
  // v(Qj) - v(Qi)
   
  vel = (Pj.vel - cross(PPj, Pj.vrot)) - ( -cross(PPi, Pi->vrot));

  return true;
}