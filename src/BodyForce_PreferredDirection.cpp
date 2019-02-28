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

#include "Rockable.hpp"
#include "BodyForce_PreferredDirection.hpp"

static Registrar<BodyForce, PreferredDirection> registrar("PreferredDirection");

PreferredDirection::PreferredDirection() { }

void PreferredDirection::read(std::istream& is) {
  is >> axisBody >> axis >> momentMax;
  axisBody.normalize();
  axis.normalize();
  Kr = 2.0 * momentMax / M_PI;
}

void PreferredDirection::write(std::ostream& os) {
  os << "PreferredDirection " << axisBody << ' ' << axis << ' ' << momentMax << '\n';
}

void PreferredDirection::getForceAndMoment(size_t ibody, vec3r & /*force*/, vec3r & moment) {

  vec3r axisB = box->Particles[ibody].Q * axisBody;
  axisB.normalize();
  vec3r u = cross(axis, axisB);
  u.normalize();
  double cosa = axisB * axis;
  double angle = acos(cosa);
  double sgn = -1.0;
  if (cosa < 0.0) {
    sgn = 1.0;
    angle = M_PI - angle;
  }
  
  //cout << "u  = " << u << '\n';
  //cout << "cos a = " << cosa << '\n';
  //cout << "angle = " << angle *180./M_PI << '\n';
  //cout << "sgn  = " << sgn << '\n';
  //cout << "sgn X u X angle = " << sgn*u*angle *180./M_PI << '\n';
  
  moment = Kr * sgn * u * angle;
}
