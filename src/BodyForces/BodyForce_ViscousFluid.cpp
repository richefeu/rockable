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

#include "BodyForce_ViscousFluid.hpp"
#include "Rockable.hpp"

//static Registrar<BodyForce, ViscousFluid> registrar("ViscousFluid");

ViscousFluid::ViscousFluid() {}

void ViscousFluid::read(std::istream& is) { is >> fluidDensity; }

void ViscousFluid::write(std::ostream& os) { os << "ViscousFluid " << fluidDensity << '\n'; }

/// @attention @emoji :construction: This solution has never been tested. So we don't know if the implementation is ok
/// @date November-2020
/// @author Vincent Richefeu
void ViscousFluid::getForceAndMoment(size_t ibody, vec3r& force, vec3r& /*moment*/) {
  static const double _2_DIV_3 = 2.0 / 3.0;
  static const double HALF_CX_SPHERE = 0.5 * 0.47;

  vec3r velSqr(box->Particles[ibody].vel.x * box->Particles[ibody].vel.x,
               box->Particles[ibody].vel.y * box->Particles[ibody].vel.y,
               box->Particles[ibody].vel.z * box->Particles[ibody].vel.z);
  double h = box->Particles[ibody].homothety;
  double pVolume = box->Particles[ibody].shape->volume * h * h * h;
  double pSurface = pow(pVolume, _2_DIV_3);
  force = HALF_CX_SPHERE * fluidDensity * pSurface * velSqr;
}
