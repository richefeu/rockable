//        Rockable, 3D-DEM with sphero-polyhedra
//        Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>
//        
//        This program is free software: you can redistribute it and/or modify
//        it under the terms of the GNU General Public License as published by
//        the Free Software Foundation, either version 3 of the License, or
//        (at your option) any later version.
//        
//        This program is distributed in the hope that it will be useful,
//        but WITHOUT ANY WARRANTY; without even the implied warranty of
//        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//        GNU General Public License for more details.
//        
//        You should have received a copy of the GNU General Public License
//        along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "factory.hpp"

#include "Rockable.hpp"
#include "BodyForce_AttractingPoint.hpp"

static Registrar<BodyForce, AttractingPoint> registrar("AttractingPoint");

AttractingPoint::AttractingPoint() { }

void AttractingPoint::read(std::istream& is) {
  is >> point >> acceleration;
}

void AttractingPoint::write(std::ostream& os) {
  os << "AttractingPoint " << point << ' ' << acceleration << '\n';
}

void AttractingPoint::getForceAndMoment(size_t ibody, vec3r & force, vec3r & /*moment*/) {
  vec3r direction = point - box->Particles[ibody].pos;
  double l2 = norm2(direction);
  
  if (l2 > 0.0) {
    direction *= 1.0 / sqrt(l2); // normalize
    force = box->Particles[ibody].mass * acceleration * direction;
  } 
}
