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

#include <limits>

#include "factory.hpp"

#include "DataExtractor_TrackBody.hpp"
#include "Rockable.hpp"

static Registrar<DataExtractor, TrackBody> registrar("TrackBody");

TrackBody::TrackBody() : ibody(0), ictrl(0) {}

void TrackBody::read(std::istream& is) {
  is >> ibody >> filename >> nrec;
  if (box->isInteractive() == false) recordFile.open(filename.c_str());
  nstep = std::numeric_limits<int>::max();

  // Documentation of the output file
  docString << "ibody = " << ibody << ", nrec = " << nrec;
  columnDoc.clear();
  columnDoc.push_back("Time");
  columnDoc.push_back("x");
  columnDoc.push_back("y");
  columnDoc.push_back("z");
  columnDoc.push_back("v.x");
  columnDoc.push_back("v.y");
  columnDoc.push_back("v.z");
  columnDoc.push_back("Q.w");
  columnDoc.push_back("Q.x");
  columnDoc.push_back("Q.y");
  columnDoc.push_back("Q.z");
  columnDoc.push_back("vrot.x");
  columnDoc.push_back("vrot.y");
  columnDoc.push_back("vrot.z");
  columnDoc.push_back("f.x - imposedForce.x");
  columnDoc.push_back("f.y - imposedForce.y");
  columnDoc.push_back("f.z - imposedForce.z");
  columnDoc.push_back("m.x");
  columnDoc.push_back("m.y");
  columnDoc.push_back("m.z");
}

void TrackBody::init() {
  ictrl = 0;
  x_ForceImposed = y_ForceImposed = z_ForceImposed = false;
  for (size_t c = 0; c < box->System.controls.size(); ++c) {
    if (box->System.controls[c].i == ibody) {
      ictrl = c;
      switch (box->System.controls[c].type) {
        case _x_For_:
          x_ForceImposed = true;
          break;
        case _y_For_:
          y_ForceImposed = true;
          break;
        case _z_For_:
          z_ForceImposed = true;
          break;
        default:
          break;
      }
      break;
    }
  }
}

void TrackBody::exec() {}

void TrackBody::record() {
  vec3r ImposedForce;
  // By default the imposed force is null,
  // but if a component is with imposed force, then this force component 
  // will be substracted in the recordFile
  if (x_ForceImposed == true) ImposedForce.x = box->System.controls[ictrl].value;
  if (y_ForceImposed == true) ImposedForce.y = box->System.controls[ictrl].value;
  if (z_ForceImposed == true) ImposedForce.z = box->System.controls[ictrl].value;

  recordFile << box->t << ' ' << box->Particles[ibody].pos << ' ' << box->Particles[ibody].vel << ' '
             << box->Particles[ibody].Q << ' ' << box->Particles[ibody].vrot << ' '
             << box->Particles[ibody].force - ImposedForce << ' ' << box->Particles[ibody].moment << std::endl
             << std::flush;
}

void TrackBody::end() {}
