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
