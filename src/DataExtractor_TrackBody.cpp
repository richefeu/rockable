// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

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
