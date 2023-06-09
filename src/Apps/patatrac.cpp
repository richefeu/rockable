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

#include "patatrac.hpp"

#define ENABLE_PROFILING
#include "profiler.hpp"

// Ctor
Patatrac::Patatrac() : epsilon_dn(1e-6) {}

// Read data from a file
void Patatrac::loadDropConfigsFromFile(const char* name) {
  if (!fileTool::fileExists(name)) {
    std::cerr << msg::warn() << "@patatrac::loadDropConfigsFromFile, file '" << name << "' not been found."
              << msg::normal() << std::endl
              << std::endl;
    return;
  }

  std::ifstream file(name);
  size_t nb = 0;
  file >> nb;
  BlockRelease B;
  for (size_t i = 0; i < nb; i++) {
    file >> B.pos >> B.vel >> B.Q >> B.vrot;
    releases.push_back(B);
  }
}

// Given a freeFlight where pos, vel, Q and vrot have been initialized,
// the method computes the time tcol of intersection of the trajectory with
// the DTM (without thickness)
// The return-value can be
//   0  OK
//   1  the trajectory did not intersect the DTM (F.duration is set so that the trajectory ends at box.tmax)
//   2  F.duration is negative (should never happen)
int Patatrac::getCollisionTime(freeFlight& F) {
  START_TIMER("getCollisionTime");
  std::cout << "> Compute raw flight duration..." << std::endl;
  F.duration = 1.0e20;  // it will be the smallest duration
  int id_face = -1;

  for (size_t iterrain = 0; iterrain < box.nDriven; iterrain++) {
    Shape* terrain = box.Particles[iterrain].shape;
    for (size_t i = 0; i < terrain->face.size(); i++) {

      size_t nb_vertices = terrain->face[i].size();
      // Get 3 distinct points of the face
      vec3r posA = box.Particles[iterrain].GlobFaceVertex(i, 0);
      vec3r posB = box.Particles[iterrain].GlobFaceVertex(i, 1);
      vec3r posC = box.Particles[iterrain].GlobFaceVertex(i, nb_vertices - 1);
      // compute a normal vector
      vec3r v1 = posB - posA;
      vec3r v2 = posC - posA;
      vec3r n = cross(v1, v2);
      n.normalize();

      // Trajectory equation: Ah^2 + Bh + C = 0
      // where h is the time elapsed since F.ti
      double A = 0.5 * (box.gravity * n);
      double B = F.vel * n;
      double C = (F.pos - posA) * n;
      double h = 0.0;    // flight duration up to next intersection with the infinit plan
      if (A == 0.0) {    // gravity perpendicular to plan normal
        if (B == 0.0) {  // velocity perpendicular to plan normal
          // C = 0 -> no intersection possible
          continue;
        } else {
          // affine relation Bh + C = 0
          if (C == 0)
            h = 0.0;
          else {
            h = fabs(-C / B);
          }
        }
      } else {
        // second order polynome Ah^2 + Bh + C = 0
        double DELTA = B * B - 4.0 * A * C;
        h = (-B + sqrt(DELTA)) / (2.0 * A);
        if (h < 0.0) h = (-B - sqrt(DELTA)) / (2.0 * A);
      }

      if (h < F.duration) {
        vec3r P = F.pos + F.vel * h + 0.5 * box.gravity * h * h;

        // Then, we check whether the point P is inside the face (a 3D
        // polygon). We use the crossing number algorithm (also known as
        // even-odd rule algorithm)
        size_t ODD = 0;
        v1.normalize();
        v2 = cross(n, v1);
        double ori1 = P * v1;
        double ori2 = P * v2;
        // size_t iva, ivb;
        for (size_t iva = 0; iva < nb_vertices; ++iva) {
          size_t ivb = iva + 1;
          if (ivb == nb_vertices) ivb = 0;
          posA = box.Particles[iterrain].GlobFaceVertex(i, iva);
          posB = box.Particles[iterrain].GlobFaceVertex(i, ivb);
          double pa1 = posA * v1;
          double pb1 = posB * v1;
          double pa2 = posA * v2;
          double pb2 = posB * v2;

          // @see http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
          // @see http://alienryderflex.com/polygon/
          if ((pa2 < ori2 && pb2 >= ori2) || (pb2 < ori2 && pa2 >= ori2)) {
            if (pa1 + (ori2 - pa2) / (pb2 - pa2) * (pb1 - pa1) < ori1) {
              ODD = 1 - ODD;
            }
          }
        }

        // ODD = 1 means that the projected point is inside the 3D-polygonal face
        if (ODD) {
          F.duration = fabs(h);
          id_face = i;
          // id_terrain = iterrain;
        }

      }  // end if (h <  F.duration)
    }    // end-for terrain faces
  }      // end-for terrains

  if (id_face < 0) {  // No face has been found
    F.duration = box.tmax - F.ti;
    return 1;
  }

  if (F.duration <= 0.0) {
    std::cerr << "!! Should never happen !! F.duration = " << F.duration << '\n';
    F.duration = 0.0;
    return 2;
  }

  return 0;
}

// Get the minimum negative (but non-zero) normal distance
// Remark: this is where most of the time is spent
double Patatrac::getDnMin() {
  START_TIMER("getDnMin");
  double dnMin = NAN;
  if (box.Interactions.empty()) return dnMin;

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    for (auto it = box.Interactions[k].begin(); it != box.Interactions[k].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      bool confidence = Interaction::UpdateDispatcher[it->type](*I, box.Particles[it->i], box.Particles[it->j]);
      if (confidence && I->dn < dnMin && I->dn != 0.0) dnMin = I->dn;
    }
  }

  return dnMin;
}

// Adjust the duration
void Patatrac::refineCollisionTime(freeFlight& F) {
  START_TIMER("refineCollisionTime");
  std::cout << "> Refine flight duration..." << std::endl;

  /*
  double h = F.duration;
  box.t = F.ti + h;
  */

  double h = F.duration;
  box.Particles[iblock].pos = F.pos + F.vel * h + 0.5 * box.gravity * h * h;
  box.Particles[iblock].vel = F.vel + box.gravity * h;
  box.Particles[iblock].vrot = F.vrot;
  box.Particles[iblock].Q = F.Q;
  box.Particles[iblock].Q += ((box.Particles[iblock].Q.dot(F.vrot)) *= h);

  box.UpdateNL();

  h *= 0.98;

  for (int k = 0; k < 50; k++) {

    box.Particles[iblock].pos = F.pos + F.vel * h + 0.5 * box.gravity * h * h;
    box.Particles[iblock].vel = F.vel + box.gravity * h;
    box.Particles[iblock].vrot = F.vrot;
    box.Particles[iblock].Q = F.Q;
    box.Particles[iblock].Q += ((box.Particles[iblock].Q.dot(F.vrot)) *= h);

    box.UpdateNL();
    box.accelerations();
    vec3r accf = box.Particles[iblock].acc - box.gravity;
    if (accf.isnull()) {
      break;
    }

    h *= 0.98;
  }

  box.t = F.ti + h;

  /*
  double dt_2 = 0.5 * box.dt;
  double interVerletC = 0.0;

  box.UpdateNL();
  double dn = getDnMin();
  if (std::isnan(dn)) dn = 0.0;

  while (h > 0.0 && dn <= 0.0) {
    h -= box.dt;
    box.t -= box.dt;

    box.Particles[iblock].pos = F.pos + F.vel * h + 0.5 * box.gravity * h * h;
    box.Particles[iblock].vel = F.vel + box.gravity * h;
    box.Particles[iblock].vrot = F.vrot;
    box.Particles[iblock].Q = F.Q;
    box.Particles[iblock].Q += ((box.Particles[iblock].Q.dot(F.vrot)) *= h);

    if (interVerletC >= box.interVerlet - dt_2) {
      box.UpdateNL();
      dn = getDnMin();
      if (std::isnan(dn)) dn = 0.0;
      __SHOW(dn);
      interVerletC = 0.0;
    }
    box.t -= box.dt;
    interVerletC += box.dt;
  }
  */

  // the adjusted duration
  F.duration = h;
}

// DEM collision starting at the end of the free flight.
// After the call of this function, the block will have
// a new position, velocity, orientation and angular velocity.
// Return value:
//   0  don't stop the trajectory, the particle will make another free flight
//   1  stop the trajectory, tmax has been reached
//   2  stop the trajectory, kinetic energy too low
//   3  stop the trajectory, ...(TODO)
int Patatrac::computeCollision(freeFlight& F) {
  START_TIMER("computeCollision");
  std::cout << "> Compute collision..." << std::endl;

  box.t = F.ti + F.duration;
  double tini = box.t;

  box.UpdateNL();

  bool firstTouch = false;

  box.accelerations();
  if (!(box.Particles[iblock].acc - box.gravity).isnull(1e-6)) {
    firstTouch = true;
    std::cout << "Block touch at the very beginning of collision!!\n";
  }

  size_t nbPackets = 0;
  size_t nsteps_after_ejection = 0;
  size_t nsteps_kin_stop = 0;
  while (box.t < box.tmax) {
    nbPackets++;

    box.t += box.dt;
    box.velocityVerletStep();
    vec3r res = (box.Particles[iblock].acc - box.gravity);

    if (nbPackets >= 5000) {
      box.UpdateNL();
      nbPackets = 0;
    }

    if (norm2(box.Particles[iblock].vel) < 0.01) {
      nsteps_kin_stop++;
    } else {
      nsteps_kin_stop = 0;
    }

    if (firstTouch == false && !res.isnull(1e-6)) {
      F.duration = box.t - tini;
      firstTouch = true;
      continue;
    }

    if (firstTouch == true) {
      if (res.isnull(1e-6)) {
        nsteps_after_ejection++;
      } else {
        nsteps_after_ejection = 0;
      }
    }

    if (nsteps_after_ejection >= 5000) {
      std::cout << "  Collision duration: " << box.t - tini << std::endl;
      std::cout << "  from " << tini << " to " << box.t << std::endl;
      return 0;
    }

    if (nsteps_kin_stop >= 1000) {
      std::cout << "  Collision duration: " << box.t - tini << std::endl;
      std::cout << "  from " << tini << " to " << box.t << std::endl;
      return 2;
    }
  }
  std::cout << "  from " << tini << " to " << box.t << std::endl;
  std::cout << "  tmax = " << box.tmax << std::endl;
  return 1;
}

// Releases for all configurations
void Patatrac::dropAll() {
  START_TIMER("dropAll");
  __SHOW(iblock);
  __SHOW(3.14 * sqrt(box.Particles[iblock].mass / (0.5e9)));
  __SHOW((3.14 * sqrt(box.Particles[iblock].mass / (0.11e8)) / box.dt));  // [pi*srqt(m/(en2*kn))] / dt

  for (size_t r = 0; r < releases.size(); r++) {
    box.Particles[iblock].pos = releases[r].pos;
    box.Particles[iblock].vel = releases[r].vel;
    box.Particles[iblock].Q = releases[r].Q;
    box.Particles[iblock].vrot = releases[r].vrot;

    std::cout << "Release number " << r << ", from position: " << releases[r].pos << std::endl;

    box.t = 0.0;
    int stopped = 0;
    do {
      freeFlight F;
      F.ti = box.t;
      F.pos = box.Particles[iblock].pos;
      F.vel = box.Particles[iblock].vel;
      F.Q = box.Particles[iblock].Q;
      F.vrot = box.Particles[iblock].vrot;

      int resCollision = getCollisionTime(F);

      std::cout << "  raw flight duration: " << F.duration << std::endl;
      box.t = F.ti + F.duration;
      std::cout << "  from " << F.ti << " to " << box.t << std::endl;

      if (resCollision == 0) {
        refineCollisionTime(F);
        std::cout << "  refined flight duration: " << F.duration << std::endl;
        std::cout << "  from " << F.ti << " to " << box.t << std::endl;
        stopped = computeCollision(F);
        std::cout << "> Add free flight..." << std::endl;
        releases[r].freeFlights.push_back(F);
      } else {
        std::cout << " !! It seems the block went out of the terrain !!" << std::endl;
        std::cout << "> Add free flight..." << std::endl;
        releases[r].freeFlights.push_back(F);
        stopped = 10;
        break;
      }
    } while (stopped == 0);

    if (stopped == 1) {
      std::cout << "> Stopped! (tmax reached)" << std::endl;
    } else if (stopped == 2) {
      std::cout << "> Stopped! (stop moving)" << std::endl;
    } else if (stopped == 10) {
      std::cout << "> Stopped! (out of terrain)" << std::endl;
    }

  }  // Loop over the releases
}

void Patatrac::saveTrajectories(const char* name) {
  std::ofstream file(name);
  file << releases.size() << std::endl;
  for (size_t r = 0; r < releases.size(); r++) {
    file << releases[r].freeFlights.size() << std::endl;
    for (size_t f = 0; f < releases[r].freeFlights.size(); f++) {
      file << releases[r].freeFlights[f].ti << ' ' << releases[r].freeFlights[f].duration << ' '
           << releases[r].freeFlights[f].pos << ' ' << releases[r].freeFlights[f].vel << ' '
           << releases[r].freeFlights[f].Q << ' ' << releases[r].freeFlights[f].vrot << std::endl;
    }
    file << std::endl;
  }
}

int main(int argc, char const* argv[]) {
	RockableProfiler::ProfilerManager prof;

  if (argc != 3) {
    std::cout << std::endl;
    std::cout << "Usage: " << argv[0] << " conf-file drop-configs" << std::endl;
    std::cout << std::endl;
    return 0;
  }

  Patatrac Ptt;
  Ptt.box.loadConf(argv[1]);
  Ptt.iblock = Ptt.box.nDriven;
  Ptt.loadDropConfigsFromFile(argv[2]);

  Ptt.dropAll();
  Ptt.saveTrajectories("traj.txt");

  return 0;
}
