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

#include "patatrac.hpp"

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
  //file >> distCol; // ???
  size_t nb = 0;
  file >> nb;
  BlockRelease B;
  for (size_t i = 0 ; i < nb ; i++) {
    file >> B.pos >> B.vel >> B.Q >> B.vrot;
    releases.push_back(B);
  }
}

// Given a freeFlight where pos, vel, Q and vrot have been initialized,
// the method computes the time tcol of intersection of the trajectory with
// the DTM (without thickness)
// The return-value can be
//   0  OK
//   1  the trajectory did not intersect the DTM (F.duration is set to 0.0)
//   2  F.duration is negative (should never happen)
int Patatrac::getCollisionTime(freeFlight& F) {
  std::cout << "> Compute raw flight duration..." << std::endl;
  F.duration = 1.0e20;  // it will be the smallest duration
  int id_face = -1;
  //int id_terrain = -1;
  
  for (size_t iterrain = 0 ; iterrain < box.nDriven ; iterrain++) {
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
      double h = 0.0; // flight duration up to next intersection with the infinit plan
      if (A == 0.0) { // gravity perpendicular to plan normal
        if (B == 0.0) { // velocity perpendicular to plan normal
          // C = 0 -> no intersection possible
          continue;
        }
        else {
          // affine relation Bh + C = 0
          if (C == 0) h = 0.0;
          else {
            h = fabs(-C / B);
          }
        }
      }
      else {
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
          F.duration = h;
          id_face = i;
          //id_terrain = iterrain;
        }

      }  // end if (h <  F.duration)
    }    // end-for terrain faces
  } // end-for terrains
  
  /*
  __SHOW(id_terrain);
  __SHOW(id_face);
  vec3r po = (F.pos + F.vel * F.duration + 0.5 * box.gravity * F.duration * F.duration);
  __SHOW(po);
  */
  
  if (id_face < 0) { // No face has been found
    F.duration = 0.0;
    return 1;
  }
  
  if (F.duration <= 0.0) {
    //F.vel.reset();
    F.duration = 0.0;
    return 2;
  }

  return 0;
}

// Given a freeFlight where pos, vel, Q and vrot have been initialized,
// the method computes the time tcol of intersection of the trajectory with
// the DTM (without thickness)
// The return-value can be
//   0  OK
//   1  the trajectory did not intersect the DTM (F.duration is set to 0.0)
//   2  F.duration is negative (should never happen)
int Patatrac::getCollisionTime_v2(freeFlight& F) {
  std::cout << "> Compute raw flight duration..." << std::endl;
  F.duration = 1.0e20;  // it will be the smallest duration
  int id_face = -1;
  Shape* terrain = box.Particles[0].shape;

  for (size_t i = 0; i < terrain->face.size(); i++) {

    size_t nb_vertices = terrain->face[i].size();
    vec3r posA = box.Particles[0].GlobFaceVertex(i, 0);
    vec3r posB = box.Particles[0].GlobFaceVertex(i, 1);
    vec3r posC = box.Particles[0].GlobFaceVertex(i, nb_vertices - 1);
    vec3r v1 = posB - posA;
    vec3r v2 = posC - posA;
    vec3r n = cross(v1, v2);
    n.normalize();

    // first normal direction
    double A = 0.5 * (box.gravity * n);
    double B = F.vel * n;
    double C = (F.pos - posA) * n - distCol;
    double DELTA = B * B - 4.0 * A * C;
    if (DELTA < 0.0) continue;
    double h1 = (-B + sqrt(DELTA)) / (2.0 * A);
    if (h1 < 0.0) h1 = (-B - sqrt(DELTA)) / (2.0 * A);

    // second normal direction
    A = -A;
    B = -B;
    C = -(F.pos - posA) * n - distCol;
    DELTA = B * B - 4.0 * A * C;
    if (DELTA < 0.0) continue;
    double h2 = (-B + sqrt(DELTA)) / (2.0 * A);
    if (h2 < 0.0) h2 = (-B - sqrt(DELTA)) / (2.0 * A);

    // get smallest positive time
    double h = h1;
    if (h2 < h1) h = h2;

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
      double pa1, pa2;
      double pb1, pb2;
      size_t iva, ivb;
      for (iva = 0; iva < nb_vertices; ++iva) {
        ivb = iva + 1;
        if (ivb == nb_vertices) ivb = 0;
        posA = box.Particles[0].GlobFaceVertex(i, iva);
        posB = box.Particles[0].GlobFaceVertex(i, ivb);
        pa1 = posA * v1;
        pb1 = posB * v1;
        pa2 = posA * v2;
        pb2 = posB * v2;

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
        F.duration = h;
        id_face = i;
      }

    }  // end if (h <  F.duration)
  }    // end-for terrain faces
  if (id_face < 0) {
    F.duration = 0.0;
    return 1;
  }
  if (F.duration <= 0.0) {
    F.vel.reset();
    return 2;
  }

  return 0;
}

// Get the minimum negative (but non-zero) normal distance
// Remark: this is where most of the time is spent
double Patatrac::getDnMin() {
  /*
  double dnMin = 1e20;
  if (box.activeInteractions.empty()) return dnMin;
  
  for (size_t k = 0 ; k < box.activeInteractions.size() ; k++) {
    Interaction* I = box.activeInteractions[k];//const_cast<Interaction*>(std::addressof(*it));
    Interaction::UpdateDispatcher[I->type](*I, box.Particles[I->i], box.Particles[I->j]);
    if (I->dn < dnMin && I->dn != 0.0) dnMin = I->dn;
  }

  return dnMin;
  */
  
  double dnMin = 1e20;
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
  std::cout << "> Refine flight duration..." << std::endl;
  
  double h = F.duration;
  box.t = F.ti + F.duration;
  
  box.Particles[iblock].pos = F.pos + F.vel * h + 0.5 * box.gravity * h * h;
  box.Particles[iblock].vel = F.vel + box.gravity * h;
  box.Particles[iblock].vrot = F.vrot;
  box.Particles[iblock].Q = F.Q;
  box.Particles[iblock].Q += ((box.Particles[iblock].Q.dot(F.vrot)) *= h);

  box.UpdateNL();
  double dn = getDnMin();

  double dt_2 = 0.5 * box.dt;
  double interVerletC = 0.0;
  
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
      interVerletC = 0.0;
    }
    box.t -= box.dt;
    interVerletC += box.dt;
  }

  // the "refined" duration
  F.duration = h;
}

/*
// Adjust the duration 
void Patatrac::refineCollisionTime(freeFlight& F) {
  std::cout << "> Refine flight duration..." << std::endl;
  
  epsilon_dn = 0.01 * box.Particles[iblock].MinskowskiRadius();
  // Range of duration:
  double tup = F.duration;
  double tdown = 0.0;
  __SHOW(tdown);
  __SHOW(tup);
  
  double t = 0.5 * tup;
  __SHOW(t);

  box.Particles[iblock].pos = F.pos + F.vel * t + 0.5 * box.gravity * t * t;
  box.Particles[iblock].vel = F.vel + box.gravity * t;
  box.Particles[iblock].vrot = F.vrot;
  box.Particles[iblock].Q = F.Q;
  box.Particles[iblock].Q += ((box.Particles[iblock].Q.dot(F.vrot)) *= t);

  __SHOW(box.Particles[iblock].pos);
  box.UpdateNL();
  double dn = getDnMin();
  __SHOW(dn);
  // std::cout << "dn = " << dn << std::endl;

  size_t attempt = 0;
  while ( attempt < 20) {

    __SHOW(attempt);
    __SHOW(dn);
    if (dn < 0.0)
      tup = t;
    else
      tdown = t;
    __SHOW(tdown);
    __SHOW(tup);
    
    t = tdown + 0.5 * (tup - tdown);
    
    __SHOW(t);
    // a new position
    attempt++;
    box.Particles[iblock].pos = F.pos + F.vel * t + 0.5 * box.gravity * t * t;
    box.Particles[iblock].vel = F.vel + box.gravity * t;
    box.Particles[iblock].vrot = F.vrot;
    box.Particles[iblock].Q = F.Q;
    box.Particles[iblock].Q += ((box.Particles[iblock].Q.dot(F.vrot)) *= t);

    // if (dn > epsilon_dn)
    __SHOW(box.Particles[iblock].pos);
    box.UpdateNL();
    dn = getDnMin();
    __SHOW(dn);
    // std::cout << "dn = " << dn << std::endl;
    if (dn > 0.0 && dn < epsilon_dn) break;
  }
  __SHOW(attempt);
  
  // the "refined" duration
  F.duration = t;
}
*/

// DEM collision starting at the end of the free flight.
// After the call of this function, the block will have
// a new position, velocity, orientation and angular velocity.
// Return value:
//   0  the particle will make another free flight
//   1  tmax has been reached
int Patatrac::computeCollision(freeFlight& F) {
  std::cout << "> Compute collision..." << std::endl;
  int nstepPacket = 600;
  
  //epsilon_dn = 0.01 * box.Particles[iblock].MinskowskiRadius();
  box.t = F.ti + F.duration;
  double tini = box.t;
  
  box.UpdateNL();
  double dn = getDnMin();
  epsilon_dn = dn;
  
  size_t nbPackets = 0;
  while (box.t < box.tmax) {
    // Packet of simuation steps
    nbPackets++;
    for (int i = 0 ; i < nstepPacket ; i++) {
      box.velocityVerletStep();
      box.t += box.dt;
    }  
    if (nbPackets >= 10) return 2;
    box.UpdateNL();
    dn = getDnMin();
    __SHOW(dn);
    if (dn > epsilon_dn) {
      std::cout << "  Collision duration: " << box.t - tini << std::endl;
      return 0;
    }
  }
  std::cout << "t = " << box.t << ", tmax = " << box.tmax << std::endl;
  return 1;
}

/*
// DEM collision starting at the end of the free flight.
// After the call of this function, the block will have
// a new position, velocity, orientation and angular velocity.
// Return value:
//   0  the particle will make another free flight
//   1  tmax has been reached
int Patatrac::computeCollision(freeFlight& F) {
  std::cout << "> Compute collision..." << std::endl;
  epsilon_dn = 0.01 * box.Particles[iblock].MinskowskiRadius();
  box.t = F.ti + F.duration;
  double tini = box.t;

  double dt_2 = 0.5 * box.dt;
  double interVerletC = 0.0;
  double dn = getDnMin();

  bool touch = false;
  while (box.t < box.tmax) {
    box.velocityVerletStep();
    if (interVerletC >= box.interVerlet - dt_2) {
      box.UpdateNL();
      interVerletC = 0.0;
    }
    box.t += box.dt;
    interVerletC += box.dt;
    
    dn = getDnMin();
    if (dn < 0.0) touch = true;
    if (touch && dn > epsilon_dn) {
      std::cout << "  Collision duration: " << box.t - tini << std::endl;
      return 0;
    }
  }
  std::cout << "t = " << box.t << ", tmax = " << box.tmax << std::endl;
  return 1;
}
*/

// Releases for all configurations
void Patatrac::dropAll() {
  //__SHOW(box.Particles[iblock].mass);
  __SHOW(iblock);
  __SHOW(3.14*sqrt(box.Particles[iblock].mass/(0.5e9)));
  __SHOW((3.14*sqrt(box.Particles[iblock].mass/(0.11e8))/box.dt));// [pi*srqt(m/(en2*kn))] / dt
  
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
      if (resCollision == 0) {
        refineCollisionTime(F);
        std::cout << "  refined flight duration: " << F.duration << std::endl;
        std::cout << "> Add free flight..." << std::endl;
        releases[r].freeFlights.push_back(F);
      } else {
        std::cout << "  !! SOMETHING WEIRD HAPPENED !!" << std::endl;
        box.t = F.ti;
        box.Particles[iblock].pos = F.pos;
        box.Particles[iblock].vel = F.vel;
        box.Particles[iblock].Q = F.Q;
        box.Particles[iblock].vrot = F.vrot;
      }

      // std::cout << "run collision..." << std::endl;
      stopped = computeCollision(F);
      if (stopped) std::cout << "> Stopped!" << std::endl;
    } while (stopped == 0);
  }  // Loop over the releases
}

void Patatrac::dropAll_v2() {
  iblock = box.nDriven;
  for (size_t r = 0; r < releases.size(); r++) {
    box.Particles[iblock].pos = releases[r].pos;
    box.Particles[iblock].vel = releases[r].vel;
    box.Particles[iblock].Q = releases[r].Q;
    box.Particles[iblock].vrot = releases[r].vrot;

    std::cout << "Release number " << r << ", from position: " << releases[r].pos << std::endl;

    box.t = 0.0;
    int resCollision = 0;
    int stopped = 0;
    do {
      freeFlight F;
      F.ti = box.t;
      F.pos = box.Particles[iblock].pos;
      F.vel = box.Particles[iblock].vel;
      F.Q = box.Particles[iblock].Q;
      F.vrot = box.Particles[iblock].vrot;
      resCollision = getCollisionTime_v2(F);
      std::cout << "  raw flight duration: " << F.duration << std::endl;
      if (resCollision == 0) {
        // refineCollisionTime(F);
        // std::cout << "  refined flight duration: " << F.duration << std::endl;
        // TODO placer à la fin du vol libre
        std::cout << "> Add free flight..." << std::endl;
        releases[r].freeFlights.push_back(F);
      } else {
        std::cout << "!! SOMETHING WEIRD HAPPENED !!" << std::endl;
        box.t = F.ti;
        box.Particles[iblock].pos = F.pos;
        box.Particles[iblock].vel = F.vel;
        box.Particles[iblock].Q = F.Q;
        box.Particles[iblock].vrot = F.vrot;
      }

      // std::cout << "run collision..." << std::endl;
      stopped = computeCollision(F);
      if (stopped) std::cout << "> Block seems stopped!" << std::endl;
    } while (stopped == 0);
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
  //Ptt.dropAll_v2();
  Ptt.saveTrajectories("traj.txt");

  return 0;
}
