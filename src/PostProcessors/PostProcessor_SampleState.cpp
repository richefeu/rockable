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

#include <set>
#include <utility>

#include "factory.hpp"
#include "kwParser.hpp"

#include "Core/Rockable.hpp"
#include "PostProcessor_SampleState.hpp"

// Registered in Rockable::ExplicitRegistrations

SampleState::SampleState() {}

void SampleState::read(std::istream& is) {
  kwParser parser;
  parser.kwMap["WallBox"] = __DO__(istr) {
    wallBox.resize(6);
    for (size_t k = 0; k < 6; k++) istr >> wallBox[k];
  };
  parser.kwMap["Volume"] = __GET__(istr, Volume);
  parser.parse(is);

  // example:
  // firstConf 0
  // lastConf 50
  // stepConf 1
  // PostProcessor SampleState
  // WallBox 0 1 2 3 4 5   # the walls xmin xmax ymin ymax zmin zmax (not needed with a periodic cell)
}

void SampleState::init() {
  file.open("SampleState.txt");
  file << "# 1:t 2:Lx 3:Ly 4:Lz 5:V 6:solidFraction 7:Z 8:Ncontacts(pairs) 9:p 10:q 11:sxx 12:syy 13:szz 14:sxy 15:sxz "
          "16:syz\n";
  file << "# stress from the contact forces, positive in compression; p = trace/3, q = von Mises\n";
}

void SampleState::end() { file.close(); }

void SampleState::exec() {
  size_t nDriven = box->nDriven;
  bool periodic = false;
#ifdef ROCKABLE_ENABLE_PERIODIC
  periodic = (box->usePeriodicCell == 1);
#endif

  // ---- Volume and sizes of the sample
  double Lx = 0.0, Ly = 0.0, Lz = 0.0, V = Volume;
  if (wallBox.size() == 6) {
    // inner faces of the walls: position of the wall -/+ its radius (Minkowski radius times homothety)
    auto face = [&](size_t k, int comp, double side) {
      const Particle& W = box->Particles[wallBox[k]];
      return W.pos[comp] + side * W.shape->radius * W.homothety;
    };
    Lx = face(1, 0, -1.0) - face(0, 0, 1.0);
    Ly = face(3, 1, -1.0) - face(2, 1, 1.0);
    Lz = face(5, 2, -1.0) - face(4, 2, 1.0);
    if (Volume <= 0.0) V = Lx * Ly * Lz;
  }
#ifdef ROCKABLE_ENABLE_PERIODIC
  else if (periodic) {
    const mat9r& h = box->Cell.h;
    Lx = norm(vec3r(h.xx, h.yx, h.zx));
    Ly = norm(vec3r(h.xy, h.yy, h.zy));
    Lz = norm(vec3r(h.xz, h.yz, h.zz));
    if (Volume <= 0.0) V = fabs(h.det());
  }
#endif
  if (V <= 0.0) {
    std::cout << "SampleState: no volume (give WallBox or Volume, or use a periodic cell)" << std::endl;
    V = 1.0;
  }

  // ---- Volume of the free particles: the volume of their shape times the cube of their homothety
  double Vs = 0.0;
  size_t nFree = 0;
  for (size_t i = nDriven; i < box->Particles.size(); i++) {
    const Particle& P = box->Particles[i];
    Vs += P.shape->volume * P.homothety * P.homothety * P.homothety;
    nFree++;
  }

  // ---- Contacts and stress
  // sigma = (1/V) sum over the free particles p and their contacts c of F_p (x) (x_p - x_c),
  // where F_p is the force received by p; it is positive in compression.
  // The normal n goes from j to i, so that fn * n + ft is the force received by i.
  mat9r S;
  // Two polyhedra in contact often touch at several points (vertex-face, edge-edge...): the coordination
  // number counts the pairs of bodies in contact, not the points
  std::set<std::pair<size_t, size_t>> pairs;
  auto lever = [&](const vec3r& xp, const vec3r& xc) {
    vec3r d = xp - xc;
#ifdef ROCKABLE_ENABLE_PERIODIC
    if (periodic) {  // nearest image
      vec3r s = box->Cell.hinv * d;
      s.x -= std::round(s.x);
      s.y -= std::round(s.y);
      s.z -= std::round(s.z);
      d = box->Cell.h * s;
    }
#endif
    return d;
  };
  for (size_t k = 0; k < box->Interactions.size(); k++) {
    for (auto it = box->Interactions[k].begin(); it != box->Interactions[k].end(); ++it) {
      const Interaction& I = *it;
      if (I.dn >= 0.0 && I.stick == nullptr) continue;  // not a contact (nor a bond)
      vec3r f = I.fn * I.n + I.ft;                     // received by i
      if (norm2(f) == 0.0) continue;
      pairs.insert(std::make_pair(I.i, I.j));
      if (I.i >= nDriven) S += dyadic_product(f, lever(box->Particles[I.i].pos, I.pos));
      if (I.j >= nDriven) S += dyadic_product(-f, lever(box->Particles[I.j].pos, I.pos));
    }
  }
  S *= (1.0 / V);

  double p = (S.xx + S.yy + S.zz) / 3.0;
  double sxy = 0.5 * (S.xy + S.yx), sxz = 0.5 * (S.xz + S.zx), syz = 0.5 * (S.yz + S.zy);
  double q = sqrt(0.5 * ((S.xx - S.yy) * (S.xx - S.yy) + (S.yy - S.zz) * (S.yy - S.zz) + (S.zz - S.xx) * (S.zz - S.xx)) +
                  3.0 * (sxy * sxy + sxz * sxz + syz * syz));
  size_t nContacts = pairs.size();
  double nEnds = 0.0;  // number of contact ends on free bodies
  for (auto& pr : pairs) {
    if (pr.first >= nDriven) nEnds += 1.0;
    if (pr.second >= nDriven) nEnds += 1.0;
  }
  double Z = (nFree > 0) ? nEnds / (double)nFree : 0.0;

  file << box->t << ' ' << Lx << ' ' << Ly << ' ' << Lz << ' ' << V << ' ' << Vs / V << ' ' << Z << ' ' << nContacts
       << ' ' << p << ' ' << q << ' ' << S.xx << ' ' << S.yy << ' ' << S.zz << ' ' << sxy << ' ' << sxz << ' ' << syz
       << std::endl;
  std::cout << "  t = " << box->t << ", solid fraction = " << Vs / V << ", Z = " << Z << ", p = " << p << std::endl;
}
