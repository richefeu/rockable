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

#include "Interaction.hpp"

Interaction::Interaction()
    : i(0),
      j(0),
      type(0),
      isub(0),
      jsub(0),
      prev_n(),
      n(),
      dn(0.0),
      prev_dn(0.0),
      pos(),
      vel(),
      fn(0.0),
      ft(),
      mom(),
      damp(0.0),
      stick(nullptr) {}

Interaction::Interaction(const Interaction& I)
    : i(I.i),
      j(I.j),
      type(I.type),
      isub(I.isub),
      jsub(I.jsub),
      prev_n(I.prev_n),
      n(I.n),
      dn(I.dn),
      prev_dn(I.prev_dn),
      pos(I.pos),
      vel(I.vel),
      fn(I.fn),
      ft(I.ft),
      mom(I.mom),
      damp(I.damp),
      stick(I.stick) {}

Interaction::Interaction(size_t I, size_t J, int Type, size_t Isub, size_t Jsub, double Damp, BreakableInterface* Stick)
    : i(I),
      j(J),
      type(Type),
      isub(Isub),
      jsub(Jsub),
      prev_n(),
      n(),
      dn(0.0),
      prev_dn(0.0),
      pos(),
      vel(),
      fn(0.0),
      ft(),
      mom(),
      damp(Damp),
      stick(Stick) {}

Interaction& Interaction::operator=(const Interaction& other) {
  if (this != &other) {  // prevent self-assignment
    i = other.i;
    j = other.j;
    type = other.type;
    isub = other.isub;
    jsub = other.jsub;
    prev_n = other.prev_n;
    n = other.n;
    dn = other.dn;
    prev_dn = other.prev_dn;
    pos = other.pos;
    vel = other.vel;
    fn = other.fn;
    ft = other.ft;
    mom = other.mom;
    damp = other.damp;
    stick = other.stick;
  }
  return *this;
}

void Interaction::deactivate() {
  dn = prev_dn = 0.0;
  fn = 0.0;
  ft.reset();
  mom.reset();
}

// Dispatching (array of lambdas)
// ===============================================================================
//  Be careful not to change I.type, I.j and I.jsub in the following lambdas,
//  because they are stored in an ordered std::set.
//  They will be called that way:
//  Interaction::UpdateDispatcher[IT->type](const_cast<Interaction&>(*IT), Pi, Pj);
//  where IT is set<Interaction>::iterator, and thus *IT is a const Interaction.
//  This is imposed to avoid the breakage of the order in std::set.
//
//  The returned value says if the interaction is valid or not.
//  "valid" means "with a valid value of dn"
// ===============================================================================

std::function<bool(Interaction&, Particle&, Particle&)> Interaction::UpdateDispatcher[4]{
    // ------ UpdateVertexVertex
    [](Interaction& I, Particle& Pi, Particle& Pj) -> bool {
      vec3r posi = Pi.GlobVertex(I.isub);
      vec3r posj = Pj.GlobVertex(I.jsub);
      double Ri = Pi.MinskowskiRadius();
      double Rj = Pj.MinskowskiRadius();

      I.prev_n = I.n;
      I.n = posi - posj;  // from j to i
      I.prev_dn = I.dn;
      I.dn = I.n.normalize() - (Ri + Rj);
      // if (I.stick == false && I.dn > 0.0) return false;

      I.pos = posi - I.n * (Ri + 0.5 * I.dn);
      // v(Qj) - v(Qi)
      I.vel = (Pj.vel - cross(I.pos - Pj.pos, Pj.vrot)) - (Pi.vel - cross(I.pos - Pi.pos, Pi.vrot));

      return true;
    },
    // ------ UpdateVertexEdge
    [](Interaction& I, Particle& Pi, Particle& Pj) -> bool {
      size_t v1 = Pj.shape->edge[I.jsub].first;
      size_t v2 = Pj.shape->edge[I.jsub].second;
      vec3r posj1 = Pj.GlobVertex(v1);
      vec3r posj2 = Pj.GlobVertex(v2);
      vec3r posi = Pi.GlobVertex(I.isub);

      vec3r E = posj2 - posj1;
      vec3r v = posi - posj1;
      double r = (v * E) / (E * E);

      if (r <= 0.0 || r >= 1.0) {
        I.deactivate();
        return false;
      }
      
      double Ri = Pi.MinskowskiRadius();
      double Rj = Pj.MinskowskiRadius();

      I.prev_n = I.n;
      I.n = posi - (posj1 + r * E);  // from j to i
      I.prev_dn = I.dn;
      I.dn = I.n.normalize() - (Ri + Rj);
      I.pos = posi - I.n * (Ri + 0.5 * I.dn);
      // v(Qj) - v(Qi)
      I.vel = (Pj.vel - cross(I.pos - Pj.pos, Pj.vrot)) - (Pi.vel - cross(I.pos - Pi.pos, Pi.vrot));

      return true;
    },
    // ------ UpdateVertexFace
    [](Interaction& I, Particle& Pi, Particle& Pj) -> bool {
      // First, we project the node position onto the face plane.
      vec3r P;
      size_t nb_vertices = Pj.shape->face[I.jsub].size();
      vec3r posNodeA_jv = Pj.GlobFaceVertex(I.jsub, 0);
      vec3r posNodeB_jv = Pj.GlobFaceVertex(I.jsub, 1);
      vec3r posNodeC_jv = Pj.GlobFaceVertex(I.jsub, nb_vertices - 1);
      vec3r pos_iv = Pi.GlobVertex(I.isub);
      vec3r v = pos_iv - posNodeA_jv;
      vec3r v1 = posNodeB_jv - posNodeA_jv;
      v1.normalize();
      vec3r v2 = posNodeC_jv - posNodeA_jv;

      I.prev_n = I.n;
      I.n = cross(v1, v2);
      I.n.normalize();
      double dist = v * I.n;
      if (dist < 0.0) {
        I.n = -I.n;
        dist = -dist;
      }  // n goes from the plan (j) to the vertex (i)
      P = pos_iv - dist * I.n;

      // Then, we check whether the projected point is inside the face (a 3D
      // polygon). We use the crossing number algorithm (also known as
      // even-odd rule algorithm)
      size_t ODD = 0;
      v2 = cross(I.n, v1);
      double ori1 = P * v1;
      double ori2 = P * v2;
      double pa1, pa2;
      double pb1, pb2;
      size_t iva, ivb;
      for (iva = 0; iva < nb_vertices; ++iva) {
        ivb = iva + 1;
        if (ivb == nb_vertices) ivb = 0;
        posNodeA_jv = Pj.GlobFaceVertex(I.jsub, iva);
        posNodeB_jv = Pj.GlobFaceVertex(I.jsub, ivb);
        pa1 = posNodeA_jv * v1;
        pb1 = posNodeB_jv * v1;
        pa2 = posNodeA_jv * v2;
        pb2 = posNodeB_jv * v2;

        // @see http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
        // @see http://alienryderflex.com/polygon/
        if ((pa2 < ori2 && pb2 >= ori2) || (pb2 < ori2 && pa2 >= ori2)) {
          if (pa1 + (ori2 - pa2) / (pb2 - pa2) * (pb1 - pa1) < ori1) {
            ODD = 1 - ODD;
          }
        }
      }

      // ODD = 1 means that the projected point is inside the 3D-polygonal
      // face
      if (ODD) {
        double Ri = Pi.MinskowskiRadius();
        double Rj = Pj.MinskowskiRadius();
        I.prev_dn = I.dn;
        I.dn = dist - (Ri + Rj);
        I.pos = pos_iv - I.n * (Ri + 0.5 * I.dn);
        // v(Qj) - v(Qi)
        I.vel = (Pj.vel - cross(I.pos - Pj.pos, Pj.vrot)) - (Pi.vel - cross(I.pos - Pi.pos, Pi.vrot));

        return true;
      }

      I.deactivate();
      return false;
    },
    // ------ UpdateEdgeEdge
    [](Interaction& I, Particle& Pi, Particle& Pj) -> bool {
#define _EPSILON_VALUE_ 1.0e-12
      // Be carreful about this small value because, if it is not
      // sufficiently small, some edges (tubes) may not see them each other.

      vec3r posi1 = Pi.GlobVertex(Pi.shape->edge[I.isub].first);
      vec3r posi2 = Pi.GlobVertex(Pi.shape->edge[I.isub].second);
      vec3r posj1 = Pj.GlobVertex(Pj.shape->edge[I.jsub].first);
      vec3r posj2 = Pj.GlobVertex(Pj.shape->edge[I.jsub].second);

      vec3r Ei = posi2 - posi1;
      vec3r Ej = posj2 - posj1;
      vec3r v = posi1 - posj1;
      double c = Ei * Ei;
      double d = Ej * Ej;
      double e = Ei * Ej;
      double f = (c * d) - (e * e);
      double s, t;

      if (fabs(f) > _EPSILON_VALUE_) {
        f = 1.0 / f;
        double a = Ei * v;
        double b = Ej * v;
        s = (e * b - a * d) * f;  // for edge i
        t = (c * b - e * a) * f;  // for edge j

        if (s <= 0.0 || s >= 1.0 || t <= 0.0 || t >= 1.0) {
          I.deactivate();
          return false;
        }

        I.prev_n = I.n;
        I.n = (posi1 + s * Ei) - (posj1 + t * Ej);  // from j to i
      } else {                                      // f = 0 means that the edges are parallel
        I.deactivate();
        return false;
        // The contact will necessary be held by a vertex (sphere)
        // (i.e., vertex-vertex or vertex-edge)
      }

      // from here n, s and t have been computed
      double Ri = Pi.MinskowskiRadius();
      double Rj = Pj.MinskowskiRadius();
      I.prev_dn = I.dn;
      I.dn = I.n.normalize() - (Ri + Rj);
      I.pos = posi1 + s * Ei - I.n * (Ri + 0.5 * I.dn);
      // v(Qj) - v(Qi)
      I.vel = (Pj.vel - cross(I.pos - Pj.pos, Pj.vrot)) - (Pi.vel - cross(I.pos - Pi.pos, Pi.vrot));
      return true;
#undef _EPSILON_VALUE_
    }};  // End of dispatching array of lambdas
