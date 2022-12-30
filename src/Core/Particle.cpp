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

#include "Particle.hpp"

// Ctor
Particle::Particle()
    : group(0),
      cluster(0),
      pos(),
      vel(),
      acc(),
      Q(),
      vrot(),
      arot(),
      shape(0),
      homothety(1.0),
      inertia(),
      mass(0.0),
      force(),
      moment(),
      obb() {}

/** @brief Given a vector expressed in the shape framework,
 *         this method returns its expression in the global framework
 *
 */
vec3r Particle::Glob(vec3r& Vertex) const { return (pos + Q * (homothety * Vertex)); }

/** @brief Given the vertex number v of the shape,
 *         this method returns the vector expression in the global framework
 *
 */
vec3r Particle::GlobVertex(size_t v) const { return (pos + Q * (homothety * shape->vertex[v])); }

/** @brief Given the vertex index v of the face f in the shape,
 *         this method returns the vector expression in the global framework
 *
 */
vec3r Particle::GlobFaceVertex(size_t f, size_t v) const {
  return (pos + Q * (homothety * shape->vertex[shape->face[f][v]]));
}

/**
 *  @brief Get the scaled Minskowski radius
 *
 */
double Particle::MinskowskiRadius() const { return homothety * shape->radius; }

/**
 * @brief Update the placement of the Oriented Bounding Box according to the current position and orientation of the
 *        particle
 *
 */
void Particle::updateObb() {
  obb = shape->obb;
  obb.rotate(Q);
  obb.extent *= homothety;
  obb.center *= homothety;
  obb.center += pos;
}

bool Particle::VertexIsNearVertex(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax,
                                  const vec3r& branchPerioCorr) {
  vec3r pos_iv = Pi.GlobVertex(isub);
  vec3r pos_jv = Pj.GlobVertex(jsub) + branchPerioCorr;
  double sum = Pi.MinskowskiRadius() + Pj.MinskowskiRadius() + dmax;
  return (norm2(pos_jv - pos_iv) <= sum * sum);
}

bool Particle::VertexIsNearEdge(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax,
                                const vec3r& branchPerioCorr) {
  vec3r pos_iv = Pi.GlobVertex(isub) - branchPerioCorr;
  size_t v1 = Pj.shape->edge[jsub].first;
  size_t v2 = Pj.shape->edge[jsub].second;
  vec3r pos0_jv = Pj.GlobVertex(v1);
  vec3r pos1_jv = Pj.GlobVertex(v2);

  vec3r E = pos1_jv - pos0_jv;
  vec3r v = pos_iv - pos0_jv;
  double r = (v * E) / (E * E);

  if (r < -dmax) return false;
  if (r > 1.0 + dmax) return false;

  double sum = Pi.MinskowskiRadius() + Pj.MinskowskiRadius() + dmax;

  return (norm2(pos_iv - (pos0_jv + r * E)) <= sum * sum);
}

bool Particle::VertexIsNearFace(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax,
                                const vec3r& branchPerioCorr) {
  // First, we project the node position onto the face plane.
  size_t nb_vertices = Pj.shape->face[jsub].size();
  vec3r posNodeA_jv = Pj.GlobFaceVertex(jsub, 0);
  vec3r posNodeB_jv = Pj.GlobFaceVertex(jsub, 1);
  vec3r posNodeC_jv = Pj.GlobFaceVertex(jsub, nb_vertices - 1);
  vec3r pos_iv = Pi.GlobVertex(isub) - branchPerioCorr;
  vec3r v = pos_iv - posNodeA_jv;
  vec3r v1 = posNodeB_jv - posNodeA_jv;
  v1.normalize();
  vec3r v2 = posNodeC_jv - posNodeA_jv;
  vec3r n = cross(v1, v2);
  n.normalize();
  double dist = v * n;
  if (dist < 0.0) {
    n = -n;
    dist = -dist;
  }  // normal to the face is NOT considered
  vec3r P = pos_iv - dist * n;

  // Then, we check whether the projected point is inside the face (a polygon).
  // We use the crossing number algorithm (also known as even-odd rule
  // algorithm)
  int ODD = 0;
  v2 = cross(n, v1);  // both n and v1 are normalized
  double ori1 = P * v1;
  double ori2 = P * v2;
  for (size_t iva = 0; iva < nb_vertices; ++iva) {
    size_t ivb = iva + 1;
    if (ivb == nb_vertices) ivb = 0;
    posNodeA_jv = Pj.GlobFaceVertex(jsub, iva);
    posNodeB_jv = Pj.GlobFaceVertex(jsub, ivb);
    double pa1 = posNodeA_jv * v1;
    double pb1 = posNodeB_jv * v1;
    double pa2 = posNodeA_jv * v2;
    double pb2 = posNodeB_jv * v2;

    // @see http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
    // @see http://alienryderflex.com/polygon/
    if ((pa2 < ori2 && pb2 >= ori2) || (pb2 < ori2 && pa2 >= ori2)) {
      if (pa1 + (ori2 - pa2) / (pb2 - pa2) * (pb1 - pa1) < ori1) {
        ODD = 1 - ODD;
      }
    }
  }

  if (ODD == 1) {
    dist -= (Pi.MinskowskiRadius() + Pj.MinskowskiRadius());
    if (dist <= dmax) return true;
  }
  return false;
}

bool Particle::EdgeIsNearEdge(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax,
                              const vec3r& branchPerioCorr) {
// Be carreful about this small value because, if it is not sufficiently small,
// some edges (tubes) couldn't see them each other.
#define _EPSILON_VALUE_ 1.0e-12
  size_t v1 = Pj.shape->edge[jsub].first;
  size_t v2 = Pj.shape->edge[jsub].second;
  vec3r pos0_jv = Pj.GlobVertex(v1);
  vec3r pos1_jv = Pj.GlobVertex(v2);
  v1 = Pi.shape->edge[isub].first;
  v2 = Pi.shape->edge[isub].second;
  vec3r pos0_iv = Pi.GlobVertex(v1);
  vec3r pos1_iv = Pi.GlobVertex(v2);

  vec3r Ei = pos1_iv - pos0_iv;
  vec3r Ej = pos1_jv - pos0_jv;
  vec3r v = pos0_iv - pos0_jv - branchPerioCorr;

  double c = Ei * Ei;
  double d = Ej * Ej;
  double e = Ei * Ej;
  double f = (c * d) - (e * e);
  vec3r n;

  if (fabs(f) > _EPSILON_VALUE_) {  // if not, the contact will necessary involve a vertex
    f = 1.0 / f;
    double a = Ei * v;
    double b = Ej * v;
    double s = (e * b - a * d) * f;  // for edge in shape i
    double t = (c * b - e * a) * f;  // for edge in shape j

    if (s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0) {
      n = (pos0_iv + s * Ei) - (pos0_jv + t * Ej);
      double sum = Pi.MinskowskiRadius() + Pj.MinskowskiRadius() + dmax;
      return (norm2(n) <= sum * sum);
    }
  }
  return false;
#undef _EPSILON_VALUE_
}
