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

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

#include "ShapeSDF.hpp"
#include "geoTool.hpp"

namespace {

/// Squared distance from x to the segment [a, b], and the closest point.
/// The returned parameter r is clamped in [0, 1]
double sqDistToSegment(const vec3r& x, const vec3r& a, const vec3r& b, vec3r& closest, double& r) {
  vec3r E = b - a;
  double EE = E * E;
  if (EE < 1.0e-30) {  // degenerate segment
    r = 0.0;
    closest = a;
    return norm2(x - a);
  }
  r = ((x - a) * E) / EE;
  if (r < 0.0) {
    r = 0.0;
  } else if (r > 1.0) {
    r = 1.0;
  }
  closest = a + r * E;
  return norm2(x - closest);
}

/// Newell normal of a polygon. Robust for non-planar and non-convex polygons,
/// and its norm is twice the projected area, so a near-zero norm flags a
/// degenerate face
vec3r newellNormal(const std::vector<vec3r>& vertex, const std::vector<size_t>& poly) {
  vec3r n;
  size_t n_v = poly.size();
  for (size_t i = 0; i < n_v; ++i) {
    const vec3r& c = vertex[poly[i]];
    const vec3r& d = vertex[poly[(i + 1) % n_v]];
    n.x += (c.y - d.y) * (c.z + d.z);
    n.y += (c.z - d.z) * (c.x + d.x);
    n.z += (c.x - d.x) * (c.y + d.y);
  }
  return n;
}

std::pair<size_t, size_t> sortedPair(size_t a, size_t b) {
  return (a < b) ? std::pair<size_t, size_t>(a, b) : std::pair<size_t, size_t>(b, a);
}

}  // namespace

ShapeSDF::ShapeSDF(const Shape& s) : shp(&s) {
  buildAdjacency();
  orientFacesOutward();
  buildBVH();
}

// Collect, for each face-boundary segment, the faces that use it, and for each
// vertex, the faces incident to it. Faces with less than 3 vertices, or with a
// vanishing Newell normal, are degenerate and take no part in any of this
void ShapeSDF::buildAdjacency() {
  size_t nf = shp->face.size();
  vertexFaces.assign(shp->vertex.size(), std::vector<size_t>());
  orientedFace.clear();
  faceNormal.assign(nf, vec3r());
  edgeFaces.clear();

  orientedFace = shp->face;

  bool allFacesUsable = (nf > 0);
  for (size_t f = 0; f < nf; ++f) {
    const std::vector<size_t>& poly = orientedFace[f];
    if (poly.size() < 3) {
      allFacesUsable = false;
      continue;
    }
    vec3r n = newellNormal(shp->vertex, poly);
    if (norm2(n) < 1.0e-30) {  // zero-area face
      allFacesUsable = false;
      continue;
    }
    n.normalize();
    faceNormal[f] = n;

    for (size_t i = 0; i < poly.size(); ++i) {
      edgeFaces[sortedPair(poly[i], poly[(i + 1) % poly.size()])].push_back(f);
      vertexFaces[poly[i]].push_back(f);
    }
  }

  // A clean closed volume needs every face-boundary segment to be shared by
  // exactly two faces. Then the exact angle-weighted pseudonormal gives the sign
  closedVolume = allFacesUsable && !shp->isSurface;
  if (closedVolume) {
    for (std::map<std::pair<size_t, size_t>, std::vector<size_t> >::const_iterator it = edgeFaces.begin();
         it != edgeFaces.end(); ++it) {
      if (it->second.size() != 2) {
        closedVolume = false;
        break;
      }
    }
  }

  // A face soup that is not a clean manifold (an imperfect STL, say) can still
  // be a solid: fall back to ray casting for the sign, exactly as Shape::inside
  // does. But this is only meaningful when the soup is essentially closed. We
  // measure that by the fraction of boundary edges (edges used by a single
  // face): a nearly closed mesh like a slightly broken STL has almost none,
  // while a genuine open sheet (a plate) is nearly all boundary. An open sheet
  // is left with an undefined sign, i.e. the unsigned offset of its skeleton,
  // which is exactly its (correct, closed) slab of thickness 2R
  rayCastSign = false;
  if (!closedVolume && !shp->isSurface && !shp->face.empty() && !edgeFaces.empty()) {
    size_t nBoundary = 0;
    for (std::map<std::pair<size_t, size_t>, std::vector<size_t> >::const_iterator it = edgeFaces.begin();
         it != edgeFaces.end(); ++it) {
      if (it->second.size() == 1) ++nBoundary;
    }
    double boundaryFraction = (double)nBoundary / (double)edgeFaces.size();
    rayCastSign = (boundaryFraction < 0.10);
  }
}

// Point-in-polyhedron by ray casting, replicating the polyhedron block of
// Shape::inside(): a ray along +x is odd-crossing when the point is inside. The
// face BVH prunes to the faces the ray could actually meet (box ahead in x, and
// straddling the ray in y and z)
bool ShapeSDF::rayCastInside(const vec3r& x) const {
  if (faceNodes.empty()) return false;

  size_t nb_intersect = 0;
  int stack[64];
  int sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const BVHNode& node = faceNodes[stack[--sp]];
    // The +x ray misses this box if the box is entirely behind, or off the ray
    // in y or z
    if (node.box.max.x < x.x) continue;
    if (x.y < node.box.min.y || x.y > node.box.max.y) continue;
    if (x.z < node.box.min.z || x.z > node.box.max.z) continue;

    if (node.count > 0) {
      for (int i = 0; i < node.count; ++i) {
        size_t f = prims[faceIds[node.start + i]].index;
        const std::vector<size_t>& poly = orientedFace[f];
        if (poly.size() < 3) continue;
        const vec3r& v1 = shp->vertex[poly[0]];
        for (size_t v = 1; v + 1 < poly.size(); ++v) {
          const vec3r& v2 = shp->vertex[poly[v]];
          const vec3r& v3 = shp->vertex[poly[v + 1]];
          if (geoTool::intersectTriangle(x, vec3r::unit_x(), v1, v2, v3) > 0) {
            ++nb_intersect;
            break;
          }
        }
      }
    } else {
      stack[sp++] = node.left;
      stack[sp++] = node.right;
    }
  }
  return (nb_intersect % 2) != 0;
}

// Walk the face adjacency graph and flip the faces that disagree with their
// neighbour on the traversal direction of the shared segment. The .shp format
// does not guarantee any orientation, and Shape::inside() does not rely on one,
// so this has to be established here. The global sign is then fixed by the
// signed volume, which must come out positive for outward normals
void ShapeSDF::orientFacesOutward() {
  skelVolume = 0.0;
  if (!closedVolume) {
    return;
  }

  size_t nf = orientedFace.size();
  std::vector<char> visited(nf, 0);
  std::queue<size_t> todo;
  todo.push(0);
  visited[0] = 1;

  while (!todo.empty()) {
    size_t f = todo.front();
    todo.pop();
    const std::vector<size_t> poly = orientedFace[f];

    for (size_t i = 0; i < poly.size(); ++i) {
      size_t a = poly[i];
      size_t b = poly[(i + 1) % poly.size()];
      const std::vector<size_t>& neighbours = edgeFaces[sortedPair(a, b)];

      for (size_t k = 0; k < neighbours.size(); ++k) {
        size_t g = neighbours[k];
        if (g == f || visited[g]) {
          continue;
        }

        // f traverses the shared segment as (a, b). A coherently oriented
        // neighbour must traverse it as (b, a)
        std::vector<size_t>& other = orientedFace[g];
        bool sameDirection = false;
        for (size_t j = 0; j < other.size(); ++j) {
          if (other[j] == a && other[(j + 1) % other.size()] == b) {
            sameDirection = true;
            break;
          }
        }
        if (sameDirection) {
          std::reverse(other.begin(), other.end());
          faceNormal[g] = -faceNormal[g];
        }

        visited[g] = 1;
        todo.push(g);
      }
    }
  }

  // A face set can be closed and manifold yet built of several connected
  // components; the traversal above only orients the one holding face 0
  for (size_t f = 0; f < nf; ++f) {
    if (!visited[f]) {
      closedVolume = false;
      return;
    }
  }

  // Signed volume by the divergence theorem. The fan triangulation is exact
  // even for a non-convex planar polygon, the spurious triangles cancel out
  for (size_t f = 0; f < nf; ++f) {
    const std::vector<size_t>& poly = orientedFace[f];
    const vec3r& v0 = shp->vertex[poly[0]];
    for (size_t i = 1; i + 1 < poly.size(); ++i) {
      const vec3r& v1 = shp->vertex[poly[i]];
      const vec3r& v2 = shp->vertex[poly[i + 1]];
      skelVolume += (v0 * cross(v1, v2)) / 6.0;
    }
  }

  if (skelVolume < 0.0) {  // the whole thing was oriented inward
    for (size_t f = 0; f < nf; ++f) {
      std::reverse(orientedFace[f].begin(), orientedFace[f].end());
      faceNormal[f] = -faceNormal[f];
    }
    skelVolume = -skelVolume;
  }
}

double ShapeSDF::sqDistToBox(const vec3r& x, const AABB& b) {
  double d2 = 0.0;
  if (x.x < b.min.x) d2 += (b.min.x - x.x) * (b.min.x - x.x);
  else if (x.x > b.max.x) d2 += (x.x - b.max.x) * (x.x - b.max.x);
  if (x.y < b.min.y) d2 += (b.min.y - x.y) * (b.min.y - x.y);
  else if (x.y > b.max.y) d2 += (x.y - b.max.y) * (x.y - b.max.y);
  if (x.z < b.min.z) d2 += (b.min.z - x.z) * (b.min.z - x.z);
  else if (x.z > b.max.z) d2 += (x.z - b.max.z) * (x.z - b.max.z);
  return d2;
}

// Recursive median split on the longest axis of the centroid spread. Leaves
// hold up to a small number of primitives
int ShapeSDF::buildNode(std::vector<int>& ids, std::vector<BVHNode>& nodes, int start, int count) const {
  BVHNode node;
  node.box = primBox[ids[start]];
  for (int i = 1; i < count; ++i) node.box.merge(primBox[ids[start + i]]);

  const int LEAF = 4;
  if (count <= LEAF) {
    node.start = start;
    node.count = count;
    nodes.push_back(node);
    return (int)nodes.size() - 1;
  }

  // Split along the axis where the centroids spread most
  vec3r lo = primCentroid[ids[start]], hi = lo;
  for (int i = 1; i < count; ++i) {
    const vec3r& c = primCentroid[ids[start + i]];
    lo = component_min(lo, c);
    hi = component_max(hi, c);
  }
  vec3r ext = hi - lo;
  int axis = (ext.x >= ext.y && ext.x >= ext.z) ? 0 : ((ext.y >= ext.z) ? 1 : 2);

  int mid = start + count / 2;
  std::nth_element(ids.begin() + start, ids.begin() + mid, ids.begin() + start + count,
                   [this, axis](int a, int b) {
                     const vec3r& ca = primCentroid[a];
                     const vec3r& cb = primCentroid[b];
                     return (axis == 0) ? (ca.x < cb.x) : ((axis == 1) ? (ca.y < cb.y) : (ca.z < cb.z));
                   });

  int myIndex = (int)nodes.size();
  nodes.push_back(node);  // reserve slot; children are appended after
  int l = buildNode(ids, nodes, start, mid - start);
  int r = buildNode(ids, nodes, mid, start + count - mid);
  nodes[myIndex].left = l;
  nodes[myIndex].right = r;
  nodes[myIndex].count = 0;
  return myIndex;
}

// Build one BVH over all primitives (for the closest-point query) and one over
// the faces only (for the sign). Brute force is fine for a few primitives, but
// a detailed STL has thousands, evaluated on millions of grid nodes
void ShapeSDF::buildBVH() {
  prims.clear();
  primBox.clear();
  primCentroid.clear();
  bvhNodes.clear();
  bvhIds.clear();
  faceNodes.clear();
  faceIds.clear();

  for (size_t i = 0; i < shp->vertex.size(); ++i) {
    prims.push_back({SDF_PRIM_VERTEX, i});
    primBox.push_back(AABB(shp->vertex[i]));
    primCentroid.push_back(shp->vertex[i]);
  }
  for (size_t e = 0; e < shp->edge.size(); ++e) {
    prims.push_back({SDF_PRIM_EDGE, e});
    AABB b(shp->vertex[shp->edge[e].first], shp->vertex[shp->edge[e].second]);
    primBox.push_back(b);
    primCentroid.push_back(0.5 * (shp->vertex[shp->edge[e].first] + shp->vertex[shp->edge[e].second]));
  }
  std::vector<int> faceOnly;
  for (size_t f = 0; f < orientedFace.size(); ++f) {
    const std::vector<size_t>& poly = orientedFace[f];
    if (poly.size() < 3) continue;
    AABB b(shp->vertex[poly[0]]);
    vec3r c;
    for (size_t k = 0; k < poly.size(); ++k) {
      b.add(shp->vertex[poly[k]]);
      c += shp->vertex[poly[k]];
    }
    c *= (1.0 / (double)poly.size());
    int primId = (int)prims.size();
    prims.push_back({SDF_PRIM_FACE, f});
    primBox.push_back(b);
    primCentroid.push_back(c);
    faceOnly.push_back(primId);
  }

  if (!prims.empty()) {
    bvhIds.resize(prims.size());
    for (size_t i = 0; i < prims.size(); ++i) bvhIds[i] = (int)i;
    buildNode(bvhIds, bvhNodes, 0, (int)bvhIds.size());
  }
  if (!faceOnly.empty()) {
    faceIds = faceOnly;
    buildNode(faceIds, faceNodes, 0, (int)faceIds.size());
  }
}

// Closest point on a polygonal face. The point is projected onto the face
// plane; if the projection falls inside the polygon the closest point is the
// projection itself, otherwise it lies on the boundary
ShapeSDF::FaceClosest ShapeSDF::closestOnFace(const vec3r& x, size_t f) const {
  FaceClosest result;
  const std::vector<size_t>& poly = orientedFace[f];
  size_t n_v = poly.size();

  if (n_v < 3 || faceNormal[f].isnull()) {  // degenerate face, push it out of the running
    result.d2 = std::numeric_limits<double>::max();
    return result;
  }

  const vec3r& n = faceNormal[f];
  const vec3r& origin = shp->vertex[poly[0]];
  double dn = (x - origin) * n;
  vec3r P = x - dn * n;

  // Crossing number (even-odd rule) in a 2D basis of the face plane
  vec3r u = shp->vertex[poly[1]] - origin;
  u.normalize();
  vec3r v = cross(n, u);
  double p1 = P * u;
  double p2 = P * v;
  bool odd = false;
  for (size_t ia = 0; ia < n_v; ++ia) {
    size_t ib = (ia + 1) % n_v;
    double a1 = shp->vertex[poly[ia]] * u;
    double a2 = shp->vertex[poly[ia]] * v;
    double b1 = shp->vertex[poly[ib]] * u;
    double b2 = shp->vertex[poly[ib]] * v;
    if ((a2 < p2 && b2 >= p2) || (b2 < p2 && a2 >= p2)) {
      if (a1 + (p2 - a2) / (b2 - a2) * (b1 - a1) < p1) {
        odd = !odd;
      }
    }
  }

  if (odd) {
    result.d2 = dn * dn;
    result.point = P;
    result.feature = SDF_FEATURE_INTERIOR;
    return result;
  }

  // Outside the polygon: the closest point sits on one of its boundary segments
  result.d2 = std::numeric_limits<double>::max();
  for (size_t ia = 0; ia < n_v; ++ia) {
    size_t ib = (ia + 1) % n_v;
    vec3r closest;
    double r;
    double d2 = sqDistToSegment(x, shp->vertex[poly[ia]], shp->vertex[poly[ib]], closest, r);
    if (d2 < result.d2) {
      result.d2 = d2;
      result.point = closest;
      if (r <= 0.0) {
        result.feature = SDF_FEATURE_VERTEX;
        result.ia = poly[ia];
      } else if (r >= 1.0) {
        result.feature = SDF_FEATURE_VERTEX;
        result.ia = poly[ib];
      } else {
        result.feature = SDF_FEATURE_EDGE;
        result.ia = poly[ia];
        result.ib = poly[ib];
      }
    }
  }
  return result;
}

double ShapeSDF::angleAtVertex(size_t f, size_t v) const {
  const std::vector<size_t>& poly = orientedFace[f];
  size_t n_v = poly.size();
  for (size_t i = 0; i < n_v; ++i) {
    if (poly[i] != v) {
      continue;
    }
    vec3r a = shp->vertex[poly[(i + n_v - 1) % n_v]] - shp->vertex[v];
    vec3r b = shp->vertex[poly[(i + 1) % n_v]] - shp->vertex[v];
    double na = norm(a);
    double nb = norm(b);
    if (na < 1.0e-15 || nb < 1.0e-15) {
      return 0.0;
    }
    double c = (a * b) / (na * nb);
    if (c < -1.0) {
      c = -1.0;
    } else if (c > 1.0) {
      c = 1.0;
    }
    return std::acos(c);
  }
  return 0.0;
}

// The pseudonormal is the face normal in the interior of a face, the sum of the
// two adjacent face normals on a segment, and the angle-weighted sum of the
// incident face normals at a vertex. It is exactly what makes the sign right
// when the closest point falls on a crease of the polyhedron
vec3r ShapeSDF::pseudoNormal(const FaceClosest& fc, size_t f) const {
  vec3r N;

  if (fc.feature == SDF_FEATURE_INTERIOR) {
    N = faceNormal[f];
  } else if (fc.feature == SDF_FEATURE_EDGE) {
    std::map<std::pair<size_t, size_t>, std::vector<size_t> >::const_iterator it =
        edgeFaces.find(sortedPair(fc.ia, fc.ib));
    if (it == edgeFaces.end()) {
      return faceNormal[f];
    }
    for (size_t k = 0; k < it->second.size(); ++k) {
      N += faceNormal[it->second[k]];
    }
  } else {  // SDF_FEATURE_VERTEX
    const std::vector<size_t>& incident = vertexFaces[fc.ia];
    for (size_t k = 0; k < incident.size(); ++k) {
      N += angleAtVertex(incident[k], fc.ia) * faceNormal[incident[k]];
    }
  }

  if (N.isnull()) {
    return faceNormal[f];
  }
  N.normalize();
  return N;
}

// Nearest face to x, through the face-only BVH. Branch-and-bound: a subtree is
// skipped when its box is already farther than the best face found so far
bool ShapeSDF::nearestFace(const vec3r& x, FaceClosest& out, size_t& faceOut) const {
  if (faceNodes.empty()) return false;

  double best_d2 = std::numeric_limits<double>::max();
  bool found = false;
  int stack[64];
  int sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const BVHNode& node = faceNodes[stack[--sp]];
    if (sqDistToBox(x, node.box) >= best_d2) continue;
    if (node.count > 0) {
      for (int i = 0; i < node.count; ++i) {
        size_t f = prims[faceIds[node.start + i]].index;
        FaceClosest fc = closestOnFace(x, f);
        if (fc.d2 < best_d2) {
          best_d2 = fc.d2;
          out = fc;
          faceOut = f;
          found = true;
        }
      }
    } else {
      stack[sp++] = node.left;
      stack[sp++] = node.right;
    }
  }
  return found;
}

double ShapeSDF::signAt(const vec3r& x) const {
  if (rayCastSign) {
    return rayCastInside(x) ? -1.0 : 1.0;
  }
  if (!closedVolume) {
    return 1.0;
  }

  FaceClosest best;
  size_t best_f = 0;
  if (!nearestFace(x, best, best_f)) return 1.0;

  vec3r N = pseudoNormal(best, best_f);
  return ((x - best.point) * N < 0.0) ? -1.0 : 1.0;
}

ClosestPoint ShapeSDF::closestPoint(const vec3r& x) const {
  ClosestPoint result;
  double best_d2 = std::numeric_limits<double>::max();
  if (bvhNodes.empty()) {
    result.distance = 0.0;
    return result;
  }

  // Branch-and-bound over the all-primitive BVH
  int stack[64];
  int sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const BVHNode& node = bvhNodes[stack[--sp]];
    if (sqDistToBox(x, node.box) >= best_d2) continue;

    if (node.count > 0) {
      for (int i = 0; i < node.count; ++i) {
        const Prim& p = prims[bvhIds[node.start + i]];
        if (p.type == SDF_PRIM_VERTEX) {
          double d2 = norm2(x - shp->vertex[p.index]);
          if (d2 < best_d2) {
            best_d2 = d2;
            result.point = shp->vertex[p.index];
            result.primType = SDF_PRIM_VERTEX;
            result.primIndex = p.index;
          }
        } else if (p.type == SDF_PRIM_EDGE) {
          vec3r closest;
          double r;
          double d2 = sqDistToSegment(x, shp->vertex[shp->edge[p.index].first],
                                      shp->vertex[shp->edge[p.index].second], closest, r);
          if (d2 < best_d2) {
            best_d2 = d2;
            result.point = closest;
            result.primType = SDF_PRIM_EDGE;
            result.primIndex = p.index;
          }
        } else {  // SDF_PRIM_FACE
          FaceClosest fc = closestOnFace(x, p.index);
          if (fc.d2 < best_d2) {
            best_d2 = fc.d2;
            result.point = fc.point;
            result.primType = SDF_PRIM_FACE;
            result.primIndex = p.index;
          }
        }
      }
    } else {
      stack[sp++] = node.left;
      stack[sp++] = node.right;
    }
  }

  result.distance = std::sqrt(best_d2);
  return result;
}

double ShapeSDF::unsignedSkeletonDistance(const vec3r& x) const { return closestPoint(x).distance; }

double ShapeSDF::signedSkeletonDistance(const vec3r& x) const { return signAt(x) * closestPoint(x).distance; }

double ShapeSDF::value(const vec3r& x) const { return signedSkeletonDistance(x) - shp->radius; }

vec3r ShapeSDF::gradient(const vec3r& x) const {
  ClosestPoint cp = closestPoint(x);

  if (cp.distance < 1.0e-15) {
    // On the skeleton itself. This only ever happens when radius == 0, since
    // otherwise the zero level set stands a distance radius away from it
    return vec3r();
  }

  return (signAt(x) / cp.distance) * (x - cp.point);
}
