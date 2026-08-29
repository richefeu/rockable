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
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <queue>
#include <vector>

#include "ShapeMesh.hpp"

namespace {

// Signed area of a 2D polygon given by its vertices in order
double signedArea2D(const std::vector<vec3r>& poly, const std::vector<int>& idx, const vec3r& u, const vec3r& v,
                    const vec3r& origin) {
  double a = 0.0;
  size_t n = idx.size();
  for (size_t i = 0; i < n; ++i) {
    const vec3r& pa = poly[idx[i]];
    const vec3r& pb = poly[idx[(i + 1) % n]];
    double ax = (pa - origin) * u, ay = (pa - origin) * v;
    double bx = (pb - origin) * u, by = (pb - origin) * v;
    a += ax * by - bx * ay;
  }
  return 0.5 * a;
}

struct V2 {
  double x, y;
};

double cross2(const V2& a, const V2& b, const V2& c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// Is p strictly inside the CCW triangle (a,b,c)? Points that merely lie on an
// edge (e.g. collinear boundary vertices along a straight side) do not count, so
// they never wrongly block an otherwise valid ear.
bool strictlyInside(const V2& p, const V2& a, const V2& b, const V2& c, double eps) {
  return cross2(a, b, p) > eps && cross2(b, c, p) > eps && cross2(c, a, p) > eps;
}

// Ear-clip a simple polygon given as 2D points in order. Returns triangles as
// triplets of indices into pts. Keeps every boundary vertex, so the shared edges
// with neighbouring regions are preserved and the mesh stays watertight.
//
// Robustness: every clipped ear is verified empty (no other vertex strictly
// inside), so a triangle can never cross the polygon --- in particular it cannot
// bridge a reentrant notch, which used to happen when the old "most convex"
// fallback clipped without an emptiness check. Marching-cubes boundaries carry
// long runs of collinear vertices; these are handled by (i) a strict emptiness
// test that ignores on-edge points and (ii) a second pass that may clip a
// non-reflex (possibly collinear) empty vertex to get past such a run. If no
// valid ear remains, we stop and let the caller keep the region's original
// triangles rather than emit anything unsafe.
std::vector<std::array<int, 3> > earClip(const std::vector<V2>& pts) {
  std::vector<std::array<int, 3> > tris;
  size_t n = pts.size();
  if (n < 3) return tris;

  std::vector<int> V(n);
  double area = 0.0;
  for (size_t i = 0; i < n; ++i) area += pts[i].x * pts[(i + 1) % n].y - pts[(i + 1) % n].x * pts[i].y;
  for (size_t i = 0; i < n; ++i) V[i] = (area >= 0.0) ? (int)i : (int)(n - 1 - i);  // make CCW

  // Scale-aware tolerance for the length^2 cross products
  double scale = 0.0;
  for (size_t i = 0; i < n; ++i) scale = std::max(scale, std::max(std::fabs(pts[i].x), std::fabs(pts[i].y)));
  double eps = 1.0e-9 * (scale > 0.0 ? scale * scale : 1.0);

  int guard = 0;
  int guardMax = 6 * (int)n + 12;
  while (V.size() > 3 && guard++ < guardMax) {
    size_t m = V.size();

    // Return the index (in V) of an ear to clip, or -1 if none. strictOnly means
    // the apex must be strictly convex; otherwise a non-reflex (collinear) apex
    // is also accepted, to advance along a straight run of vertices.
    auto findEar = [&](bool strictOnly) -> int {
      int best = -1;
      double bestCr = 0.0;
      for (size_t i = 0; i < m; ++i) {
        int ip = V[(i + m - 1) % m], ic = V[i], in = V[(i + 1) % m];
        const V2 &a = pts[ip], &b = pts[ic], &c = pts[in];
        double cr = cross2(a, b, c);
        if (strictOnly ? (cr <= eps) : (cr < -eps)) continue;  // reflex (or, in pass A, collinear)
        bool empty = true;
        for (size_t k = 0; k < m; ++k) {
          int iv = V[k];
          if (iv == ip || iv == ic || iv == in) continue;
          if (strictlyInside(pts[iv], a, b, c, eps)) { empty = false; break; }
        }
        if (!empty) continue;
        if (strictOnly) {
          if (best < 0 || cr > bestCr) { best = (int)i; bestCr = cr; }  // largest ear first
        } else {
          return (int)i;  // any safe ear, to get past a collinear run
        }
      }
      return best;
    };

    int i = findEar(true);
    if (i < 0) i = findEar(false);
    if (i < 0) break;  // no safe ear: stop (caller keeps the original triangles)

    size_t mm = V.size();
    int ip = V[(i + mm - 1) % mm], ic = V[i], in = V[(i + 1) % mm];
    tris.push_back({ip, ic, in});
    V.erase(V.begin() + i);
  }
  if (V.size() == 3) tris.push_back({V[0], V[1], V[2]});
  return tris;
}

}  // namespace

ShapeMesh simplifyCoplanar(const ShapeMesh& in, const SimplifyOptions& opt) {
  ShapeMesh out;
  size_t nt = in.tri.size();
  if (nt == 0) return in;

  // Bounding box diagonal, for the absolute plane tolerance
  vec3r lo = in.P[0], hi = in.P[0];
  for (size_t i = 1; i < in.P.size(); ++i) {
    lo = component_min(lo, in.P[i]);
    hi = component_max(hi, in.P[i]);
  }
  double diag = norm(hi - lo);
  // Bound how far a dropped vertex may sit from the region plane by a fraction of
  // the mesh's own sag, so the merge error scales with the discretisation error
  // instead of staying fixed while the mesh refines. Fall back to a bounding-box
  // fraction only for a mesh that carries no sag.
  double planeTol = (in.epsilon > 0.0) ? opt.planeToleranceSag * in.epsilon
                                       : opt.planeToleranceRel * (diag > 0.0 ? diag : 1.0);
  double cosTol = std::cos(opt.angleToleranceDeg * M_PI / 180.0);

  // Per-triangle unit normal
  std::vector<vec3r> triN(nt);
  for (size_t t = 0; t < nt; ++t) {
    const vec3r& a = in.P[in.tri[t].a];
    const vec3r& b = in.P[in.tri[t].b];
    const vec3r& c = in.P[in.tri[t].c];
    vec3r n = cross(b - a, c - a);
    double ln = norm(n);
    triN[t] = (ln > 0.0) ? (1.0 / ln) * n : vec3r(0, 0, 1);
  }

  // Triangle adjacency across shared edges
  std::map<std::pair<size_t, size_t>, std::vector<int> > edgeTris;
  for (size_t t = 0; t < nt; ++t) {
    size_t v[3] = {in.tri[t].a, in.tri[t].b, in.tri[t].c};
    for (int k = 0; k < 3; ++k) {
      size_t p = v[k], q = v[(k + 1) % 3];
      if (p > q) std::swap(p, q);
      edgeTris[std::make_pair(p, q)].push_back((int)t);
    }
  }
  std::vector<std::array<int, 3> > neigh(nt, {-1, -1, -1});
  for (size_t t = 0; t < nt; ++t) {
    size_t v[3] = {in.tri[t].a, in.tri[t].b, in.tri[t].c};
    for (int k = 0; k < 3; ++k) {
      size_t p = v[k], q = v[(k + 1) % 3];
      if (p > q) std::swap(p, q);
      const std::vector<int>& ts = edgeTris[std::make_pair(p, q)];
      for (size_t j = 0; j < ts.size(); ++j) {
        if (ts[j] != (int)t) {
          neigh[t][k] = ts[j];
          break;
        }
      }
    }
  }

  // Grow planar regions. The seed plane (point + normal) is fixed so a chain of
  // slightly tilted triangles cannot drift away from the plane
  std::vector<int> region(nt, -1);
  std::vector<vec3r> regNormal, regOrigin;
  for (size_t s = 0; s < nt; ++s) {
    if (region[s] != -1) continue;
    int rid = (int)regNormal.size();
    vec3r rN = triN[s];
    vec3r rO = in.P[in.tri[s].a];
    regNormal.push_back(rN);
    regOrigin.push_back(rO);

    std::queue<int> todo;
    todo.push((int)s);
    region[s] = rid;
    while (!todo.empty()) {
      int t = todo.front();
      todo.pop();
      for (int k = 0; k < 3; ++k) {
        int g = neigh[t][k];
        if (g < 0 || region[g] != -1) continue;
        if (triN[g] * rN < cosTol) continue;  // normals not aligned
        // All three vertices of g must lie on the seed plane
        double d0 = std::fabs((in.P[in.tri[g].a] - rO) * rN);
        double d1 = std::fabs((in.P[in.tri[g].b] - rO) * rN);
        double d2 = std::fabs((in.P[in.tri[g].c] - rO) * rN);
        if (d0 > planeTol || d1 > planeTol || d2 > planeTol) continue;
        region[g] = rid;
        todo.push(g);
      }
    }
  }

  // Assemble the output: re-triangulate each multi-triangle region from its
  // boundary loop, keep single-triangle regions as they are
  std::vector<std::vector<int> > regionTris(regNormal.size());
  for (size_t t = 0; t < nt; ++t) regionTris[region[t]].push_back((int)t);

  std::vector<Tri> newTris;

  for (size_t r = 0; r < regionTris.size(); ++r) {
    const std::vector<int>& rts = regionTris[r];
    if (rts.size() == 1) {
      newTris.push_back(in.tri[rts[0]]);
      continue;
    }

    // Directed boundary edges of the region: a directed edge (u->v) of a region
    // triangle whose opposite triangle is outside the region
    std::map<size_t, size_t> nextOf;
    bool simple = true;
    size_t nBoundary = 0;
    for (size_t j = 0; j < rts.size() && simple; ++j) {
      int t = rts[j];
      size_t v[3] = {in.tri[t].a, in.tri[t].b, in.tri[t].c};
      for (int k = 0; k < 3; ++k) {
        int g = neigh[t][k];
        if (g >= 0 && region[g] == (int)r) continue;  // interior edge
        size_t a = v[k], b = v[(k + 1) % 3];
        if (nextOf.count(a)) {
          simple = false;  // a vertex with two outgoing boundary edges: not a simple loop
          break;
        }
        nextOf[a] = b;
        ++nBoundary;
      }
    }

    // Follow the directed edges into a single closed loop
    std::vector<int> loop;
    if (simple && nBoundary >= 3) {
      size_t start = nextOf.begin()->first;
      size_t cur = start;
      for (size_t step = 0; step <= nBoundary; ++step) {
        loop.push_back((int)cur);
        std::map<size_t, size_t>::iterator it = nextOf.find(cur);
        if (it == nextOf.end()) {
          simple = false;
          break;
        }
        cur = it->second;
        if (cur == start) break;
      }
      if (cur != start || loop.size() != nBoundary) simple = false;  // several loops (holes)
    } else {
      simple = false;
    }

    if (!simple) {  // safe fallback: keep the region's original triangles
      for (size_t j = 0; j < rts.size(); ++j) newTris.push_back(in.tri[rts[j]]);
      continue;
    }

    // Re-triangulate the loop in the region plane
    const vec3r& n = regNormal[r];
    vec3r u = (std::fabs(n.x) < 0.9) ? cross(n, vec3r(1, 0, 0)) : cross(n, vec3r(0, 1, 0));
    u.normalize();
    vec3r vv = cross(n, u);  // u x vv aligned with n, so a CCW loop faces outward
    const vec3r& o = in.P[loop[0]];
    std::vector<V2> pts(loop.size());
    for (size_t i = 0; i < loop.size(); ++i) {
      pts[i].x = (in.P[loop[i]] - o) * u;
      pts[i].y = (in.P[loop[i]] - o) * vv;
    }
    (void)signedArea2D;  // orientation is handled inside earClip

    std::vector<std::array<int, 3> > local = earClip(pts);
    if (local.size() != loop.size() - 2) {
      // earClip could not fully triangulate this loop: keep the original
      // triangles rather than emit an incomplete (leaky) region.
      for (size_t j = 0; j < rts.size(); ++j) newTris.push_back(in.tri[rts[j]]);
      continue;
    }
    for (size_t i = 0; i < local.size(); ++i) {
      Tri tr;
      tr.a = (size_t)loop[local[i][0]];
      tr.b = (size_t)loop[local[i][1]];
      tr.c = (size_t)loop[local[i][2]];
      newTris.push_back(tr);
    }
  }

  // Compact the vertices actually referenced by the new triangles
  std::vector<int> remap(in.P.size(), -1);
  for (size_t t = 0; t < newTris.size(); ++t) {
    size_t idx[3] = {newTris[t].a, newTris[t].b, newTris[t].c};
    for (int k = 0; k < 3; ++k) {
      if (remap[idx[k]] == -1) {
        remap[idx[k]] = (int)out.P.size();
        out.P.push_back(in.P[idx[k]]);
        out.N.push_back(in.N[idx[k]]);
        if (idx[k] < in.primType.size()) out.primType.push_back(in.primType[idx[k]]);
        if (idx[k] < in.primIndex.size()) out.primIndex.push_back(in.primIndex[idx[k]]);
      }
    }
  }
  out.tri.resize(newTris.size());
  for (size_t t = 0; t < newTris.size(); ++t) {
    out.tri[t].a = (size_t)remap[newTris[t].a];
    out.tri[t].b = (size_t)remap[newTris[t].b];
    out.tri[t].c = (size_t)remap[newTris[t].c];
  }

  out.epsilon = in.epsilon;
  out.volume = out.computeVolume();
  out.area = out.computeArea();
  return out;
}
