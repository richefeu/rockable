// Phase 1 validation of the skin mesher (buildShapeMesh).
//
// Checks:
//   (a) watertight   : every triangle edge is shared by exactly two triangles,
//                      and the Euler characteristic matches the expected genus
//   (b) on-surface   : every mesh vertex lies on f = 0 to rounding, and its
//                      stored normal equals the SDF gradient there
//   (c) convergence  : the mesh volume converges to the exact / Monte-Carlo
//                      volume, at the O(eps^2) rate the sag argument predicts

#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "Core/Shape.hpp"
#include "Core/ShapeMesh.hpp"
#include "Core/ShapeSDF.hpp"
#include "testShapes.hpp"

static int g_failed = 0;
static int g_passed = 0;

static void check(bool ok, const std::string& what, double measured = 0.0, double tol = 0.0) {
  if (ok) {
    ++g_passed;
  } else {
    ++g_failed;
    std::printf("  FAIL  %-52s  measured %.4e  tol %.4e\n", what.c_str(), measured, tol);
  }
}

// Shape builders live in testShapes.hpp, shared with the SDF test.

// ---------------------------------------------------------------- exact volumes

static double volSphere(double R) { return 4.0 / 3.0 * M_PI * R * R * R; }

static double volCapsule(double R, double L) {  // cylinder + two half-balls
  return M_PI * R * R * L + 4.0 / 3.0 * M_PI * R * R * R;
}

// Volume of a cube of half-side a dilated by R: box + 6 slabs + 12 quarter
// cylinders + 8 eighth-spheres
static double volRoundedCube(double a, double R) {
  double s = 2.0 * a;
  double box = s * s * s;
  double faces = 6.0 * (s * s) * R;
  double edges = 12.0 * (0.25 * M_PI * R * R) * s;
  double corners = 4.0 / 3.0 * M_PI * R * R * R;
  return box + faces + edges + corners;
}

// ---------------------------------------------------------------- mesh checks

// A watertight triangle mesh has every undirected edge shared by exactly two
// triangles. Reports the worst offender count
static void checkWatertight(const ShapeMesh& m, const std::string& label) {
  std::map<std::pair<size_t, size_t>, int> edgeCount;
  for (size_t t = 0; t < m.tri.size(); ++t) {
    size_t idx[3] = {m.tri[t].a, m.tri[t].b, m.tri[t].c};
    for (int k = 0; k < 3; ++k) {
      size_t u = idx[k], v = idx[(k + 1) % 3];
      if (u > v) std::swap(u, v);
      ++edgeCount[std::make_pair(u, v)];
    }
  }
  int bad = 0;
  for (std::map<std::pair<size_t, size_t>, int>::const_iterator it = edgeCount.begin(); it != edgeCount.end(); ++it) {
    if (it->second != 2) ++bad;
  }
  check(bad == 0, label + ": watertight (every edge shared by 2 triangles)", (double)bad);

  // Euler characteristic V - E + F = 2 for a genus-0 closed surface
  long V = (long)m.P.size();
  long E = (long)edgeCount.size();
  long F = (long)m.tri.size();
  check(V - E + F == 2, label + ": Euler characteristic == 2", (double)(V - E + F));
  std::printf("  %-20s V=%ld E=%ld F=%ld  chi=%ld\n", label.c_str(), V, E, F, V - E + F);
}

// Every mesh vertex must sit on f = 0, and its stored normal must be the SDF
// gradient there
static void checkOnSurface(const Shape& s, const ShapeMesh& m, const std::string& label) {
  ShapeSDF sdf(s);
  double maxF = 0.0;
  double maxNdiff = 0.0;
  for (size_t i = 0; i < m.P.size(); ++i) {
    double f = std::fabs(sdf.value(m.P[i]));
    if (f > maxF) maxF = f;
    vec3r g = sdf.gradient(m.P[i]);
    double nd = norm(g - m.N[i]);
    if (nd > maxNdiff) maxNdiff = nd;
  }
  // Newton converges quadratically; two iterations from a linear guess on a
  // grid of step h leave a residual well below the tolerance
  check(maxF < 1.0e-6, label + ": vertices lie on f = 0", maxF, 1.0e-6);
  check(maxNdiff < 1.0e-9, label + ": stored normals equal SDF gradient", maxNdiff, 1.0e-9);
  std::printf("  %-20s max|f| = %.3e  max|dN| = %.3e\n", label.c_str(), maxF, maxNdiff);
}

// The headline phase-1 result: the mesh volume converges to the exact volume,
// at second order in the grid step h. After projection the mesh is a polyhedron
// inscribed in the smooth surface, whose volume deficit is O(h^2); since
// h ~ sqrt(8 eps R) this is O(eps) in the tolerance. We drive the grid by h
// directly, halving it between levels, and expect the error to drop about 4x
static void checkVolumeConvergence(const Shape& s, double exactVol, const std::string& label) {
  double R = s.radius;
  // Choose eps so that h halves at each level (h ~ sqrt(8 eps R), so eps /= 4),
  // and stay fine enough that the thin-feature step cap (h <= R/2, active for
  // eps > R/32 on non-solid shapes) never clips the coarsest level
  double eps[3] = {4.0e-3, 1.0e-3, 2.5e-4};
  double err[3];
  std::printf("  %s: convergence of the mesh volume (exact = %.6f)\n", label.c_str(), exactVol);
  for (int i = 0; i < 3; ++i) {
    ShapeMeshOptions opt;
    opt.epsilon = eps[i];
    ShapeMesh m = buildShapeMesh(s, opt);
    err[i] = std::fabs(m.volume - exactVol) / exactVol;
    std::printf("      eps=%.5f  h~%.4f  nv=%7zu  nt=%7zu  vol=%.6f  relErr=%.3e\n", eps[i],
                std::sqrt(8.0 * eps[i] * R), m.nbVertices(), m.nbTriangles(), m.volume, err[i]);
  }

  check(err[1] < err[0] && err[2] < err[1], label + ": volume error decreases monotonically", err[2]);

  // Order in h: halving h (eps /= 4) should roughly quarter the error, i.e.
  // order 2 in h. Measured order in eps ~ log-ratio / log(4); expect ~1
  double pEps = std::log(err[0] / err[2]) / std::log(eps[0] / eps[2]);
  double pH = 2.0 * pEps;  // h ratio is the square root of the eps ratio
  std::printf("      observed order: %.2f in eps, %.2f in h\n", pEps, pH);
  check(pH > 1.6, label + ": second-order convergence in h", pH, 1.6);
}

// A quiet Monte-Carlo estimate of the r-shape volume, straight from
// Shape::inside(), used as an independent reference where no closed form exists.
// Independent of ShapeSDF, so it is a genuine cross-check of the mesher
static double monteCarloVolume(Shape& s, size_t nstep) {
  AABB box;
  s.getAABB(box);  // already enlarged by radius
  double Vbox = (box.max.x - box.min.x) * (box.max.y - box.min.y) * (box.max.z - box.min.z);
  std::mt19937 gen(20260710);
  std::uniform_real_distribution<double> Ux(box.min.x, box.max.x);
  std::uniform_real_distribution<double> Uy(box.min.y, box.max.y);
  std::uniform_real_distribution<double> Uz(box.min.z, box.max.z);
  size_t nin = 0;
  for (size_t i = 0; i < nstep; ++i) {
    if (s.inside(vec3r(Ux(gen), Uy(gen), Uz(gen)))) ++nin;
  }
  return Vbox * (double)nin / (double)nstep;
}

// Convergence for shapes with no closed-form volume (concave L, loaded .shp):
// refine the grid twice, take the finest mesh as the reference, and check that
// the two coarser meshes approach it monotonically at second order in h. Also
// cross-check the finest mesh against a Monte-Carlo estimate
static void checkVolumeConvergenceRef(Shape& s, const std::string& label) {
  double R = s.radius;
  double eps[3] = {1.6e-2, 4.0e-3, 1.0e-3};
  double vol[3];
  double refVol = 0.0;
  std::printf("  %s: convergence toward the finest mesh\n", label.c_str());
  for (int i = 0; i < 3; ++i) {
    ShapeMeshOptions opt;
    opt.epsilon = eps[i];
    ShapeMesh m = buildShapeMesh(s, opt);
    vol[i] = m.volume;
    std::printf("      eps=%.5f  h~%.4f  nv=%7zu  nt=%7zu  vol=%.6f\n", eps[i], std::sqrt(8.0 * eps[i] * R),
                m.nbVertices(), m.nbTriangles(), m.volume);
  }
  // A finer reference, half the step of the finest sampled level
  {
    ShapeMeshOptions opt;
    opt.epsilon = 2.5e-4;
    refVol = buildShapeMesh(s, opt).volume;
  }
  double e0 = std::fabs(vol[0] - refVol);
  double e2 = std::fabs(vol[2] - refVol);
  check(std::fabs(vol[1] - refVol) < e0 && e2 < std::fabs(vol[1] - refVol),
        label + ": volume converges monotonically to the fine reference", e2);
  double pH = 2.0 * std::log(e0 / e2) / std::log(eps[0] / eps[2]);
  std::printf("      ref vol (eps=2.5e-4) = %.6f   observed order ~ %.2f in h\n", refVol, pH);
  check(pH > 1.5, label + ": second-order convergence in h", pH, 1.5);

  double mc = monteCarloVolume(s, 4000000);
  double relDiff = std::fabs(refVol - mc) / mc;
  std::printf("      Monte-Carlo vol = %.6f   |ref - MC|/MC = %.3e\n", mc, relDiff);
  check(relDiff < 5.0e-3, label + ": fine mesh agrees with Monte-Carlo volume", relDiff, 5.0e-3);
}

// Coplanar simplification must keep the surface: watertight, same volume and
// area to rounding. reduceExpected says whether flat regions should shrink the
// triangle count (false for a purely curved shape like the sphere)
static void checkSimplify(Shape& s, bool reduceExpected, const std::string& label) {
  ShapeMeshOptions o;
  o.epsilon = 4.0e-3;
  ShapeMesh full = buildShapeMesh(s, o);
  ShapeMesh simp = simplifyCoplanar(full);

  // watertight
  std::map<std::pair<size_t, size_t>, int> ec;
  for (size_t t = 0; t < simp.tri.size(); ++t) {
    size_t v[3] = {simp.tri[t].a, simp.tri[t].b, simp.tri[t].c};
    for (int k = 0; k < 3; ++k) {
      size_t a = v[k], b = v[(k + 1) % 3];
      if (a > b) std::swap(a, b);
      ++ec[std::make_pair(a, b)];
    }
  }
  int bad = 0;
  for (std::map<std::pair<size_t, size_t>, int>::const_iterator it = ec.begin(); it != ec.end(); ++it) {
    if (it->second != 2) ++bad;
  }
  long chi = (long)simp.P.size() - (long)ec.size() + (long)simp.tri.size();
  check(bad == 0 && chi == 2, label + ": simplified mesh stays watertight (chi=2)", (double)chi);

  double dv = std::fabs(simp.volume - full.volume) / full.volume;
  double da = std::fabs(simp.area - full.area) / full.area;
  check(dv < 1.0e-3, label + ": volume preserved by simplification", dv, 1.0e-3);
  check(da < 1.0e-3, label + ": area preserved by simplification", da, 1.0e-3);

  double ratio = (double)simp.nbTriangles() / (double)full.nbTriangles();
  if (reduceExpected) {
    check(ratio < 0.85, label + ": flat regions reduce the triangle count", ratio, 0.85);
  } else {
    check(ratio <= 1.0 + 1.0e-9, label + ": curved shape is not inflated", ratio);
  }
  std::printf("  %-20s tri %zu -> %zu (%.0f%%), dVol=%.2e dArea=%.2e chi=%ld\n", label.c_str(), full.nbTriangles(),
              simp.nbTriangles(), 100.0 * ratio, dv, da, chi);
}

// ---------------------------------------------------------------- main

int main(int argc, char** argv) {
  Shape sph = makeSphere(0.5);
  Shape cap = makeCapsule(0.3, 1.0);
  Shape cub = makeCube(0.15, 0.5);
  Shape ell = makeL(0.1);

  std::printf("\n=== watertightness and Euler characteristic ===\n");
  {
    ShapeMeshOptions opt;
    opt.epsilon = 1.0e-2;
    ShapeMesh ms = buildShapeMesh(sph, opt);
    ShapeMesh mc = buildShapeMesh(cap, opt);
    ShapeMesh mb = buildShapeMesh(cub, opt);
    ShapeMesh ml = buildShapeMesh(ell, opt);
    checkWatertight(ms, "sphere");
    checkWatertight(mc, "capsule");
    checkWatertight(mb, "cube");
    checkWatertight(ml, "L (concave)");

    std::printf("\n=== vertices on the surface, exact normals ===\n");
    checkOnSurface(sph, ms, "sphere");
    checkOnSurface(cap, mc, "capsule");
    checkOnSurface(cub, mb, "cube");
    checkOnSurface(ell, ml, "L (concave)");
  }

  std::printf("\n=== volume convergence (closed form) ===\n");
  checkVolumeConvergence(sph, volSphere(0.5), "sphere");
  checkVolumeConvergence(cap, volCapsule(0.3, 1.0), "capsule");
  checkVolumeConvergence(cub, volRoundedCube(0.5, 0.15), "cube");
  // Same solid as the cube, but meshed through the ray-cast sign path
  Shape soup = makeCubeSoup(0.15, 0.5);
  checkVolumeConvergence(soup, volRoundedCube(0.5, 0.15), "cubeSoup (ray-cast)");

  std::printf("\n=== volume convergence (no closed form) ===\n");
  checkVolumeConvergenceRef(ell, "L (concave)");

  std::printf("\n=== coplanar simplification ===\n");
  checkSimplify(cub, true, "cube");
  checkSimplify(ell, true, "L (concave)");
  checkSimplify(sph, false, "sphere");

  // A real polyhedron loaded from a .shp file. The mesh is validated against a
  // fresh Monte-Carlo estimate of the same r-shape (an independent oracle). The
  // volume field stored in the file is only reported, not asserted: it can be
  // stale, as it is for POLYR0 in test/input/518_poly/shape.shp, where the
  // stored value matches neither the r-shape nor the bare polyhedron
  if (argc > 1) {
    Shape poly = loadFirstShape(argv[1]);
    if (!poly.name.empty()) {
      std::printf("\n=== real shape from %s (name '%s') ===\n", argv[1], poly.name.c_str());
      ShapeMeshOptions opt;
      opt.epsilon = 1.0e-1 * poly.radius;  // target sag ~ 10% of the Minkowski radius
      if (opt.epsilon <= 0.0) opt.epsilon = 1.0e-3;
      ShapeMesh m = buildShapeMesh(poly, opt);
      checkWatertight(m, poly.name);
      checkOnSurface(poly, m, poly.name);

      double mc = monteCarloVolume(poly, 8000000);
      double relDiff = std::fabs(m.volume - mc) / mc;
      std::printf("  mesh vol = %.6e   fresh MC vol = %.6e   relDiff = %.3e   (stored in file = %.6e)\n", m.volume, mc,
                  relDiff, poly.volume);
      // Inscribed in a convex body, so it must not exceed the true volume (bar
      // Monte-Carlo noise), and must land within ~2% of it
      check(m.volume <= mc * 1.005, poly.name + ": inscribed mesh does not exceed the MC volume");
      check(relDiff < 2.0e-2, poly.name + ": mesh volume within 2% of a fresh MC volume", relDiff, 2.0e-2);
    } else {
      std::printf("\n(could not load a shape from %s, skipping the real-shape test)\n", argv[1]);
    }
  }

  std::printf("\n%d passed, %d failed\n\n", g_passed, g_failed);
  return (g_failed == 0) ? 0 : 1;
}
