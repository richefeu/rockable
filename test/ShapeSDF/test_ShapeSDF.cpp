// Phase 0 validation of ShapeSDF.
//
// Three families of checks:
//   (a) analytic  : the SDF of shapes whose exact distance function is known
//                   by hand (sphere, capsule, plate, cube)
//   (b) consistency: sign(value(x)) must agree with Shape::inside(x)
//   (c) gradient  : the analytic gradient must match a centered finite
//                   difference of value(), and be of unit length

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "Core/Shape.hpp"
#include "Core/ShapeSDF.hpp"
#include "testShapes.hpp"

static int g_failed = 0;
static int g_passed = 0;

static void check(bool ok, const std::string& what, double measured = 0.0, double tol = 0.0) {
  if (ok) {
    ++g_passed;
  } else {
    ++g_failed;
    std::printf("  FAIL  %-52s  measured %.3e  tol %.3e\n", what.c_str(), measured, tol);
  }
}

static void checkNear(double a, double b, double tol, const std::string& what) {
  check(std::fabs(a - b) <= tol, what, std::fabs(a - b), tol);
}

// Shape builders (makeSphere, makeCapsule, makePlate, makeCube, makeL) live in
// testShapes.hpp, shared with the ShapeMesh test.

// ---------------------------------------------------------------- exact SDFs

static double sdfSphere(const vec3r& x, double R) { return norm(x) - R; }

static double sdfCapsule(const vec3r& x, double R, double L) {
  vec3r a(-0.5 * L, 0, 0);
  vec3r b(0.5 * L, 0, 0);
  vec3r E = b - a;
  double r = ((x - a) * E) / (E * E);
  r = (r < 0.0) ? 0.0 : ((r > 1.0) ? 1.0 : r);
  return norm(x - (a + r * E)) - R;
}

// Exact SDF of a box of half-extent h, dilated by R (rounded box)
static double sdfRoundedBox(const vec3r& x, double h, double R) {
  vec3r q(std::fabs(x.x) - h, std::fabs(x.y) - h, std::fabs(x.z) - h);
  vec3r qp(std::max(q.x, 0.0), std::max(q.y, 0.0), std::max(q.z, 0.0));
  double outside = norm(qp);
  double inside = std::min(std::max(q.x, std::max(q.y, q.z)), 0.0);
  return outside + inside - R;
}

// Exact SDF of the square plate [-a,a]^2 in z=0, dilated by R
static double sdfPlate(const vec3r& x, double a, double R) {
  double dx = std::fabs(x.x) - a;
  double dy = std::fabs(x.y) - a;
  double px = std::max(dx, 0.0);
  double py = std::max(dy, 0.0);
  double inPlane = std::sqrt(px * px + py * py);
  if (dx < 0.0 && dy < 0.0) {
    inPlane = 0.0;
  }
  return std::sqrt(inPlane * inPlane + x.z * x.z) - R;
}

// ---------------------------------------------------------------- test drivers

typedef double (*ExactSDF)(const vec3r&);
static double g_R = 0.0, g_L = 0.0, g_a = 0.0;
static double exSphere(const vec3r& x) { return sdfSphere(x, g_R); }
static double exCapsule(const vec3r& x) { return sdfCapsule(x, g_R, g_L); }
static double exBox(const vec3r& x) { return sdfRoundedBox(x, g_a, g_R); }
static double exPlate(const vec3r& x) { return sdfPlate(x, g_a, g_R); }

// Compare the SDF against a closed-form expression on a random cloud
static void testAgainstExact(Shape& s, ExactSDF exact, double box, const std::string& label) {
  ShapeSDF sdf(s);
  std::mt19937 gen(12345);
  std::uniform_real_distribution<double> U(-box, box);

  double maxErr = 0.0;
  for (int i = 0; i < 20000; ++i) {
    vec3r x(U(gen), U(gen), U(gen));
    double err = std::fabs(sdf.value(x) - exact(x));
    if (err > maxErr) {
      maxErr = err;
    }
  }
  check(maxErr < 1.0e-12, label + ": value vs exact SDF", maxErr, 1.0e-12);
  std::printf("  %-24s max |f - f_exact| = %.3e\n", label.c_str(), maxErr);
}

// The whole point of the SDF is to reproduce Shape::inside() exactly, except
// within rounding distance of the surface where both are entitled to disagree
static void testAgainstInside(Shape& s, double box, const std::string& label) {
  ShapeSDF sdf(s);
  std::mt19937 gen(6789);
  std::uniform_real_distribution<double> U(-box, box);

  int mismatch = 0;
  int nearSurface = 0;
  double worst = 0.0;
  for (int i = 0; i < 50000; ++i) {
    vec3r x(U(gen), U(gen), U(gen));
    bool a = sdf.inside(x);
    bool b = s.inside(x);
    if (a != b) {
      double f = std::fabs(sdf.value(x));
      if (f < 1.0e-9) {
        ++nearSurface;  // a tie on the surface, nobody is wrong
      } else {
        ++mismatch;
        if (f > worst) {
          worst = f;
        }
      }
    }
  }
  check(mismatch == 0, label + ": inside() agrees with Shape::inside()", worst, 0.0);
  std::printf("  %-24s mismatches = %d (ties on surface: %d)\n", label.c_str(), mismatch, nearSurface);
}

// grad f must be of unit length and match a centered finite difference
static void testGradient(Shape& s, double box, const std::string& label) {
  ShapeSDF sdf(s);
  std::mt19937 gen(2024);
  std::uniform_real_distribution<double> U(-box, box);
  const double h = 1.0e-6;

  double maxNormErr = 0.0;
  double maxFdErr = 0.0;
  int nsample = 0;
  for (int i = 0; i < 5000; ++i) {
    vec3r x(U(gen), U(gen), U(gen));
    vec3r g = sdf.gradient(x);
    if (g.isnull()) {
      continue;
    }
    ++nsample;

    double ne = std::fabs(norm(g) - 1.0);
    if (ne > maxNormErr) {
      maxNormErr = ne;
    }

    // Finite differences straddle the creases of the distance function, where
    // it is not differentiable. Skip the samples where the two one-sided slopes
    // disagree, they say nothing about the analytic gradient
    vec3r fd;
    fd.x = (sdf.value(x + vec3r(h, 0, 0)) - sdf.value(x - vec3r(h, 0, 0))) / (2 * h);
    fd.y = (sdf.value(x + vec3r(0, h, 0)) - sdf.value(x - vec3r(0, h, 0))) / (2 * h);
    fd.z = (sdf.value(x + vec3r(0, 0, h)) - sdf.value(x - vec3r(0, 0, h))) / (2 * h);
    if (std::fabs(norm(fd) - 1.0) > 1.0e-4) {
      continue;  // on a crease
    }
    double e = norm(fd - g);
    if (e > maxFdErr) {
      maxFdErr = e;
    }
  }
  check(maxNormErr < 1.0e-12, label + ": |grad f| == 1", maxNormErr, 1.0e-12);
  check(maxFdErr < 1.0e-5, label + ": grad f vs finite differences", maxFdErr, 1.0e-5);
  std::printf("  %-24s |grad|-1 = %.3e, grad vs FD = %.3e (%d samples)\n", label.c_str(), maxNormErr, maxFdErr,
              nsample);
}

// ---------------------------------------------------------------- main

int main() {
  std::printf("\n=== analytic SDFs ===\n");

  g_R = 0.3;
  Shape sph = makeSphere(g_R);
  testAgainstExact(sph, exSphere, 1.0, "sphere");

  g_R = 0.2;
  g_L = 1.0;
  Shape cap = makeCapsule(g_R, g_L);
  testAgainstExact(cap, exCapsule, 1.2, "capsule");

  g_R = 0.15;
  g_a = 0.5;
  Shape cub = makeCube(g_R, g_a);
  testAgainstExact(cub, exBox, 1.2, "cube (rounded box)");

  g_R = 0.1;
  g_a = 0.6;
  Shape pla = makePlate(g_R, g_a);
  testAgainstExact(pla, exPlate, 1.2, "plate (surface)");

  std::printf("\n=== orientation and closedness ===\n");
  {
    ShapeSDF sdfCube(cub);
    check(sdfCube.isClosedVolume(), "cube: detected as a closed volume");
    // The cube faces were fed in with inconsistent winding on purpose
    checkNear(sdfCube.skeletonVolume(), 1.0, 1.0e-12, "cube: skeleton volume == 1");
    std::printf("  cube skeleton volume = %.15f (expected 1)\n", sdfCube.skeletonVolume());

    ShapeSDF sdfPlate(pla);
    check(!sdfPlate.isClosedVolume(), "plate: not a closed volume");
    ShapeSDF sdfCap(cap);
    check(!sdfCap.isClosedVolume(), "capsule: not a closed volume");
  }

  // A concave solid: the sign must stay right on the reentrant edge, which is
  // exactly where the angle-weighted pseudonormal earns its keep
  Shape ell = makeL(0.1);
  {
    ShapeSDF sdfL(ell);
    check(sdfL.isClosedVolume(), "L: detected as a closed volume");
    checkNear(sdfL.skeletonVolume(), 3.0, 1.0e-12, "L: skeleton volume == 3");
  }

  // A non-manifold face soup: not a clean closed volume, so the sign comes from
  // ray casting rather than the pseudonormal. This is the STL-soup path
  Shape soup = makeCubeSoup(0.15, 0.5);
  {
    ShapeSDF sdfSoup(soup);
    check(!sdfSoup.isClosedVolume(), "cubeSoup: not a clean closed volume (ray-cast sign)");
  }

  std::printf("\n=== consistency with Shape::inside() ===\n");
  g_R = 0.3;
  testAgainstInside(sph, 1.0, "sphere");
  g_R = 0.2;
  testAgainstInside(cap, 1.2, "capsule");
  testAgainstInside(cub, 1.2, "cube");
  testAgainstInside(pla, 1.2, "plate");
  testAgainstInside(ell, 1.5, "L (concave)");
  testAgainstInside(soup, 1.2, "cubeSoup (ray-cast)");

  std::printf("\n=== gradient ===\n");
  testGradient(sph, 1.0, "sphere");
  testGradient(cap, 1.2, "capsule");
  testGradient(cub, 1.2, "cube");
  testGradient(pla, 1.2, "plate");
  testGradient(ell, 1.5, "L (concave)");

  std::printf("\n%d passed, %d failed\n\n", g_passed, g_failed);
  return (g_failed == 0) ? 0 : 1;
}
