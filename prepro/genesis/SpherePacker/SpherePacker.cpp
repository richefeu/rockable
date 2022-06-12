#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "addParticle.hpp"
#include "generateShape_cylindricalMold.hpp"
#include "generateShape_sphere.hpp"
#include "generateShape_xyz_walls.hpp"

#define CONTAINER_BOX

using namespace std;

unsigned int seed = 143502;

// =========================================================
// =========================================================
double Rw = 0.5e-3;  // wall radius
double baseLength = 20.e-3; // inner width or diameter of the container

double Dmax = 1.2e-3;
double Dmin = 0.8e-3;

size_t nb_spheres = 4000;
size_t nb_trials_max = nb_spheres * 30;

double SolidFractionTarget = 0.29;
// =========================================================
// =========================================================

struct sphere {
  double x, y, z;
  double R;

  sphere(double _x, double _y, double _z, double _R) : x(_x), y(_y), z(_z), R(_R) {}
  sphere() : x(0.0), y(0.0), z(0.0), R(0.0) {}
};

struct BoxContainer {
  double L;
  double H;
};

struct CylContainer {
  double diam;
  double H;
};

void save_pack(const char* fname, vector<sphere>& g, double xmin, double ymin, double zmin, double xmax, double ymax,
               double zmax) {
  ofstream fog(fname, ios::out);
  if (fog) {
    fog.precision(5);
    fog << scientific;
    //fog << xmin << ' ' << ymin << ' ' << zmin << "\n";
    //fog << xmax << ' ' << ymax << ' ' << zmax << "\n";
    for (int i = 0; i < g.size(); i++) {
      fog << g[i].x << " " << g[i].y << " " << g[i].z << " " << g[i].R << endl;
    }
  }
}

void save_rockable(const char* fInputName, const char* fShapesName, vector<sphere>& g, double xmin, double ymin,
                   double zmin, double xmax, double ymax, double zmax) {
  ofstream fshape(fShapesName, ios::out);
  if (fshape) {
#if defined(CONTAINER_BOX)
    generateShape_xyz_walls(fshape, 1.1 * (xmax - xmin), 1.1 * (ymax - ymin), 1.1 * (zmax - zmin), Rw);
#else
    generateShape_cylindricalMold(fshape, "cylBox", Rw, ymax - ymin, 36, 0.5 * (xmax - xmin));
#endif
    generateShape_sphere(fshape, "sphere", 1.0);
  }

  ofstream fog(fInputName, ios::out);
  if (fog) {

    vec3r position;
    quat Q;
    int wallGroup = 1;
    int clusterId = 0;
#if defined(CONTAINER_BOX)
    fog << "Particles " << g.size() + 6 << std::endl;
    Q.reset();
    position.set(xmin - Rw, 0.5*(ymin+ymax), 0.0);
    addParticle(fog, "x-wall", wallGroup, clusterId++, 1.0, position, Q);
    position.set(xmax + Rw, 0.5*(ymin+ymax), 0.0);
    addParticle(fog, "x-wall", wallGroup, clusterId++, 1.0, position, Q);

    position.set(0.0, ymin - Rw, 0.0);
    addParticle(fog, "y-wall", wallGroup, clusterId++, 1.0, position, Q);
    position.set(0.0, ymax + Rw, 0.0);
    addParticle(fog, "y-wall", wallGroup, clusterId++, 1.0, position, Q);

    position.set(0.0, 0.5*(ymin+ymax), zmin - Rw);
    addParticle(fog, "z-wall", wallGroup, clusterId++, 1.0, position, Q);
    position.set(0.0, 0.5*(ymin+ymax), zmax + Rw);
    addParticle(fog, "z-wall", wallGroup, clusterId++, 1.0, position, Q);
#else
    fog << "Particles " << g.size() + 1 << std::endl;

    position.reset();
    Q.reset();
    addParticle(fog, "cylBox", wallGroup, clusterId++, 1.0, position, Q);
#endif

    int sphereGroup = 0;
    Q.reset();
    for (int i = 0; i < g.size(); i++) {
      position.set(g[i].x, g[i].y, g[i].z);
      addParticle(fog, "sphere", sphereGroup, clusterId++, g[i].R, position, Q);
    }
  }
}

bool intersect(sphere& si, sphere& sj) {
  double dx = sj.x - si.x;
  double dy = sj.y - si.y;
  double dz = sj.z - si.z;
  double distSqr = dx * dx + dy * dy + dz * dz;
  double sumR = si.R + sj.R;
  double sumRSqr = sumR * sumR;
  if (distSqr <= sumRSqr) return true;
  return false;
}

int main(int argc, char* argv[]) {

  vector<sphere> g;

#if defined(CONTAINER_BOX)
  std::cout << "CONTAINER = 6-WALLS-BOX\n";
  BoxContainer box;
#else
  std::cout << "CONTAINER = CYLINDER\n";
  CylContainer box;
#endif

#if defined(CONTAINER_BOX)
  box.L = baseLength;
#else
  box.diam = baseLength;
#endif

  // distribution uniforme
  cout << "generation of diameters...\n";
  double dD = (Dmax - Dmin) / (double)(nb_spheres);
  double Rtab[nb_spheres];
  double currentD = Dmax;
  double Vs = 0.0;
  for (size_t i = 0; i < nb_spheres; i++) {
    Rtab[i] = 0.5 * currentD;
    Vs += (4.0 / 3.0) * M_PI * Rtab[i] * Rtab[i] * Rtab[i];
    currentD -= dD;
  }

  g.clear();
  cout << "Targeted solid fraction = " << SolidFractionTarget << endl;
  double Vbox = Vs / SolidFractionTarget;

#if defined(CONTAINER_BOX)
  box.H = Vbox / (box.L * box.L);
#else
  box.H = Vbox / (M_PI * box.diam * box.diam / 4.0);
#endif

  default_random_engine generator(seed);
  uniform_real_distribution<double> random01(0.0, 1.0);
  uniform_real_distribution<double> random02PI(0.0, 2 * M_PI);
  uniform_real_distribution<double> random0PI(0.0, M_PI);
  double sumV = 0.0;
  size_t cumul_trials = 0;

#if defined(CONTAINER_BOX)
  sphere g1(0.0, 0.5 * box.H, 0.0, Rtab[0]);
#else
  sphere g1(0.0, 0.5 * box.H, 0.0, Rtab[0]);
#endif

  g.push_back(g1);
  sumV += (4.0 / 3.0) * M_PI * Rtab[0] * Rtab[0] * Rtab[0];
  cumul_trials++;

  cout << "Running...\n" << flush;
  for (size_t i = 1; i < nb_spheres; i++) {
    sphere gi;
    double Vi = (4.0 / 3.0) * M_PI * Rtab[i] * Rtab[i] * Rtab[i];

    size_t nb_trials = 1;
  retry:

#if defined(CONTAINER_BOX)
    double distBound = Rtab[i];
    double widthFill = box.L - 2.0 * distBound;
    double heightFill = box.H - 2.0 * distBound;

    gi.x = distBound + widthFill * random01(generator) - 0.5 * box.L;
    gi.y = distBound + heightFill * random01(generator);
    gi.z = distBound + widthFill * random01(generator) - 0.5 * box.L;
    gi.R = Rtab[i];
#else
    double distBound = Rtab[i];
    double Rmax = 0.5 * box.diam - distBound;
    double angle = random02PI(generator);
    double heightFill = box.H - 2.0 * distBound;
    gi.x = cos(angle) * Rmax * random01(generator);
    gi.y = distBound + heightFill * random01(generator);
    gi.z = sin(angle) * Rmax * random01(generator);
    gi.R = Rtab[i];
#endif

    for (long k = 0; k < g.size(); k++) {

      if (intersect(gi, g[k])) {
        nb_trials++;
        if (nb_trials >= nb_trials_max) {
          cout << "Maximum number of trials reached !" << endl;
          cout << "Number of shapes packed = " << i << "\n";
          cout << "Solid fraction = " << sumV / Vbox << "\n";
#if defined(CONTAINER_BOX)
          save_pack("pack.txt", g, -0.5 * box.L, 0, -0.5 * box.L, 0.5 * box.L, box.H, 0.5 * box.L);
#else
          save_pack("pack.txt", g, -0.5 * box.diam, 0, -0.5 * box.diam, 0.5 * box.diam, box.H, 0.5 * box.diam);
#endif
          return 0;
        }
        goto retry;
      }
    }
    g.push_back(gi);
    cumul_trials += nb_trials;
    sumV += Vi;
  }
  cout << "Okay ! The " << nb_spheres << " shapes have been packed\n";
  cout << "Solid fraction = " << sumV / Vbox << "\n";
#if defined(CONTAINER_BOX)
  save_pack("pack.txt", g, -0.5 * box.L, 0, -0.5 * box.L, 0.5 * box.L, box.H, 0.5 * box.L);
  save_rockable("input.txt", "shapes.txt", g, -0.5 * box.L, 0, -0.5 * box.L, 0.5 * box.L, box.H, 0.5 * box.L);
#else
  save_pack("pack.txt", g, -0.5 * box.diam, 0, -0.5 * box.diam, 0.5 * box.diam, box.H, 0.5 * box.diam);
  save_rockable("input.txt", "shapes.txt", g, -0.5 * box.diam, 0, -0.5 * box.diam, 0.5 * box.diam, box.H,
                0.5 * box.diam);
#endif
}
