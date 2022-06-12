#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <algorithm>

#include "addParticle.hpp"
#include "consoleProgressBar.hpp"
#include "generateShape_cylindricalMold.hpp"
#include "generateShape_diskPlate.hpp"
#include "generateShape_stick.hpp"
#include "generateShape_xyz_walls.hpp"

const int CONTAINER_BOX = 0;
const int CONTAINER_CYL = 1;

const int DISTRI_MODE_UNIF = 0;
const int DISTRI_MODE_TAB = 1;

const int SORT_NEVER = 0;
const int SORT_BEFORE = 1;
const int SORT_AFTER = 2;

using namespace std;

unsigned int seed = 3503; //  3503 143502  29292983  26

// =========================================================
// =========================================================
int containerType = CONTAINER_CYL;
int nSecCyl = 18;                  // number of sectors of the cylindrical container
double baseLength = 7e-2;        // diameter of cylindrical container or squared baselength of rectangular container
double R = 1e-3;                 // radius of the sticks
int distriMode = DISTRI_MODE_UNIF;  // 0: uniform in [Lmin Lmax]; 1: length defined in LDefTab[] and Lrate[]

// for uniformly distributed lengths
double Lmax = 2.5e-2;
double Lmin = 2.5e-2;

// for tabulated lengths
int nbLen = 1;                // number of different lengths
double LDefTab[] = {1.5e-3};  // lengths (largest first)
int Lrate[] = {100};          // percentages

size_t nb_tubes = 1000;
size_t nb_trials_max = nb_tubes * 100;

double SolidFractionTarget = 0.08;  // the algorithm will do its best for obtaining this solid fraction
int sortBottomTop = SORT_AFTER;
int nbPostCompact = 0;
// =========================================================
// =========================================================

struct tube {
  double x0, y0, z0;
  double x1, y1, z1;
  double R;
  int shapeId;
  tube() { this->set(0, 0, 0, 0, 0, 0, 0); }
  tube(double x, double y, double z, double anglez, double angley, double len, double rad) {
    this->set(x, y, z, anglez, angley, len, rad);
  }

  void set(double x, double y, double z, double anglez, double angley, double len, double rad) {
    double cy = cos(angley);
    double sy = sin(angley);
    double cz = cos(anglez);
    double sz = sin(anglez);
    double ux = cy * cz;
    double uy = sz;
    double uz = -cz * sy;

    double halfLen = 0.5 * len;
    x0 = x - halfLen * ux;
    y0 = y - halfLen * uy;
    z0 = z - halfLen * uz;
    x1 = x + halfLen * ux;
    y1 = y + halfLen * uy;
    z1 = z + halfLen * uz;
    R = rad;
  }
};

struct BoxContainer {
  double L;
  double H;
};

struct CylContainer {
  double diam;
  double H;
};

void save_pack(const char* fname, vector<tube>& g, double xmin, double ymin, double zmin, double xmax, double ymax,
               double zmax) {
  ofstream fog(fname, ios::out);
  if (fog) {
    fog.precision(5);
    fog << scientific;
    //fog << xmin << ' ' << ymin << ' ' << zmin << "\n";
    //fog << xmax << ' ' << ymax << ' ' << zmax << "\n";
    for (int i = 0; i < g.size(); i++) {
      fog << g[i].x0 << " " << g[i].y0 << " " << g[i].z0 << " " << g[i].x1 << " " << g[i].y1 << " " << g[i].z1 << " "
          << g[i].R << endl;
    }
  }
}

void save_rockable(const char* fInputName, const char* fShapesName, vector<tube>& g, double xmin, double ymin,
                   double zmin, double xmax, double ymax, double zmax) {

  double Rw = g[0].R;

  ofstream fshape(fShapesName, ios::out);
  fshape << setprecision(15);
  if (fshape) {

    if (containerType == CONTAINER_BOX) {
      generateShape_xyz_walls(fshape, 1.1 * (xmax - xmin), 1.1 * (ymax - ymin), 1.1 * (zmax - zmin), Rw);
    } else if (containerType == CONTAINER_CYL) {
      double RCylIn = 0.5 * (xmax - xmin);
      double beta = M_PI / (double)nSecCyl;
      double RCylInPlus = RCylIn / cos(beta);
      generateShape_cylindricalMold(fshape, "cylBox", Rw, ymax - ymin, nSecCyl, RCylInPlus);
      generateShape_diskPlate(fshape, "plate", Rw, nSecCyl, RCylInPlus);
    }

    char partName[256];
    if (distriMode == DISTRI_MODE_UNIF) {  // one shape per particle
      for (int i = 0; i < g.size(); i++) {
        double Lx = g[i].x1 - g[i].x0;
        double Ly = g[i].y1 - g[i].y0;
        double Lz = g[i].z1 - g[i].z0;
        double L = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
        sprintf(partName, "stick%d", g[i].shapeId);
        generateShape_stick(fshape, partName, L, g[i].R);
      }
    } else if (distriMode == DISTRI_MODE_TAB) {  // per class
      for (int i = 0; i < nbLen; i++) {
        sprintf(partName, "stick%d", i);
        generateShape_stick(fshape, partName, LDefTab[i], R);
      }
    }
  }

  ofstream fog(fInputName, ios::out);
  fog << setprecision(15);
  if (fog) {

    vec3r position;
    quat Q;
    int wallGroup = 1;
    int clusterId = 0;

    if (containerType == CONTAINER_BOX) {
      fog << "Particles " << g.size() + 6 << std::endl;
      Q.reset();
      position.set(xmin - Rw, 0.5 * (ymin + ymax), 0.0);
      addParticle(fog, "x-wall", wallGroup, clusterId++, 1.0, position, Q);
      position.set(xmax + Rw, 0.5 * (ymin + ymax), 0.0);
      addParticle(fog, "x-wall", wallGroup, clusterId++, 1.0, position, Q);

      position.set(0.0, ymin - Rw, 0.0);
      addParticle(fog, "y-wall", wallGroup, clusterId++, 1.0, position, Q);
      position.set(0.0, ymax + Rw, 0.0);
      addParticle(fog, "y-wall", wallGroup, clusterId++, 1.0, position, Q);

      position.set(0.0, 0.5 * (ymin + ymax), zmin - Rw);
      addParticle(fog, "z-wall", wallGroup, clusterId++, 1.0, position, Q);
      position.set(0.0, 0.5 * (ymin + ymax), zmax + Rw);
      addParticle(fog, "z-wall", wallGroup, clusterId++, 1.0, position, Q);
    } else if (containerType == CONTAINER_CYL) {
      fog << "Particles " << g.size() + 2 << std::endl;

      position.reset();
      Q.reset();
      addParticle(fog, "cylBox", wallGroup, clusterId++, 1.0, position, Q);
      position.set(0.0, ymax + Rw, 0.0);
      addParticle(fog, "plate", wallGroup, clusterId++, 1.0, position, Q);
    }

    vec3r xAxis(1.0, 0.0, 0.0);
    vec3r tubeAxis;
    int stickGroup = 0;
    char partName[256];
    for (int i = 0; i < g.size(); i++) {
      position.set(0.5 * (g[i].x0 + g[i].x1), 0.5 * (g[i].y0 + g[i].y1), 0.5 * (g[i].z0 + g[i].z1));
      tubeAxis.set(g[i].x1 - g[i].x0, g[i].y1 - g[i].y0, g[i].z1 - g[i].z0);
      tubeAxis.normalize();
      Q.set_from_to(xAxis, tubeAxis);
      sprintf(partName, "stick%d", g[i].shapeId);
      addParticle(fog, partName, stickGroup, clusterId++, 1.0, position, Q);
    }
  }
}

#define _EPSILON_VALUE_ 1.0e-12
bool intersectCylCyl(tube& ti, tube& tj) {

  double Eix = ti.x1 - ti.x0;
  double Eiy = ti.y1 - ti.y0;
  double Eiz = ti.z1 - ti.z0;

  double Ejx = tj.x1 - tj.x0;
  double Ejy = tj.y1 - tj.y0;
  double Ejz = tj.z1 - tj.z0;

  double vx = ti.x0 - tj.x0;
  double vy = ti.y0 - tj.y0;
  double vz = ti.z0 - tj.z0;

  double c = Eix * Eix + Eiy * Eiy + Eiz * Eiz;
  double d = Ejx * Ejx + Ejy * Ejy + Ejz * Ejz;
  double e = Eix * Ejx + Eiy * Ejy + Eiz * Ejz;
  double f = (c * d) - (e * e);
  double s, t;

  double lx, ly, lz;
  if (fabs(f) > _EPSILON_VALUE_) {
    f = 1.0 / f;
    double a = Eix * vx + Eiy * vy + Eiz * vz;
    double b = Ejx * vx + Ejy * vy + Ejz * vz;
    s = (e * b - a * d) * f;  // for stick i
    t = (c * b - e * a) * f;  // for stick j

    if (s < 0.0) s = 0.0;
    if (s > 1.0) s = 1.0;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    lx = (ti.x0 + s * Eix) - (tj.x0 + t * Ejx);
    ly = (ti.y0 + s * Eiy) - (tj.y0 + t * Ejy);
    lz = (ti.z0 + s * Eiz) - (tj.z0 + t * Ejz);
  } else {         // f = 0 means that the sticks are parallel
    return false;  // it will be managed by other tests
  }

  // from here lx, ly, lz have been computed
  double l2 = lx * lx + ly * ly + lz * lz;
  double sum = ti.R + tj.R;
  if (l2 < sum * sum) return true;

  return false;
}
#undef _EPSILON_VALUE_

bool intersectSphCyl(tube& t, double x, double y, double z, double R) {
  double Ex = t.x1 - t.x0;
  double Ey = t.y1 - t.y0;
  double Ez = t.z1 - t.z0;

  double vx = x - t.x0;
  double vy = y - t.y0;
  double vz = z - t.z0;

  double r = (vx * Ex + vy * Ey + vz * Ez) / (Ex * Ex + Ey * Ey + Ez * Ez);

  if (r < 0.0) r = 0.0;
  if (r > 1.0) r = 1.0;

  double lx = t.x0 + r * Ex - x;
  double ly = t.y0 + r * Ey - y;
  double lz = t.z0 + r * Ez - z;

  double l2 = lx * lx + ly * ly + lz * lz;
  double sum = t.R + R;

  if (l2 < sum * sum) return true;

  return false;
}

bool intersectSphSph(double x0, double y0, double z0, double R0, double x1, double y1, double z1, double R1) {
  double lx = x0 - x1;
  double ly = y0 - y1;
  double lz = z0 - z1;

  double l2 = lx * lx + ly * ly + lz * lz;
  double sum = R0 + R1;

  if (l2 < sum * sum) return true;

  return false;
}

bool intersect(tube& ti, tube& tj) {
  if (intersectCylCyl(ti, tj)) return true;
  if (intersectSphCyl(ti, tj.x0, tj.y0, tj.z0, tj.R)) return true;
  if (intersectSphCyl(ti, tj.x1, tj.y1, tj.z1, tj.R)) return true;
  if (intersectSphCyl(tj, ti.x0, ti.y0, ti.z0, ti.R)) return true;
  if (intersectSphCyl(tj, ti.x1, ti.y1, ti.z1, ti.R)) return true;
  if (intersectSphSph(ti.x0, ti.y0, ti.z0, ti.R, tj.x0, tj.y0, tj.z0, tj.R)) return true;
  if (intersectSphSph(ti.x0, ti.y0, ti.z0, ti.R, tj.x1, tj.y1, tj.z1, tj.R)) return true;
  if (intersectSphSph(ti.x1, ti.y1, ti.z1, ti.R, tj.x0, tj.y0, tj.z0, tj.R)) return true;
  if (intersectSphSph(ti.x1, ti.y1, ti.z1, ti.R, tj.x1, tj.y1, tj.z1, tj.R)) return true;
  return false;
}

// ======================================================================
// ======================================================================
int main(int argc, char* argv[]) {
  vector<tube> g;

  BoxContainer box;
  CylContainer cyl;
  if (containerType == CONTAINER_BOX) {
    std::cout << "CONTAINER = 6-WALLS-BOX\n";
    box.L = baseLength;
  } else if (containerType == CONTAINER_CYL) {
    std::cout << "CONTAINER = CYLINDER\n";
    cyl.diam = baseLength;
  }

  double Ltab[nb_tubes];
  int shapeIdtab[nb_tubes];
  double Vs = 0.0;
  if (distriMode == DISTRI_MODE_UNIF) {  // distribution uniforme
    cout << "generation of uniformly distributed lengths...\n";
    double dL = (Lmax - Lmin) / (double)(nb_tubes);

    double currentL = Lmax;
    Vs = 0.0;
    for (size_t i = 0; i < nb_tubes; i++) {
      Ltab[i] = currentL;
      shapeIdtab[i] = i;
      Vs += (4.0 / 3.0) * M_PI * R * R * R + M_PI * R * R * currentL;
      currentL -= dL;
    }
  } else if (distriMode == DISTRI_MODE_TAB) {
    cout << "generation of tabulated (per class) lengths...\n";
    size_t i = 0;
    int nbAdded = 0;
    for (int shapeId = 0; shapeId < nbLen; shapeId++) {
      double currentL = LDefTab[shapeId];
      if (shapeId == nbLen - 1) {  // last
        nbAdded = nb_tubes - i;
      } else {
        nbAdded = (int)floor(0.01 * Lrate[shapeId] * nb_tubes);
      }
      Vs = 0.0;
      for (int a = 0; a < nbAdded; a++) {
        Ltab[i] = currentL;
        shapeIdtab[i] = shapeId;
        Vs += (4.0 / 3.0) * M_PI * R * R * R + M_PI * R * R * currentL;
        i++;
      }
    }
  }

  g.clear();
  cout << "Targeted solid fraction = " << SolidFractionTarget << endl;
  double Vbox = Vs / SolidFractionTarget;

  if (containerType == CONTAINER_BOX) {
    box.H = Vbox / (box.L * box.L);
  } else if (containerType == CONTAINER_CYL) {
    cyl.H = Vbox / (M_PI * cyl.diam * cyl.diam / 4.0);
  }

  default_random_engine generator(seed);
  uniform_real_distribution<double> random01(0.0, 1.0);
  uniform_real_distribution<double> random02PI(0.0, 2.0 * M_PI);
  uniform_real_distribution<double> random0PI(0.0, M_PI);
  double sumV = 0.0;
  size_t cumul_trials = 0;

  tube g1;
  if (containerType == CONTAINER_BOX) {
    g1.set(0.0, 0.5 * box.H, 0.0, random02PI(generator), random0PI(generator), Ltab[0], R);
  } else if (containerType == CONTAINER_CYL) {
    g1.set(0.0, 0.5 * cyl.H, 0.0, random02PI(generator), random0PI(generator), Ltab[0], R);
  }

  g1.shapeId = shapeIdtab[0];
  g.push_back(g1);
  sumV += (4.0 / 3.0) * M_PI * R * R * R + M_PI * R * R * Ltab[0];
  cumul_trials++;

  ConsoleProgressBar progressBar(nb_tubes - 1);
  progressBar.setTitle("Placement: ");
  progressBar.update(1, std::cerr);

  for (size_t i = 1; i < nb_tubes; i++) {
    tube gi;
    double Vi = (4.0 / 3.0) * M_PI * R * R * R + M_PI * R * R * Ltab[i];

    size_t nb_trials = 1;
  retry:

    if (containerType == CONTAINER_BOX) {
      double distBound = 0.5 * Ltab[i] + R;
      double widthFill = box.L - 2.0 * distBound;
      double heightFill = box.H - 2.0 * distBound;
      gi.set(distBound + widthFill * random01(generator) - 0.5 * box.L, distBound + heightFill * random01(generator),
             distBound + widthFill * random01(generator) - 0.5 * box.L, random02PI(generator), random0PI(generator),
             Ltab[i], R);
      gi.shapeId = shapeIdtab[i];
    } else if (containerType == CONTAINER_CYL) {
      double distBound = 0.5 * Ltab[i] + R;
      double Rmax = 0.5 * cyl.diam - distBound;
      double angle = random02PI(generator);
      double heightFill = cyl.H - 2.0 * distBound;
      gi.set(cos(angle) * Rmax * random01(generator), distBound + heightFill * random01(generator),
             sin(angle) * Rmax * random01(generator), random02PI(generator), random0PI(generator), Ltab[i], R);
      gi.shapeId = shapeIdtab[i];
    }

    bool inters = false;
    for (long k = 0; k < g.size(); k++) {

      if (intersect(gi, g[k])) {
        nb_trials++;
        if (nb_trials >= nb_trials_max) {
          cout << "\nMaximum number of trials reached !" << endl;
          cout << "Number of shapes packed = " << g.size() << "\n";
          cout << "Solid fraction = " << sumV / Vbox << "\n";
          cout << "input files for Rockable have NOT been created (try with a lower packing fraction target)\n";
          if (containerType == CONTAINER_BOX) {
            save_pack("pack.txt", g, -0.5 * box.L, 0, -0.5 * box.L, 0.5 * box.L, box.H, 0.5 * box.L);
          } else if (containerType == CONTAINER_CYL) {
            save_pack("pack.txt", g, -0.5 * cyl.diam, 0, -0.5 * cyl.diam, 0.5 * cyl.diam, cyl.H, 0.5 * cyl.diam);
          }
          return 0;
        }

        inters = true;
        break;
      }
    }

    if (inters == true) goto retry;

    g.push_back(gi);
    progressBar.update(i, std::cerr);
    cumul_trials += nb_trials;
    sumV += Vi;
  }

  cout << "\nOkay ! The " << nb_tubes << " particles have been packed\n";
  cout << "Solid fraction = " << sumV / Vbox << "\n";

  if (sortBottomTop == SORT_BEFORE) {
    sort(g.begin(), g.end(), [](const tube& a, const tube& b) -> bool {
      double aymin = std::min(a.y0, a.y1);
      double bymin = std::min(b.y0, b.y1);
      return (aymin < bymin);
    });
    cout << "\nBodies ordered from bottom to top.\n" << flush;
  } else if (sortBottomTop == SORT_NEVER) {
    cout << "No ordering of the bodies from bottom to top.\n" << flush;
  }

  if (nbPostCompact > 0) {
	  progressBar.setTitle("pre-compaction: ");
    progressBar.setMax(nbPostCompact - 1);
    double dy = R * 0.01;
    for (int c = 0; c < nbPostCompact; c++) {
      for (size_t i = 0; i < g.size(); i++) {
        bool ok = true;
        while (ok) {
          g[i].y0 -= dy;
          g[i].y1 -= dy;
          if (g[i].y0 < R || g[i].y1 < R) {
            ok = false;
          } else {
            for (size_t j = 0; j < g.size(); j++) {
              if (j == i) continue;
              if (intersect(g[i], g[j])) {
                ok = false;
                break;
              }
            }
          }
        }
        g[i].y0 += dy;
        g[i].y1 += dy;
      }
      progressBar.update(c, std::cerr);
    }
  }

  if (sortBottomTop == SORT_AFTER) {
    sort(g.begin(), g.end(), [](const tube& a, const tube& b) -> bool {
      double aymin = std::min(a.y0, a.y1);
      double bymin = std::min(b.y0, b.y1);
      return (aymin < bymin);
    });
    cout << "\nBodies ordered from bottom to top.\n" << flush;
  }

  if (containerType == CONTAINER_BOX) {
    save_pack("pack.txt", g, -0.5 * box.L, 0, -0.5 * box.L, 0.5 * box.L, box.H, 0.5 * box.L);
    save_rockable("input.txt", "shapes.txt", g, -0.5 * box.L, 0.0, -0.5 * box.L, 0.5 * box.L, box.H, 0.5 * box.L);
  } else if (containerType == CONTAINER_CYL) {
    save_pack("pack.txt", g, -0.5 * cyl.diam, 0, -0.5 * cyl.diam, 0.5 * cyl.diam, cyl.H, 0.5 * cyl.diam);
    save_rockable("input.txt", "shapes.txt", g, -0.5 * cyl.diam, 0.0, -0.5 * cyl.diam, 0.5 * cyl.diam, cyl.H,
                  0.5 * cyl.diam);
  }
  cout << "Files 'input.txt' and 'shapes.txt' have been created.\n";
}
