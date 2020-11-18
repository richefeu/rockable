#include <cmath>
#include <fstream>
#include <iostream>

#include "generatePacking_wallBox.hpp"
#include "generateShape_xyz_walls.hpp"
#include "generateShape_cube.hpp"

#include "message.hpp"
#include "quat.hpp"

void generateSystem_triAxialCubes(std::ostream& os, const char* cubeName, int wallGroup, int group, int nx, int ny,
                                   int nz, int nbCubesPerBranch, double cubeSize, double Rw) {
  using namespace std;
  int Np = nx * ny * nz;  // number of clusters in fact

  os << "Particles " << Np * (nbCubesPerBranch * 6 + 1) + 6 << endl;
  quat Qclust;

  int clustID = 0;
  double dstep = sqrt(2.0) * (1 + 2 * nbCubesPerBranch) * cubeSize;
  vec3r orig(0.5 * dstep, 0.5 * dstep, 0.5 * dstep);

  double LX = dstep * nx;
  double LY = dstep * ny;
  double LZ = dstep * nz;
  generatePacking_wallBox(os, wallGroup, LX, LY, LZ, Rw);

  ofstream shpFile("cube.shp");
  generateShape_cube(shpFile, cubeName, 0.1 * cubeSize, cubeSize);
  generateShape_xyz_walls(shpFile, LX, LY, LZ, Rw);

  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      for (int iz = 0; iz < nz; iz++) {

        vec3r gridPos = orig + vec3r(ix * dstep, iy * dstep, iz * dstep);
        Qclust.randomize();

        msg::bestPrecision(os);
        os << cubeName << " " << group << " " << clustID << " 1 " << gridPos << "  0 0 0  0 0 0  " << Qclust
           << "  0 0 0  0 0 0 " << endl;
        for (int i = 0; i < nbCubesPerBranch; i++) {
          os << cubeName << " " << group << " " << clustID << " 1 "
             << gridPos + (i + 1) * cubeSize * (Qclust * vec3r::unit_x()) << "  0 0 0  0 0 0  " << Qclust
             << "  0 0 0  0 0 0 " << endl;
          os << cubeName << " " << group << " " << clustID << " 1 "
             << gridPos - (i + 1) * cubeSize * (Qclust * vec3r::unit_x()) << "  0 0 0  0 0 0  " << Qclust
             << "  0 0 0  0 0 0 " << endl;

          os << cubeName << " " << group << " " << clustID << " 1 "
             << gridPos + (i + 1) * cubeSize * (Qclust * vec3r::unit_y()) << "  0 0 0  0 0 0  " << Qclust
             << "  0 0 0  0 0 0 " << endl;
          os << cubeName << " " << group << " " << clustID << " 1 "
             << gridPos - (i + 1) * cubeSize * (Qclust * vec3r::unit_y()) << "  0 0 0  0 0 0  " << Qclust
             << "  0 0 0  0 0 0 " << endl;

          os << cubeName << " " << group << " " << clustID << " 1 "
             << gridPos + (i + 1) * cubeSize * (Qclust * vec3r::unit_z()) << "  0 0 0  0 0 0  " << Qclust
             << "  0 0 0  0 0 0 " << endl;
          os << cubeName << " " << group << " " << clustID << " 1 "
             << gridPos - (i + 1) * cubeSize * (Qclust * vec3r::unit_z()) << "  0 0 0  0 0 0  " << Qclust
             << "  0 0 0  0 0 0 " << endl;
        }
        clustID++;
        msg::normalPrecision(os);
      }
    }
  }
}


