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

#ifndef GENERATESYSTEM_TRIAXIALCUBES_HPP
#define GENERATESYSTEM_TRIAXIALCUBES_HPP

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

#endif /* end of include guard: GENERATESYSTEM_TRIAXIALCUBES_HPP */ 
