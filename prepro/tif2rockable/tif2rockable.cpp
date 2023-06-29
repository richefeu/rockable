/*
g++-12 -o test tif2rockable.cpp -I/usr/X11R6/include
-I/usr/local/Cellar/libtiff/4.4.0_1/include
-L/usr/local/Cellar/libtiff/4.4.0_1/lib -L/usr/X11R6/lib -lX11 -lpthread -ltiff
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

#include "GrainShape.hpp"

#define cimg_use_tiff
#include "CImg-3.2.1/CImg.h"
using namespace cimg_library;

struct Center {
  double x, y, z;
  unsigned int N;
  bool isROI;
  Center() : x(0.0), y(0.0), N(0), isROI(false) {}
};

std::vector<Center> findCenters(CImg<> &img, size_t nbMax = 10000) {

  std::cout << "Compute..." << std::endl;
  std::vector<Center> centers(nbMax);
  for (int z = 0; z < img.depth(); z++) {
    for (int y = 0; y < img.height(); y++) {
      for (int x = 0; x < img.width(); x++) {
        size_t id = (size_t)(img(x, y, z));
        centers[id].x += (double)x;
        centers[id].y += (double)y;
        centers[id].z += (double)z;
        centers[id].N += 1;
      }
    }
  }

  int nbG = 0;
  for (size_t i = 0; i < centers.size(); i++) {
    if (centers[i].N > 0) {
      nbG++;
      double NN = (double)(centers[i].N);
      centers[i].x /= NN;
      centers[i].y /= NN;
      centers[i].z /= NN;
    }
  }

  std::cout << "Nombre de zones segmentees = " << nbG << std::endl;

  return centers;
}

void getROIs(const char *filename, std::vector<Center> &centers) {
  std::ifstream file(filename);

  size_t number;
  int found = 0;
  while (file) {
    file >> number;
    if (file.eof())
      break;
    if (number < centers.size()) {
      centers[number].isROI = true;
      found++;
    }
  }

  std::cout << "found = " << found << std::endl;
}

void wrap(CImg<> &img, GrainShape &G) {

  double Rmax = std::max({img.width(), img.height(), img.depth()});
  Rmax /= 2;

  for (size_t d = 0; d < G.points.size(); d++) {
    double R = 0.0;
    vec3r n = G.points[d] - G.pos;
    n.normalize();
    int x, y, z;
    vec3r prevPos = G.pos;

    while (R < Rmax) {
      R += 1.0;
      x = (int)round(G.pos.x + R * n.x);
      y = (int)round(G.pos.y + R * n.y);
      z = (int)round(G.pos.z + R * n.z);

      // il faudra checker qu'on tape pas en dehors de l'image
      vec3r boundCorrection;
      if (x < 0) {
        G.points[d] = G.pos + R * n;
        boundCorrection = G.points[d];
        G.points[d].x = 0.0;
        boundCorrection -= G.points[d];
        x = 0;
      } else if (x >= (int)img.width()) {
        G.points[d] = G.pos + R * n;
        boundCorrection = G.points[d];
        G.points[d].x = (double)img.width() - 1.0;
        boundCorrection -= G.points[d];
        x = img.width() - 1;
      }
      if (y < 0) {
        G.points[d] = G.pos + R * n;
        boundCorrection = G.points[d];
        G.points[d].y = 0.0;
        boundCorrection -= G.points[d];
        y = 0;
      } else if (y >= (int)img.height()) {
        G.points[d] = G.pos + R * n;
        boundCorrection = G.points[d];
        G.points[d].y = (double)img.height() - 1.0;
        y = img.height() - 1;
        boundCorrection -= G.points[d];
      }
      if (z < 0) {
        G.points[d] = G.pos + R * n;
        boundCorrection = G.points[d];
        G.points[d].z = 0.0;
        boundCorrection -= G.points[d];
        z = 0;
      } else if (z >= (int)img.depth()) {
        G.points[d] = G.pos + R * n;
        boundCorrection = G.points[d];
        G.points[d].z = (double)img.depth() - 1.0;
        z = img.depth() - 1;
        boundCorrection -= G.points[d];
      }

      int col = (int)(img(x, y, z));
      if (G.colorId != col) {
        G.points[d] = prevPos;
        break;
      }
      prevPos = G.pos + R * n - boundCorrection;
    }
  }
}

int main(int argc, char const *argv[]) {
  CImg<> img("Data/LENGP00/02_Results/01_Label/LENGP00_Label/LENGP00_00.tif");

  std::vector<Center> centers = findCenters(img, 10000);
  // stgetROIs("Data/LENGP00/02_Results/03_ROI/LENGP00_LabelsInside.txt",
  // centers);

  std::vector<GrainShape> grains;
  for (size_t i = 0 + 1; i < centers.size(); i++) {
    if (centers[i].N > 0) {
      GrainShape G(1); // 0 = icosaèdre, 1 = 4x plus de facettes, 2 = x16... 
      G.colorId = (int)i;
      G.pos.x = round(centers[i].x);
      G.pos.y = round(centers[i].y);
      G.pos.z = round(centers[i].z);
      G.setRadius(1);
      G.move(G.pos.x, G.pos.y, G.pos.z);
      grains.push_back(G);
    }
  }

  // domaine pour les grains exportés
	vec3r minBox(200, 200, 335);
	vec3r maxBox(265, 265, 400);

  std::ofstream sfile("shapes.txt");
  for (size_t i = 0; i < grains.size(); i++) {
    //if (grains[i].pos.z < 275 || grains[i].pos.z > 325)
    //  continue;
		if ( grains[i].pos.x < minBox.x || grains[i].pos.x > maxBox.x
			|| grains[i].pos.y < minBox.y || grains[i].pos.y > maxBox.y 
			|| grains[i].pos.z < minBox.z || grains[i].pos.z > maxBox.z ) continue;	
		
    std::cout << "Building Grains " << grains[i].colorId << std::endl;
    wrap(img, grains[i]);
    grains[i].R = 0.5;
    grains[i].saveRockableShape(sfile, grains[i].colorId, 180e-6);
  }

  /*
  Procédure dans shapeSurvey:
  ouvrir le fichier shapes.txt
  appuyer sur 'C' pour pré-calculer tous les propriétés de masse (c'est long)
  appuyer sur 's' pour sauvegarder les formes avec les valeurs pré-calculées
  (nom = mod_shapes.txt) appuyer sur 'p' pour sauvegarder les positions au
  format rockable (à insérer dans le fichier input.txt)
  */

  return 0;
}