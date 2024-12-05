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
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <locale>
#include <string>

#include "OBBpacker.hpp"
#include "kwParser.hpp"
#include "message.hpp"
#include "nanoExprParser.hpp"
#include "quat.hpp"
#include "transformation.hpp"
#include "vec3.hpp"

Transformation<double> globalTransformation;
quat individualParticleRotation;

int NumberOfParticles = 0;
std::streampos Np_pos = std::streampos(0);

#include "addParticle.hpp"
#include "generatePacking_grid.hpp"
#include "generatePacking_wallBox.hpp"
#include "generateShape_cube.hpp"
#include "generateShape_cuboid.hpp"
#include "generateShape_cuboid_container.hpp"
#include "generateShape_pyramid3.hpp"
#include "generateShape_rectangle.hpp"
#include "generateShape_rhombicuboctahedron.hpp"
#include "generateShape_sphere.hpp"
#include "generateShape_thinCylinder.hpp"
#include "generateShape_xyz_walls.hpp"

std::ostream* outputStream = &std::cout;
std::ofstream fileStream;

void inplace_trim(std::string& line) {
  // trim from start (in place)
  line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) { return !std::isspace(ch); }));
  // trim from end (in place)
  line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(),
             line.end());
}

// This is a very basic kind of language to generate rockable initial conf-file
void readCommands(const char* name) {
  std::ifstream com(name);
  if (!com.is_open()) {
    std::cerr << msg::warn() << "@readCommands, Cannot read " << name << msg::normal() << std::endl;
    return;
  }

  kwParser parser;

  parser.kwMap["open"] = __DO__(com) {
    std::string filename;
    com >> filename;
    fileStream.open(filename);
    if (!fileStream.is_open()) {
      std::cerr << "Cannot open output file: " << filename << std::endl;
    } else {
      outputStream = &fileStream;
    }
  };

  parser.kwMap["close"] = __DO__(com) {
    if (fileStream.is_open()) {
      if (Np_pos != std::streampos(0)) {
        fileStream.seekp(Np_pos);
        fileStream << NumberOfParticles;
      }
      fileStream.close();
      outputStream = &std::cout;
    }
  };

  parser.kwMap["print"] = __DO__(com) {
    std::string line;
    getline(com, line);
    inplace_trim(line);
    *outputStream << line << '\n';
  };

  parser.kwMap["print>"] = __DO__(com) {
    std::string line;
    getline(com, line);
    inplace_trim(line);
    *outputStream << line;
  };

  parser.kwMap["compute"] = __DO__(com) {
    std::string preword;
    com >> preword;

    std::string line;
    getline(com, line);

    nanoExprParser<double> parser;
    double result;
    if (parser.parse(line, result)) {
      *outputStream << preword << ' ' << result << '\n';
    } else {
      std::cerr << "Error parsing expression: " << line << std::endl;
      *outputStream << preword << " ???_parse_error_???\n";
    }
  };

  parser.kwMap["<compute"] = __DO__(com) {
    std::string line;
    getline(com, line);

    nanoExprParser<double> parser;
    double result;
    if (parser.parse(line, result)) {
      *outputStream << " " << result << '\n';
    } else {
      std::cerr << "Error parsing expression: " << line << std::endl;
      *outputStream << " ???_parse_error_???\n";
    }
  };

  parser.kwMap["incrementNoParticles"] = __DO__(com) {
    int inc;
    com >> inc;
    NumberOfParticles += inc;
  };

  parser.kwMap["Particles:placeHolder"] = __DO__() {
    if (fileStream.is_open()) {
      fileStream << "Particles ";
      Np_pos = fileStream.tellp();
      fileStream << "            \n";
    }
  };

  parser.kwMap["tranformation:add_translation"] = __DO__(com) {
    double xtrans, ytrans, ztrans;
    com >> xtrans >> ytrans >> ztrans;
    globalTransformation.translate(xtrans, ytrans, ztrans);
  };

  parser.kwMap["tranformation:add_rotation"] = __DO__(com) {
    double xaxe, yaxe, zaxe, angleDeg;
    com >> xaxe >> yaxe >> zaxe >> angleDeg;
    double angle = angleDeg * M_PI / 180.0;
    globalTransformation.rotate(xaxe, yaxe, zaxe, angle);
    vec3r axis(xaxe, yaxe, zaxe);
    axis.normalize();
    individualParticleRotation.set_axis_angle(axis, angle);
  };

  parser.kwMap["tranformation:reset"] = __DO__() {
    globalTransformation.reset();
    individualParticleRotation.reset();
  };

  parser.kwMap["tranformation:reset_translation"] = __DO__() { globalTransformation.reset_translation(); };

  parser.kwMap["tranformation:reset_rotation"] = __DO__() {
    globalTransformation.reset_rotation();
    individualParticleRotation.reset();
  };

  parser.kwMap["addParticle"] = __DO__(com) {
    std::string name;
    int group;
    int cluster;
    double homothety;
    vec3r position;
    quat angularPosition;
    com >> name >> group >> cluster >> homothety >> position >> angularPosition;
    addParticle(*outputStream, name.c_str(), group, cluster, homothety, position, angularPosition);
    NumberOfParticles++;
  };

  parser.kwMap["generateShape:sphere"] = __DO__(com) {
    std::string name;
    double radius;
    com >> name >> radius;
    generateShape_sphere(*outputStream, name.c_str(), radius);
    std::cerr << "A spherical shape has been generated\n";
  };

  parser.kwMap["generateShape:cube"] = __DO__(com) {
    std::string name;
    double radius;
    double sideSize;
    com >> name >> radius >> sideSize;
    generateShape_cube(*outputStream, name.c_str(), radius, sideSize);
    std::cerr << "A cubic shape has been generated\n";
  };

  parser.kwMap["generateShape:cuboid"] = __DO__(com) {
    std::string name;
    double radius;
    vec3r sideSize;
    com >> name >> radius >> sideSize;
    generateShape_cuboid(*outputStream, name.c_str(), radius, sideSize);
    std::cerr << "A cuboid shape has been generated\n";
  };

  parser.kwMap["generateShape:pyramid3"] = __DO__(com) {
    std::string name;
    double radius;
    double sideSize;
    com >> name >> radius >> sideSize;
    generateShape_pyramid3(*outputStream, name.c_str(), radius, sideSize);
    std::cerr << "A pyramid3 shape has been generated\n";
  };

  parser.kwMap["generateShape:thin_cylinder"] = __DO__(com) {
    std::string name;
    double Rin, Rout, H;
    int nbSectors;
    com >> name >> Rin >> Rout >> H >> nbSectors;
    generateShape_thinCylinder(*outputStream, name.c_str(), Rin, Rout, H, nbSectors);
    std::cerr << "A thin cylinder shape has been generated\n";
  };

  parser.kwMap["generateShape:cuboid_container"] = __DO__(com) {
    std::string name;
    double radius;
    vec3r sideSize;
    int hasTop = 1, hasBottom = 1;
    com >> name >> radius >> sideSize >> hasTop >> hasBottom;
    bool hasTopBool = (hasTop == 0) ? false : true;
    bool hasBottomBool = (hasBottom == 0) ? false : true;
    generateShape_cuboid_container(*outputStream, name.c_str(), radius, sideSize, hasTopBool, hasBottomBool);
    std::cerr << "A cuboid container shape has been generated\n";
  };

  parser.kwMap["generateShape:rectangle_xz"] = __DO__(com) {
    std::string name;
    double radius;
    double sidex;
    double sidez;
    com >> name >> radius >> sidex >> sidez;
    generateShape_rectangle_xz(*outputStream, name.c_str(), radius, sidex, sidez);
    std::cerr << "A rectangle shape (in plane xz) has been generated\n";
  };

  parser.kwMap["generateShape:rectangle_xy"] = __DO__(com) {
    std::string name;
    double radius;
    double sidex;
    double sidey;
    com >> name >> radius >> sidex >> sidey;
    generateShape_rectangle_xy(*outputStream, name.c_str(), radius, sidex, sidey);
    std::cerr << "A rectangle shape (in plane xy) has been generated\n";
  };

  parser.kwMap["generateShape:rhombicuboctahedron"] = __DO__(com) {
    std::string name;
    double radius;
    vec3r sideSize;
    com >> name >> radius >> sideSize;
    generateShape_rhombicuboctahedron(*outputStream, name.c_str(), radius, sideSize);
    std::cerr << "A rhombicuboctahedron shape has been generated\n";
  };

  parser.kwMap["generateShape:xyz_walls"] = __DO__(com) {
    vec3r size;
    double Rw;
    com >> size >> Rw;
    generateShape_xyz_walls(*outputStream, size, Rw);
    std::cerr << "3 wall-shapes have been generated\n";
  };

  parser.kwMap["generatePacking:wallBox"] = __DO__(com) {
    int group;
    double LX, LY, LZ;
    double Rw;
    com >> group >> LX >> LY >> LZ >> Rw;
    int nbAdded = generatePacking_wallBox(*outputStream, group, LX, LY, LZ, Rw);
    NumberOfParticles += nbAdded;
    std::cerr << "@generatePacking:wallBox, Number of added bodies: " << nbAdded << '\n';
  };

  parser.kwMap["generatePacking:grid"] = __DO__(com) {
    std::string name;
    vec3r origBox;
    vec3r boxSize;
    vec3i n;
    int group;
    int clust;
    double homothety;
    int randQ;
    com >> name >> origBox >> boxSize >> n >> group >> clust >> homothety >> randQ;
    int nbAdded =
        generatePacking_grid(*outputStream, name.c_str(), origBox, boxSize, n, group, clust, homothety, randQ);
    NumberOfParticles += nbAdded;
    std::cerr << "@generatePacking:grid, Number of added bodies: " << nbAdded << '\n';
  };

  parser.kwMap["generatePacking:grid_clust"] = __DO__(com) {
    std::string name;
    vec3r origBox;
    vec3r boxSize;
    vec3i n;
    int group;
    int clust;
    double homothety;
    int randQ;
    com >> name >> origBox >> boxSize >> n >> group >> clust >> homothety >> randQ;
    int nbAdded =
        generatePacking_grid_clust(*outputStream, name.c_str(), origBox, boxSize, n, group, clust, homothety, randQ);
    NumberOfParticles += nbAdded;
    std::cerr << "@generatePacking:grid_clust, Number of added bodies: " << nbAdded << '\n';
  };

  parser.kwMap["generatePacking:RandomClosePacking"] = __DO__(com) {
    std::string name;
    std::string boxShape;
    std::string direction;
    vec3r boxSize;
    double xOBB;
    double yOBB;
    double zOBB;
    double homothety_min;
    double homothety_max;
    size_t nbTarget;
    double solidFractionTarget;
    int group;
    int clust;

    com >> name >> boxShape >> direction >> boxSize >> xOBB >> yOBB >> zOBB >> homothety_min >> homothety_max >>
        nbTarget >> solidFractionTarget >> group >> clust;

    OBBPacker packer;
    packer.setOBBDimensions(xOBB, yOBB, zOBB);
    packer.setHomothetyRange(homothety_min, homothety_max);
    packer.setNbOBBsTarget(nbTarget);
    packer.setSolidFractionTarget(solidFractionTarget);
    packer.setContainerLength(boxSize.x, boxSize.y, boxSize.z);
    if (boxShape == "CYLINDER") {
      packer.setContainerType(OBBPacker::CYLINDER);
    } else {
      packer.setContainerType(OBBPacker::CUBOID);
    }
    if (direction == "X") {
      packer.setPackingDirection(OBBPacker::X);
    } else if (direction == "Y") {
      packer.setPackingDirection(OBBPacker::Y);
    } else if (direction == "Z") {
      packer.setPackingDirection(OBBPacker::Z);
    }

    std::vector<OBB> obbs;
    packer.pack(obbs);
    packer.sort(obbs);

    for (size_t i = 0; i < obbs.size(); i++) {
      double homothety = 2.0 * obbs[i].extent.x / xOBB;
      quat angularPosition = obbs[i].getQuaternion();
      vec3r position = obbs[i].center;  // Ok only if the particle center is the obb center
      addParticle(*outputStream, name.c_str(), group, clust, homothety, position, angularPosition);
    }

    std::cerr << "@generatePacking:RandomClosePacking, Number of added bodies: " << obbs.size() << '\n';
    NumberOfParticles += (int)(obbs.size());
  };

  // This single line actually parses the file
  parser.parse(com);
}

int main(int argc, char* argv[]) {

  if (argc != 2) {
    std::cerr << "usage " << argv[0] << " <command-file>\n";
    return 0;
  } else {
    readCommands(argv[1]);
  }

  // Close file stream if it is still open
  if (fileStream.is_open()) {
    fileStream.close();
  }

  return 0;
}
