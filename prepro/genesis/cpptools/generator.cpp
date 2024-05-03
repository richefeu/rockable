#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "kwParser.hpp"
#include "message.hpp"
#include "quat.hpp"
#include "vec3.hpp"

#include "addParticle.hpp"
#include "generatePacking_grid.hpp"
#include "generatePacking_wallBox.hpp"
#include "generateShape_cube.hpp"
#include "generateShape_cuboid.hpp"
#include "generateShape_rectangle.hpp"
#include "generateShape_rhombicuboctahedron.hpp"
#include "generateShape_sphere.hpp"
#include "generateShape_xyz_walls.hpp"

std::ostream* outputStream = &std::cout;
std::ofstream fileStream;

void readCommands(const char* name) {
  std::ifstream com(name);
  if (!com.is_open()) {
    std::cerr << msg::warn() << "@readCommands, Cannot read " << name << msg::normal() << std::endl;
    return;
  }

  kwParser parser;

  // Map "open" command to a function
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

  // Map "close" command to a function
  parser.kwMap["close"] = __DO__(com) {
    if (fileStream.is_open()) {
      fileStream.close();
      outputStream = &std::cout;
    }
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

  parser.kwMap["generateShape:rectangle_xz"] = __DO__(com) {
    std::string name;
    double radius;
    double sidex;
    double sidez;
    com >> name >> radius >> sidex >> sidez;
    generateShape_rectangle_xz(*outputStream, name.c_str(), radius, sidex, sidez);
    std::cerr << "A rectangle shape (in plane xz) has been generated\n";
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
    std::cerr << "@wallBox, Number of added bodies: " << nbAdded << '\n';
  };

  parser.kwMap["generatePacking:grid"] = __DO__(com) {
    std::string name;
    vec3r origBox;
    vec3r boxSize;
    vec3i n;
    int group;
    double homothety;
    int randQ;
    com >> name >> origBox >> boxSize >> n >> group >> homothety >> randQ;
    int nbAdded = generatePacking_grid(*outputStream, name.c_str(), origBox, boxSize, n, group, homothety, randQ);
    std::cerr << "@grid, Number of added bodies: " << nbAdded << '\n';
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
