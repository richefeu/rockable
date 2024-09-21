#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <locale>
#include <string>

#include "kwParser.hpp"
#include "message.hpp"
#include "nanoExprParser.hpp"
#include "quat.hpp"
#include "transformation.hpp"
#include "vec3.hpp"

Transformation<double> globalTransformation;
quat individualParticleRotation;

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

void inplace_trim(std::string& line) {
  // trim from start (in place)
  line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) { return !std::isspace(ch); }));
  // trim from end (in place)
  line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(),
             line.end());
}

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

  parser.kwMap["tranformation:reset"] = __DO__() { globalTransformation.reset(); };

  parser.kwMap["addParticle"] = __DO__(com) {
    std::string name;
    int group;
    int cluster;
    double homothety;
    vec3r position;
    quat angularPosition;
    com >> name >> group >> cluster >> homothety >> position >> angularPosition;
    angularPosition = angularPosition * individualParticleRotation;
    globalTransformation.apply(position);
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
    int clust;
    double homothety;
    int randQ;
    com >> name >> origBox >> boxSize >> n >> group >> clust >> homothety >> randQ;
    int nbAdded =
        generatePacking_grid(*outputStream, name.c_str(), origBox, boxSize, n, group, clust, homothety, randQ);
    std::cerr << "@grid, Number of added bodies: " << nbAdded << '\n';
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
