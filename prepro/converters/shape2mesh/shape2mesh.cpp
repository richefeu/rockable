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
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

// shape2mesh: turn the r-shapes of a .shp file into surface meshes.
//
// Each r-shape (the Minkowski sum of its skeleton with a ball of radius R) is
// meshed on its exact signed distance function and written as an OBJ or PLY
// file, one file per shape, named <basename>_<shapeName>.<ext>.
//
// Usage:
//   shape2mesh <input.shp> [options]
// Options:
//   -f, --format obj|ply   output format          (default: obj)
//   -e, --epsilon <value>  absolute sag tolerance  (default: 1% of radius)
//   -o, --outdir <dir>     output directory        (default: alongside input)
//   -h, --help

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Core/Shape.hpp"
#include "Core/ShapeMesh.hpp"

static void usage(const char* prog) {
  std::printf(
      "Usage: %s <input.shp> [options]\n"
      "  -f, --format obj|ply   output format (default obj)\n"
      "  -e, --epsilon <value>  absolute sag tolerance (default 1%% of radius)\n"
      "  -s, --simplify [deg]    merge coplanar triangles (default angle 0.25 deg)\n"
      "  -r, --rmsh             write one <stem>.rmsh companion (all shapes) instead of per-shape files\n"
      "  -o, --outdir <dir>     output directory (default alongside input)\n"
      "  -h, --help\n",
      prog);
}

// The stem of a path: drop the directory and the last extension
static std::string baseStem(const std::string& path) {
  size_t slash = path.find_last_of("/\\");
  std::string file = (slash == std::string::npos) ? path : path.substr(slash + 1);
  size_t dot = file.find_last_of('.');
  return (dot == std::string::npos) ? file : file.substr(0, dot);
}

// Keep the output file names well-behaved
static std::string sanitize(const std::string& name) {
  std::string out;
  for (size_t i = 0; i < name.size(); ++i) {
    char c = name[i];
    out += (std::isalnum((unsigned char)c) || c == '-' || c == '_') ? c : '_';
  }
  return out.empty() ? std::string("shape") : out;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }

  std::string input;
  std::string format = "obj";
  std::string outdir;
  double epsilon = 0.0;  // 0 asks the mesher for its default (1% of the radius)
  bool simplify = false;
  double simplifyAngle = SimplifyOptions().angleToleranceDeg;  // the merge's own default
  bool rmshMode = false;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a == "-r" || a == "--rmsh") {
      rmshMode = true;
    } else if ((a == "-f" || a == "--format") && i + 1 < argc) {
      format = argv[++i];
    } else if ((a == "-e" || a == "--epsilon") && i + 1 < argc) {
      epsilon = std::atof(argv[++i]);
    } else if (a == "-s" || a == "--simplify") {
      simplify = true;
      // An optional numeric argument sets the coplanarity angle in degrees
      if (i + 1 < argc && argv[i + 1][0] != '-') {
        char* end = nullptr;
        double val = std::strtod(argv[i + 1], &end);
        if (end != argv[i + 1] && *end == '\0') simplifyAngle = val, ++i;
      }
    } else if ((a == "-o" || a == "--outdir") && i + 1 < argc) {
      outdir = argv[++i];
    } else if (a[0] != '-') {
      input = a;
    } else {
      std::fprintf(stderr, "Unknown or incomplete option: %s\n", a.c_str());
      usage(argv[0]);
      return 1;
    }
  }

  if (input.empty()) {
    std::fprintf(stderr, "No input .shp file given.\n");
    return 1;
  }
  if (format != "obj" && format != "ply") {
    std::fprintf(stderr, "Unknown format '%s' (use obj or ply).\n", format.c_str());
    return 1;
  }

  std::ifstream is(input.c_str());
  if (!is) {
    std::fprintf(stderr, "Cannot open '%s'.\n", input.c_str());
    return 1;
  }

  std::string dir = outdir;
  if (dir.empty()) {
    size_t slash = input.find_last_of("/\\");
    dir = (slash == std::string::npos) ? std::string(".") : input.substr(0, slash);
  }
  std::string stem = baseStem(input);

  // Collect the shapes first, so the output naming can depend on how many there
  // are: a single shape simply takes the input file name, while several shapes
  // are disambiguated by their own name (or an index)
  std::vector<Shape> shapes;
  std::string tok;
  while (is >> tok) {
    if (tok != "<") continue;  // Shape::read consumes the block up to '>'
    Shape shape;
    shape.read(is);
    if (!shape.vertex.empty()) shapes.push_back(shape);
  }

  if (shapes.empty()) {
    std::fprintf(stderr, "No shape found in '%s'.\n", input.c_str());
    return 1;
  }

  ShapeMeshOptions opt;
  opt.epsilon = epsilon;
  bool single = (shapes.size() == 1);

  // ".rmsh" companion mode: mesh every shape and write them all into one file
  // named after the input stem, so the viewer can find it next to the shapes.
  if (rmshMode) {
    std::string out = dir + "/" + stem + ".rmsh";
    std::ofstream os(out.c_str());
    if (!os) {
      std::fprintf(stderr, "Cannot write '%s'.\n", out.c_str());
      return 1;
    }
    os << "rmsh 1\n";
    for (size_t i = 0; i < shapes.size(); ++i) {
      // Auto sag proportional to each shape's size, so a large flat wall and a
      // small detailed grain both get a sensible mesh from one command.
      ShapeMeshOptions o = opt;
      if (epsilon <= 0.0) {
        AABB box;
        shapes[i].getAABB(box);
        o.epsilon = 0.004 * norm(box.max - box.min);
      }
      ShapeMesh mesh = buildShapeMesh(shapes[i], o);
      if (simplify) {
        SimplifyOptions sopt;
        sopt.angleToleranceDeg = simplifyAngle;
        mesh = simplifyCoplanar(mesh, sopt);
      }
      mesh.writeRmsh(os, shapes[i].name);
      std::printf("%-24s meshed (%zu vertices, %zu triangles)\n", shapes[i].name.c_str(), mesh.nbVertices(),
                  mesh.nbTriangles());
    }
    std::printf("Wrote %s (%zu shape(s)).\n", out.c_str(), shapes.size());
    return 0;
  }

  for (size_t i = 0; i < shapes.size(); ++i) {
    Shape& shape = shapes[i];
    ShapeMesh mesh = buildShapeMesh(shape, opt);
    size_t triFull = mesh.nbTriangles();
    if (simplify) {
      SimplifyOptions sopt;
      sopt.angleToleranceDeg = simplifyAngle;
      mesh = simplifyCoplanar(mesh, sopt);
    }

    // One shape: mirror the input file name. Several: append the shape name, or
    // an index when the shape carries no usable name
    std::string base = stem;
    if (!single) {
      std::string tag = shape.name.empty() ? ("shape" + std::to_string(i)) : sanitize(shape.name);
      base += "_" + tag;
    }
    std::string out = dir + "/" + base + "." + format;

    std::ofstream os(out.c_str());
    if (!os) {
      std::fprintf(stderr, "Cannot write '%s'.\n", out.c_str());
      return 1;
    }
    if (format == "obj") {
      mesh.writeOBJ(os);
    } else {
      mesh.writePLY(os);
    }

    if (simplify) {
      std::printf("%-24s -> %s  (%zu vertices, %zu triangles from %zu, vol %.6g)\n", shape.name.c_str(), out.c_str(),
                  mesh.nbVertices(), mesh.nbTriangles(), triFull, mesh.volume);
    } else {
      std::printf("%-24s -> %s  (%zu vertices, %zu triangles, vol %.6g)\n", shape.name.c_str(), out.c_str(),
                  mesh.nbVertices(), mesh.nbTriangles(), mesh.volume);
    }
  }

  std::printf("Meshed %zu shape(s).\n", shapes.size());
  return 0;
}
