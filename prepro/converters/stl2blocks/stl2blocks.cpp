// stl2blocks: cut a closed STL mesh into blocks, along families of parallel planes,
// and write the blocks as Rockable shapes (sphero-polyhedra) and particles.
//
// Each block is cut along its own direction, with its own spacing of the planes, and
// the slabs obtained can be cut again (a tree of cuts). The Minkowski radius R is given
// first: the faces of the mesh are moved inwards by R, and each cutting plane removes a
// slice of thickness 2R, so that the blocks, rounded by R, fill the volume of the mesh
// and touch each other.
//
// Author: Vincent.Richefeu@univ-grenoble-alpes.fr
// Lab 3SR, Grenoble
// 2026

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "polyhTool.hpp"  // the copy in this folder (namespace stl2blocks)

#include "kwParser.hpp"
#include "message.hpp"

#include "Core/Shape.hpp"

using stl2blocks::polyhTool;

// A rule of cut: the blocks whose path matches 'pattern' are cut along 'direction'
struct CutRule {
  std::vector<std::string> pattern;  // root, then slab numbers or '*'
  vec3r direction;
  std::string mode;  // number, spacing, positions, fractions, random
  std::vector<double> values;
  int nbWildcards{0};
};

// A face of a cut: the plane x.u = c of the core, on one side of a cutting plane (side +1:
// the block is below the cutting plane, -1: above); the two sides of a cutting plane have
// the same id
struct CapPlane {
  int id{0};
  int side{0};
  vec3r u;
  double c{0.0};
};

// A final block (one connected part of a leaf of the tree of cuts)
struct Block {
  std::string path;
  int part{0};
  polyhTool::poly core;
  double coreVolume{0.0};
  int openSides{0};
  std::vector<CapPlane> caps;  // the planes of the cuts around the block
  std::vector<vec3r> glue;     // vertices added in its faces of the cuts
};

class Stl2Blocks {
 public:
  void readCommands(const char* name);
  void exec();

 private:
  bool readStl(const std::string& name, polyhTool::poly& mesh);
  const CutRule* findRule(const std::vector<std::string>& path) const;
  std::vector<double> planePositions(const CutRule& rule, double emin, double emax) const;
  void addLeaf(const polyhTool::poly& core, const std::vector<std::string>& path, const std::vector<CapPlane>& caps);
  void cut(const polyhTool::poly& core, const std::vector<std::string>& path, const std::vector<CapPlane>& caps);
  void addGluePoints(Block& A, const CapPlane& ca, Block& B, const CapPlane& cb);
  void writeFiles();

  // read in the command file
  std::string stlFileName{"mesh.stl"};
  std::string shapesFileName{"blocks.shp"};
  std::string inputFileName{"blocks.txt"};
  std::string shapeName{"block"};
  double radius{0.0};
  double homothety{1.0};
  vec3r translation;
  std::vector<CutRule> rules;
  int ParticlesGroup{0};
  int ParticlesCluster{0};
  size_t MCnstep{20000};
  double minVolume{0.0};

  std::vector<Block> blocks;
  std::vector<Shape> shapes;
  double meshVolume{0.0};
  int nbBadCaps{0};
  int nbOpenSides{0};
  int nbMergedParts{0};
  int nbUncut{0};
  int nextPlaneId{0};
  int gluePoints{1};
  size_t nbGluePoints{0};
};

static std::string pathString(const std::vector<std::string>& path) {
  std::string s;
  for (size_t i = 0; i < path.size(); i++) {
    if (i > 0) s += '.';
    s += path[i];
  }
  return s;
}

void Stl2Blocks::readCommands(const char* name) {
  std::ifstream com(name);
  if (!com.is_open()) {
    std::cerr << msg::warn() << "@Stl2Blocks::readCommands, cannot read " << name << msg::normal() << std::endl;
    exit(1);
  }

  kwParser parser;
  parser.kwMap["stlFile"] = __GET__(com, stlFileName);
  parser.kwMap["shapesFile"] = __GET__(com, shapesFileName);
  parser.kwMap["inputFile"] = __GET__(com, inputFileName);
  parser.kwMap["shapeName"] = __GET__(com, shapeName);
  parser.kwMap["MinkowskiRadius"] = __GET__(com, radius);
  parser.kwMap["homothety"] = __GET__(com, homothety);
  parser.kwMap["translation"] = __GET__(com, translation);
  parser.kwMap["ParticlesGroup"] = __GET__(com, ParticlesGroup);
  parser.kwMap["ParticlesCluster"] = __GET__(com, ParticlesCluster);
  parser.kwMap["MCnstep"] = __GET__(com, MCnstep);
  parser.kwMap["minVolume"] = __GET__(com, minVolume);
  parser.kwMap["gluePoints"] = __GET__(com, gluePoints);

  // cut <path> <dx dy dz> number <n>
  //                       spacing <s> <offset>
  //                       positions <k> <p1> ... <pk>
  //                       fractions <k> <f1> ... <fk>
  //                       random <k> <seed>
  parser.kwMap["cut"] = __DO__(is) {
    CutRule rule;
    std::string path;
    is >> path >> rule.direction >> rule.mode;
    std::stringstream ss(path);
    std::string seg;
    while (std::getline(ss, seg, '.')) {
      rule.pattern.push_back(seg);
      if (seg == "*") rule.nbWildcards++;
    }
    if (rule.pattern.empty() || rule.pattern[0] != "root") {
      std::cerr << msg::warn() << "@cut, the path '" << path << "' must start with 'root'" << msg::normal()
                << std::endl;
      exit(1);
    }
    if (rule.direction.length() == 0.0) {
      std::cerr << msg::warn() << "@cut, the direction of '" << path << "' is zero" << msg::normal() << std::endl;
      exit(1);
    }
    rule.direction.normalize();
    double v;
    if (rule.mode == "number") {
      is >> v;
      rule.values.push_back(v);
    } else if (rule.mode == "spacing" || rule.mode == "random") {
      double w;
      is >> v >> w;
      rule.values = {v, w};
    } else if (rule.mode == "positions" || rule.mode == "fractions") {
      size_t k;
      is >> k;
      for (size_t i = 0; i < k; i++) {
        is >> v;
        rule.values.push_back(v);
      }
    } else {
      std::cerr << msg::warn() << "@cut, unknown way of placing the planes: " << rule.mode
                << " (number, spacing, positions, fractions or random)" << msg::normal() << std::endl;
      exit(1);
    }
    rules.push_back(rule);
  };

  parser.parse(com);
}

// Read a STL file, binary or ascii. The nodes that have the same coordinates are merged.
bool Stl2Blocks::readStl(const std::string& name, polyhTool::poly& mesh) {
  std::ifstream file(name, std::ios::binary);
  if (!file) {
    std::cerr << msg::warn() << "Cannot read " << name << msg::normal() << std::endl;
    return false;
  }
  std::vector<vec3r> corners;  // 3 per triangle

  // A binary file has 84 + 50 n bytes; an ascii file begins with "solid" and has "facet"
  file.seekg(0, std::ios::end);
  std::streamoff fileSize = file.tellg();
  file.seekg(0, std::ios::beg);
  char head[80];
  uint32_t nbTriangles = 0;
  bool binary = false;
  if (fileSize >= 84) {
    file.read(head, 80);
    file.read((char*)&nbTriangles, 4);
    binary = (fileSize == 84 + 50 * (std::streamoff)nbTriangles);
  }

  if (binary) {
    for (uint32_t t = 0; t < nbTriangles; t++) {
      float buf[12];
      uint16_t attr;
      file.read((char*)buf, 12 * sizeof(float));
      file.read((char*)&attr, sizeof(uint16_t));
      for (int c = 0; c < 3; c++) corners.push_back(vec3r(buf[3 + 3 * c], buf[4 + 3 * c], buf[5 + 3 * c]));
    }
  } else {
    file.seekg(0, std::ios::beg);
    std::string token;
    while (file >> token) {
      if (token == "vertex") {
        vec3r p;
        file >> p;
        corners.push_back(p);
      }
    }
    corners.resize(3 * (corners.size() / 3));
  }

  std::map<vec3r, int, std::less<vec3r>> ids;
  mesh = polyhTool::poly();
  int nbDegenerate = 0;
  for (size_t t = 0; t < corners.size() / 3; t++) {
    int id[3];
    for (int c = 0; c < 3; c++) {
      auto it = ids.find(corners[3 * t + c]);
      if (it == ids.end()) {
        id[c] = (int)mesh.nodes.size();
        ids[corners[3 * t + c]] = id[c];
        polyhTool::pnode N;
        N.pos = corners[3 * t + c];
        mesh.nodes.push_back(N);
      } else {
        id[c] = it->second;
      }
    }
    if (id[0] == id[1] || id[1] == id[2] || id[2] == id[0]) {
      nbDegenerate++;
      continue;
    }
    polyhTool::pface F;
    F.nodes = {id[0], id[1], id[2]};
    mesh.faces.push_back(F);
  }
  std::cout << "> Read '" << name << "' (" << (binary ? "binary" : "ascii") << "): " << mesh.faces.size()
            << " triangles, " << mesh.nodes.size() << " nodes";
  if (nbDegenerate > 0) std::cout << " (" << nbDegenerate << " degenerate triangles removed)";
  std::cout << std::endl;
  return !mesh.faces.empty();
}

// The rule for a block: the one that matches its path with the fewest wildcards (the last
// one given if several match)
const CutRule* Stl2Blocks::findRule(const std::vector<std::string>& path) const {
  const CutRule* found = nullptr;
  for (const CutRule& r : rules) {
    if (r.pattern.size() != path.size()) continue;
    bool match = true;
    for (size_t i = 0; i < path.size() && match; i++) {
      if (r.pattern[i] != "*" && r.pattern[i] != path[i]) match = false;
    }
    if (match && (found == nullptr || r.nbWildcards <= found->nbWildcards)) found = &r;
  }
  return found;
}

// The positions of the planes along the direction, strictly inside the extent [emin, emax]
std::vector<double> Stl2Blocks::planePositions(const CutRule& rule, double emin, double emax) const {
  std::vector<double> p;
  double L = emax - emin;
  if (rule.mode == "number") {
    int n = (int)rule.values[0];
    for (int k = 1; k < n; k++) p.push_back(emin + k * L / n);
  } else if (rule.mode == "spacing") {
    double s = rule.values[0], offset = rule.values[1];
    if (s > 0.0) {
      for (double x = offset + std::ceil((emin - offset) / s) * s; x < emax; x += s) p.push_back(x);
    }
  } else if (rule.mode == "positions") {
    p = rule.values;
  } else if (rule.mode == "fractions") {
    for (double f : rule.values) p.push_back(emin + f * L);
  } else if (rule.mode == "random") {
    std::mt19937 gen((unsigned)rule.values[1]);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int k = 0; k < (int)rule.values[0]; k++) p.push_back(emin + U(gen) * L);
  }
  std::sort(p.begin(), p.end());
  std::vector<double> inside;
  for (double x : p) {
    if (x > emin && x < emax) inside.push_back(x);
  }
  return inside;
}

// A block that is not cut any more: its connected parts are the final blocks
void Stl2Blocks::addLeaf(const polyhTool::poly& core, const std::vector<std::string>& path,
                         const std::vector<CapPlane>& caps) {
  std::vector<polyhTool::poly> parts;
  polyhTool::split_components(core, parts);

  // Two parts closer than 2R (a neck of the mesh thinner than 2R, cut by the planes)
  // would overlap once rounded: they are kept together, in one block
  std::vector<int> group(parts.size());
  std::iota(group.begin(), group.end(), 0);
  std::function<int(int)> find = [&](int i) { return group[i] == i ? i : (group[i] = find(group[i])); };
  for (size_t i = 0; i < parts.size(); i++) {
    for (size_t j = i + 1; j < parts.size(); j++) {
      if (find(i) != find(j) && polyhTool::distance(parts[i], parts[j]) < 2.0 * radius) group[find(j)] = find(i);
    }
  }
  std::map<int, polyhTool::poly> merged;
  for (size_t i = 0; i < parts.size(); i++) {
    if (find(i) != (int)i) nbMergedParts++;
    polyhTool::append(merged[find(i)], parts[i]);
  }

  int ipart = 0;
  for (auto& m : merged) {
    Block B;
    B.path = pathString(path);
    B.part = ipart++;
    B.caps = caps;
    B.core = m.second;
    B.coreVolume = polyhTool::volume(m.second);
    B.openSides = polyhTool::open_sides(m.second);
    nbOpenSides += B.openSides;
    blocks.push_back(B);
  }
}

// Cut a block (its core: the block shrunk by the radius) with the rule that matches its
// path, then its slabs, recursively. A block without rule is a leaf: its connected parts
// are the final blocks.
void Stl2Blocks::cut(const polyhTool::poly& core, const std::vector<std::string>& path,
                     const std::vector<CapPlane>& caps) {
  const CutRule* rule = findRule(path);
  std::string indent(2 * path.size(), ' ');

  if (rule == nullptr) {
    addLeaf(core, path, caps);
    return;
  }

  const vec3r& u = rule->direction;
  double cmin = 1e300, cmax = -1e300;
  for (const auto& N : core.nodes) {
    cmin = std::min(cmin, N.pos * u);
    cmax = std::max(cmax, N.pos * u);
  }
  // The extent of the block itself (its core enlarged by the radius)
  std::vector<double> p = planePositions(*rule, cmin - radius, cmax + radius);

  std::cout << indent << pathString(path) << ": direction " << u << ", " << p.size() << " plane"
            << (p.size() > 1 ? "s" : "") << ", slabs";
  std::vector<std::pair<int, polyhTool::poly>> slabs;
  int firstId = nextPlaneId;
  nextPlaneId += (int)p.size();
  for (size_t k = 0; k <= p.size(); k++) {
    polyhTool::poly slab = core, sub;
    int bad = 0, status = 0;
    if (k > 0) {
      polyhTool::plan pl;
      pl.normal = u;
      pl.pos = (p[k - 1] + radius) * u;
      status = polyhTool::cut_mesh(slab, pl, sub, bad);
      nbBadCaps += bad;
      if (status == 2) continue;
      slab = sub;
    }
    if (k < p.size()) {
      polyhTool::plan pl;
      pl.normal = -u;
      pl.pos = (p[k] - radius) * u;
      status = polyhTool::cut_mesh(slab, pl, sub, bad);
      nbBadCaps += bad;
      if (status == 2) continue;
      slab = sub;
    }
    std::cout << ' ' << k;
    slabs.push_back(std::make_pair((int)k, slab));
  }
  std::cout << std::endl;
  if (slabs.empty()) {
    // All the slabs are thinner than 2R: the block is kept whole
    std::cout << msg::warn() << indent << "no slab is thicker than 2R: " << pathString(path) << " is not cut"
              << msg::normal() << std::endl;
    nbUncut++;
    slabs.push_back(std::make_pair(-1, core));
  }

  for (auto& s : slabs) {
    if (s.first < 0) {
      addLeaf(s.second, path, caps);
      continue;
    }
    std::vector<std::string> sub = path;
    sub.push_back(std::to_string(s.first));
    std::vector<CapPlane> subCaps = caps;
    size_t k = (size_t)s.first;
    if (k > 0) subCaps.push_back(CapPlane{firstId + (int)k - 1, -1, u, p[k - 1] + radius});
    if (k < p.size()) subCaps.push_back(CapPlane{firstId + (int)k, +1, u, p[k] - radius});
    cut(s.second, sub, subCaps);
  }
}

// The faces of a block in a plane of cut, in 2D
struct CapRegion {
  std::vector<std::vector<vec2r>> faces;        // counterclockwise
  std::vector<std::pair<vec2r, vec2r>> border;  // the sides of the region
  std::vector<vec2r> nodes;

  bool contains(const vec2r& q, double tol) const {
    for (const auto& F : faces) {
      bool in = true;
      for (size_t k = 0; k < F.size() && in; k++) {
        const vec2r& a = F[k];
        const vec2r& b = F[(k + 1) % F.size()];
        vec2r ab = b - a;
        double l = ab.length();
        if (l > 0.0 && polyhTool::cross2(a, b, q) < -tol * l) in = false;
      }
      if (in) return true;
    }
    return false;
  }
};

static CapRegion capRegion(const polyhTool::poly& P, const CapPlane& cap, const vec3r& e1, const vec3r& e2,
                           double tol) {
  CapRegion R;
  std::set<std::pair<int, int>> sides;
  std::vector<std::vector<int>> capFaces;
  for (const auto& F : P.faces) {
    bool inPlane = true;
    for (int n : F.nodes) {
      if (std::fabs(P.nodes[n].pos * cap.u - cap.c) > tol) inPlane = false;
    }
    if (!inPlane) continue;
    std::vector<int> v(F.nodes.begin(), F.nodes.end());
    capFaces.push_back(v);
    for (size_t k = 0; k < v.size(); k++) sides.insert(std::make_pair(v[k], v[(k + 1) % v.size()]));
  }
  auto to2 = [&](int n) { return vec2r(P.nodes[n].pos * e1, P.nodes[n].pos * e2); };
  std::set<int> nodes;
  for (const auto& v : capFaces) {
    std::vector<vec2r> f;
    for (int n : v) {
      f.push_back(to2(n));
      nodes.insert(n);
    }
    double area = 0.0;
    for (size_t k = 0; k < f.size(); k++) area += f[k].x * f[(k + 1) % f.size()].y - f[(k + 1) % f.size()].x * f[k].y;
    if (area < 0.0) std::reverse(f.begin(), f.end());
    R.faces.push_back(f);
  }
  for (const auto& s : sides) {
    if (!sides.count(std::make_pair(s.second, s.first)))
      R.border.push_back(std::make_pair(to2(s.first), to2(s.second)));
  }
  for (int n : nodes) R.nodes.push_back(to2(n));
  return R;
}

// Where the faces of two blocks face each other across a cutting plane, both blocks get
// the corners of the part they share: the nodes of each face that are in front of the
// other face, and the crossings of their sides. Each of these points is then a vertex of
// both blocks, 2R apart, that stickBCM can bond.
void Stl2Blocks::addGluePoints(Block& A, const CapPlane& ca, Block& B, const CapPlane& cb) {
  const vec3r& u = ca.u;
  vec3r e1 = (std::fabs(u.x) < 0.9) ? cross(u, vec3r(1, 0, 0)) : cross(u, vec3r(0, 1, 0));
  e1.normalize();
  vec3r e2 = cross(u, e1);
  double L = polyhTool::size(A.core) + polyhTool::size(B.core);
  double tol = 1e-9 * L;

  CapRegion RA = capRegion(A.core, ca, e1, e2, tol);
  CapRegion RB = capRegion(B.core, cb, e1, e2, tol);
  if (RA.faces.empty() || RB.faces.empty()) return;

  std::vector<vec2r> pts;
  for (const vec2r& q : RA.nodes) {
    if (RB.contains(q, tol)) pts.push_back(q);
  }
  for (const vec2r& q : RB.nodes) {
    if (RA.contains(q, tol)) pts.push_back(q);
  }
  for (const auto& sa : RA.border) {
    for (const auto& sb : RB.border) {
      vec2r r = sa.second - sa.first, s = sb.second - sb.first, w = sb.first - sa.first;
      double den = r.x * s.y - r.y * s.x;
      if (std::fabs(den) <= 1e-300) continue;  // parallel sides: their ends are in pts already
      double t = (w.x * s.y - w.y * s.x) / den;
      double v = (w.x * r.y - w.y * r.x) / den;
      if (t > 0.0 && t < 1.0 && v > 0.0 && v < 1.0) pts.push_back(sa.first + t * r);
    }
  }

  // Add the points that are not already vertices of the block
  double tolPoint = 1e-6 * radius;
  auto add = [&](Block& K, const CapRegion& RK, const CapPlane& cap) {
    std::vector<vec2r> have = RK.nodes;
    for (const vec2r& q : pts) {
      bool found = false;
      for (const vec2r& h : have) {
        if ((h - q).length() < tolPoint) {
          found = true;
          break;
        }
      }
      if (found) continue;
      have.push_back(q);
      K.glue.push_back(q.x * e1 + q.y * e2 + cap.c * u);
      nbGluePoints++;
    }
  };
  add(A, RA, ca);
  add(B, RB, cb);
}

void Stl2Blocks::writeFiles() {
  std::ofstream shp(shapesFileName);
  for (Shape& S : shapes) S.write(shp);
  std::cout << "> Wrote '" << shapesFileName << "'" << std::endl;

  std::ofstream conf(inputFileName);
  conf << "Rockable 29-11-2018\n";
  conf << "shapeFile " << shapesFileName << '\n';
  conf << "Particles " << shapes.size() << '\n';
  msg::bestPrecision(conf);
  for (Shape& S : shapes) {
    conf << S.name << ' ' << ParticlesGroup << ' ' << ParticlesCluster << " 1 " << S.position << " 0 0 0 0 0 0 "
         << S.orientation << " 0 0 0 0 0 0\n";
  }
  std::cout << "> Wrote '" << inputFileName << "'" << std::endl;
}

void Stl2Blocks::exec() {
  if (radius <= 0.0) {
    std::cerr << msg::warn() << "MinkowskiRadius must be positive" << msg::normal() << std::endl;
    exit(1);
  }
  polyhTool::poly mesh;
  if (!readStl(stlFileName, mesh)) exit(1);

  // Repair: dangling triangles, orientation of the faces, holes
  int nbDangling = polyhTool::remove_dangling_faces(mesh);
  if (nbDangling > 0) std::cout << "  " << nbDangling << " dangling triangles removed" << std::endl;
  int nbReversed = polyhTool::orient_faces(mesh);
  int nbNonManifold = 0;
  int nbHoles = polyhTool::close_holes(mesh, nbNonManifold);
  if (nbHoles > 0) nbReversed += polyhTool::orient_faces(mesh);
  if (nbReversed > 0) std::cout << "  " << nbReversed << " triangles reversed (orientation)" << std::endl;
  if (nbHoles > 0) std::cout << "  " << nbHoles << " holes closed" << std::endl;
  if (nbNonManifold > 0) {
    std::cout << msg::warn() << "  " << nbNonManifold << " sides shared by more than two triangles: the cuts may fail"
              << msg::normal() << std::endl;
  }

  for (auto& N : mesh.nodes) N.pos = homothety * N.pos + translation;
  meshVolume = polyhTool::volume(mesh);
  std::cout << "  volume " << meshVolume << " (homothety " << homothety << ")" << std::endl;

  // The core: the faces moved inwards by the radius
  double maxError = 0.0;
  int nbFlipped = 0, nbClamped = 0, nbReduced = 0;
  polyhTool::offset_inward(mesh, radius, maxError, nbFlipped, nbClamped, nbReduced);
  std::cout << "> Faces moved inwards by the Minkowski radius " << radius << std::endl;
  std::cout << "  largest gap to the moved planes: " << maxError << " R" << std::endl;
  if (nbFlipped > 0) {
    std::cout << msg::warn() << "  " << nbFlipped << " triangles turned over: parts of the mesh thinner than 2R"
              << msg::normal() << std::endl;
  }
  if (nbClamped > 0) std::cout << "  " << nbClamped << " nodes moved by 4R only (sharp corners)" << std::endl;
  if (nbReduced > 0) {
    std::cout << "  " << nbReduced << " nodes moved less, so that no triangle is turned over (triangles smaller than R)"
              << std::endl;
  }

  std::cout << "> Cuts (path: direction, planes, slabs kept)" << std::endl;
  cut(mesh, std::vector<std::string>(1, "root"), std::vector<CapPlane>());
  if (nbBadCaps > 0) {
    std::cout << msg::warn() << "  " << nbBadCaps
              << " faces of the cuts have a crossed outline, or could not be closed (the mesh moved\n             "
                 "inwards crosses itself: details of the mesh smaller than R)"
              << msg::normal() << std::endl;
  }
  if (nbMergedParts > 0) {
    std::cout << "  " << nbMergedParts << " parts closer than 2R to another part of the same block are kept with it"
              << std::endl;
  }
  if (nbUncut > 0) std::cout << "  " << nbUncut << " blocks not cut (all their slabs thinner than 2R)" << std::endl;
  if (nbOpenSides > 0) {
    std::cout << msg::warn() << "  " << nbOpenSides << " sides of the blocks are open (blocks not closed)"
              << msg::normal() << std::endl;
  }

  // The final blocks, with simple faces
  size_t nbRemoved = 0;
  std::vector<Block> kept;
  for (Block& B : blocks) {
    if (B.core.nodes.size() < 4 || B.coreVolume <= minVolume) {
      nbRemoved++;
      continue;
    }
    polyhTool::simplify_faces(B.core);
    kept.push_back(B);
  }
  blocks = kept;
  if (nbRemoved > 0) std::cout << "  " << nbRemoved << " tiny parts removed (core volume <= minVolume)" << std::endl;

  // The vertices for the glue, on both sides of each cutting plane
  if (gluePoints) {
    std::map<int, std::vector<std::pair<size_t, size_t>>> onPlane;  // id -> (block, cap)
    for (size_t b = 0; b < blocks.size(); b++) {
      for (size_t c = 0; c < blocks[b].caps.size(); c++) onPlane[blocks[b].caps[c].id].push_back(std::make_pair(b, c));
    }
    for (auto& pl : onPlane) {
      for (auto& a : pl.second) {
        const CapPlane& ca = blocks[a.first].caps[a.second];
        if (ca.side != 1) continue;
        for (auto& b : pl.second) {
          const CapPlane& cb = blocks[b.first].caps[b.second];
          if (cb.side != -1) continue;
          addGluePoints(blocks[a.first], ca, blocks[b.first], cb);
        }
      }
    }
    std::cout << "> " << nbGluePoints << " vertices added in the faces of the cuts, for the glue (stickBCM)"
              << std::endl;
    double capArea = 0.0;  // the faces of the cuts, counted once (from below the cutting planes)
    for (const Block& B : blocks) {
      for (const CapPlane& cap : B.caps) {
        if (cap.side != 1) continue;
        vec3r e1 = (std::fabs(cap.u.x) < 0.9) ? cross(cap.u, vec3r(1, 0, 0)) : cross(cap.u, vec3r(0, 1, 0));
        e1.normalize();
        CapRegion R = capRegion(B.core, cap, e1, cross(cap.u, e1), 1e-9 * polyhTool::size(B.core));
        for (const auto& f : R.faces) {
          for (size_t k = 0; k < f.size(); k++)
            capArea += 0.5 * (f[k].x * f[(k + 1) % f.size()].y - f[(k + 1) % f.size()].x * f[k].y);
        }
      }
    }
    std::cout << "  area of the faces of the cuts: " << capArea << std::endl;
  }

  std::cout << "> Mass properties (Monte Carlo, " << MCnstep << " points per block)" << std::endl;
  double totalVolume = 0.0;
  for (Block& B : blocks) {
    Shape S;
    S.name = shapeName + "_" + std::to_string(shapes.size());
    S.radius = radius;
    S.MCnstep = MCnstep;
    for (const auto& N : B.core.nodes) S.vertex.push_back(N.pos);
    for (const vec3r& g : B.glue) S.vertex.push_back(g);
    for (const auto& E : B.core.edges) S.edge.push_back(std::make_pair((size_t)E.node0, (size_t)E.node1));
    for (const auto& F : B.core.faces) S.face.push_back(std::vector<size_t>(F.nodes.begin(), F.nodes.end()));

    std::streambuf* coutBuf = std::cout.rdbuf();  // massProperties is talkative
    std::ostringstream quiet;
    std::cout.rdbuf(quiet.rdbuf());
    S.massProperties();
    std::cout.rdbuf(coutBuf);
    S.preCompDone = 'y';
    totalVolume += S.volume;

    std::cout << "  " << std::left << std::setw(12) << S.name << " " << std::setw(16) << B.path << " part "
              << std::setw(3) << B.part << " nv " << std::setw(4) << S.vertex.size() << " nf " << std::setw(4)
              << S.face.size() << " volume " << S.volume << std::right;
    if (B.openSides > 0) std::cout << " (" << B.openSides << " open sides)";
    std::cout << std::endl;
    shapes.push_back(S);
  }
  double vmin = 1e300;
  std::string smallest;
  for (const Shape& S : shapes) {
    if (S.volume < vmin) {
      vmin = S.volume;
      smallest = S.name;
    }
  }
  if (!shapes.empty())
    std::cout << "  smallest block: " << smallest << ", volume " << vmin << " (it sets the time step)" << std::endl;
  std::cout << "> " << shapes.size() << " blocks, total volume " << totalVolume << " ("
            << 100.0 * totalVolume / meshVolume << "% of the mesh)" << std::endl;

  writeFiles();
}

int main(int argc, char const* argv[]) {
  if (argc != 2) {
    std::cout << "usage: stl2blocks <command file>" << std::endl;
    return 0;
  }
  Stl2Blocks S2B;
  S2B.readCommands(argv[1]);
  S2B.exec();
  return 0;
}
