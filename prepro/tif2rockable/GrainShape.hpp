#ifndef GRAINSHAPE_HPP
#define GRAINSHAPE_HPP

#include <set>
#include <vector>
#include <utility>

#include "vec3.hpp"
#include "fastSort3.hpp"

struct GrainShape {
  int colorId;
  size_t volume;
  double R;
  vec3r pos;

  std::vector<vec3r> points;
  std::vector<vec3<size_t>> facets;

  GrainShape(size_t nsubs = 0) : colorId(-1), volume(0), R(0.0), pos() {
    subdivide(nsubs);
  }

  void saveRockableShape(std::ostream &os, size_t Id, double voxSize = 1.0) {

    std::set<std::pair<size_t, size_t>> edges;
		size_t v0,v1,v2;
    for (size_t i = 0; i < facets.size(); i++) {
			v0 = facets[i][0];
			v1 = facets[i][1];
			v2 = facets[i][2];
			fastSort3<size_t>(v0, v1, v2);
			std::pair<size_t, size_t> p01 = std::make_pair(v0, v1);
			std::pair<size_t, size_t> p02 = std::make_pair(v0, v2);
			std::pair<size_t, size_t> p12 = std::make_pair(v1, v2);
			edges.insert(p01);
			edges.insert(p02);
			edges.insert(p12);
    }

    os << "<\n";
		char name[256];
		snprintf(name, 256, "GRAIN%zu", Id);
    os << "name " << name << '\n';
    os << "radius " << R * voxSize << '\n';
    os << "preCompDone n\n";
		os << "MCnstep 200000\n"; // c'est une valeur ni trop grande nin trop petite
    os << "nv " << points.size() << '\n';
    for (size_t i = 0; i < points.size(); i++) {
      os << points[i] * voxSize << '\n';
    }
		os << "ne " << edges.size() << '\n';
		for (auto p : edges) {
			os << p.first << ' ' << p.second << '\n';
		}
    os << "nf " << facets.size() << '\n';
    for (size_t i = 0; i < facets.size(); i++) {
      os << "3 " << facets[i] << '\n';
    }
    os << "position " << pos * voxSize << '\n';
    os << "orientation 1 0 0 0" << '\n';
    os << ">\n";
  }

  void save(std::ostream &file) {
    file << pos << '\n';
    file << volume << '\n';
    file << colorId << '\n';
    file << R << '\n';
    file << points.size() << '\n';
    for (size_t i = 0; i < points.size(); i++) {
      file << points[i] << '\n';
    }
    file << facets.size() << '\n';
    for (size_t i = 0; i < facets.size(); i++) {
      file << facets[i] << '\n';
    }
  }

  void read(std::istream &file) {
    points.clear();
    facets.clear();
    file >> pos;
    file >> volume;
    file >> colorId;
    file >> R;
    size_t nbp;
    file >> nbp;
    vec3r ppos;
    for (size_t i = 0; i < nbp; i++) {
      file >> ppos;
      points.push_back(ppos);
    }
    size_t nbf;
    file >> nbf;
    vec3<size_t> tri;
    for (size_t f = 0; f < nbf; f++) {
      file >> tri;
      facets.push_back(tri);
    }
  }

  void setRadius(double R) {
    for (size_t p = 0; p < points.size(); p++) {
      points[p].normalize();
      points[p] *= R;
    }
  }

  void move(double dx, double dy, double dz) {
    vec3r depl(dx, dy, dz);
    for (size_t p = 0; p < points.size(); p++) {
      points[p] += depl;
    }
  }

  void subdivide(size_t nsubs = 0) {
    setIcosahedron();

    for (size_t s = 0; s < nsubs; s++) {
      size_t nf = facets.size();
      for (size_t f = 0; f < nf; f++) {
        size_t i0 = facets[f][0];
        size_t i1 = facets[f][1];
        size_t i2 = facets[f][2];

        vec3r p01 = 0.5 * (points[i0] + points[i1]);
        p01.normalize();
        size_t i01 = points.size();
        points.push_back(p01);

        vec3r p12 = 0.5 * (points[i1] + points[i2]);
        p12.normalize();
        size_t i12 = points.size();
        points.push_back(p12);

        vec3r p02 = 0.5 * (points[i0] + points[i2]);
        p02.normalize();
        size_t i02 = points.size();
        points.push_back(p02);

        facets[f].set(i01, i12, i02);
        facets.push_back(vec3<size_t>(i0, i01, i02));
        facets.push_back(vec3<size_t>(i1, i12, i01));
        facets.push_back(vec3<size_t>(i2, i02, i12));
      }
    }
  }

  void setIcosahedron() {
    points.clear();
    facets.clear();

#define X 0.525731112119133606
#define Z 0.850650808352039932

    points.push_back(vec3r(-X, 0.0, Z));
    points.push_back(vec3r(X, 0.0, Z));
    points.push_back(vec3r(-X, 0.0, -Z));
    points.push_back(vec3r(X, 0.0, -Z));
    points.push_back(vec3r(0.0, Z, X));
    points.push_back(vec3r(0.0, Z, -X));
    points.push_back(vec3r(0.0, -Z, X));
    points.push_back(vec3r(0.0, -Z, -X));
    points.push_back(vec3r(Z, X, 0.0));
    points.push_back(vec3r(-Z, X, 0.0));
    points.push_back(vec3r(Z, -X, 0.0));
    points.push_back(vec3r(-Z, -X, 0.0));

#undef X
#undef Z

    facets.push_back(vec3<size_t>(0, 4, 1));
    facets.push_back(vec3<size_t>(0, 9, 4));
    facets.push_back(vec3<size_t>(9, 5, 4));
    facets.push_back(vec3<size_t>(4, 5, 8));
    facets.push_back(vec3<size_t>(4, 8, 1));
    facets.push_back(vec3<size_t>(8, 10, 1));
    facets.push_back(vec3<size_t>(8, 3, 10));
    facets.push_back(vec3<size_t>(5, 3, 8));
    facets.push_back(vec3<size_t>(5, 2, 3));
    facets.push_back(vec3<size_t>(2, 7, 3));
    facets.push_back(vec3<size_t>(7, 10, 3));
    facets.push_back(vec3<size_t>(7, 6, 10));
    facets.push_back(vec3<size_t>(7, 11, 6));
    facets.push_back(vec3<size_t>(11, 0, 6));
    facets.push_back(vec3<size_t>(0, 1, 6));
    facets.push_back(vec3<size_t>(6, 1, 10));
    facets.push_back(vec3<size_t>(9, 0, 11));
    facets.push_back(vec3<size_t>(9, 11, 2));
    facets.push_back(vec3<size_t>(9, 2, 5));
    facets.push_back(vec3<size_t>(7, 2, 11));
  }
};

#endif /* end of include guard: GRAINSHAPE_HPP */
