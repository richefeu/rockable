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

#ifndef SHAPESDF_HPP
#define SHAPESDF_HPP

#include <map>
#include <vector>

#include "Shape.hpp"
#include "vec3.hpp"

/// Type of skeleton primitive that carries a closest point
#define SDF_PRIM_VERTEX 0
#define SDF_PRIM_EDGE 1
#define SDF_PRIM_FACE 2

/// Which feature of a polygon holds the closest point
#define SDF_FEATURE_INTERIOR 0
#define SDF_FEATURE_EDGE 1
#define SDF_FEATURE_VERTEX 2

/// Result of a closest-point query on the skeleton
struct ClosestPoint {
  double distance{0.0};   ///< Unsigned distance from the query point to the skeleton
  vec3r point;            ///< Closest point, lying on the skeleton
  int primType{0};        ///< SDF_PRIM_VERTEX, SDF_PRIM_EDGE or SDF_PRIM_FACE
  size_t primIndex{0};    ///< Index of that primitive in the Shape
};

/**
   Exact signed distance function of a r-shape (sphero-polyhedron).

   A r-shape is the Minkowski sum S = K + B(R) of its skeleton K (the vertices,
   edges and faces held by the Shape) with a ball of radius R = Shape::radius.
   Its signed distance function is therefore

       f(x) = s(x) - R

   where s(x) is the signed distance to the skeleton, negative strictly inside
   the polyhedron when the faces enclose a volume, and unsigned otherwise
   (open shell, free edges, isolated vertices). This holds exactly everywhere,
   not only near the surface.

   The gradient is analytic:

       grad f(x) = sign(x) * (x - p*(x)) / ||x - p*(x)||

   where p*(x) is the closest point of the skeleton. Surface normals obtained
   this way are exact, which is what makes a coarse skin mesh render well.

   The queries use BVH.
*/
class ShapeSDF {
 public:
  explicit ShapeSDF(const Shape& s);

  /// Closest point of the skeleton, and the unsigned distance to it
  ClosestPoint closestPoint(const vec3r& x) const;

  /// Unsigned distance to the skeleton
  double unsignedSkeletonDistance(const vec3r& x) const;

  /// Signed distance to the skeleton (negative strictly inside a closed polyhedron)
  double signedSkeletonDistance(const vec3r& x) const;

  /// The signed distance function of the r-shape itself: f(x) = s(x) - R
  double value(const vec3r& x) const;

  /// Analytic gradient of value(). Unit length, except where it is undefined
  /// (on the skeleton itself, which only happens when radius == 0)
  vec3r gradient(const vec3r& x) const;

  /// Same predicate as Shape::inside(), obtained from the sign of value()
  bool inside(const vec3r& x) const { return value(x) < 0.0; }

  /// True when the faces were found to enclose a volume, so that the sign is meaningful
  bool isClosedVolume() const { return closedVolume; }

  /// True when the shape is a filled solid with a thick interior: a closed
  /// volume, or a face soup read as a solid by ray casting. A shape that is NOT
  /// a solid (a sphere, a capsule, a plate, an open surface) is entirely at the
  /// R scale, so the grid step must be capped to resolve it
  bool isSolid() const { return closedVolume || rayCastSign; }

  /// Signed volume of the skeleton polyhedron, after coherent orientation.
  /// Zero when the shape is not a closed volume. This is the volume of K, not of S.
  double skeletonVolume() const { return skelVolume; }

 private:
  const Shape* shp;

  bool closedVolume{false};       ///< The faces form a closed orientable surface, and !isSurface
  bool rayCastSign{false};        ///< Faces exist and !isSurface but are not a clean manifold:
                                  ///< the sign is taken by ray casting, as Shape::inside() does
  double skelVolume{0.0};         ///< Signed volume enclosed by the oriented faces
  std::vector<vec3r> faceNormal;  ///< Outward unit normal of each face (coherently oriented)
  std::vector<std::vector<size_t> > orientedFace;  ///< Faces, re-oriented outward

  /// Faces sharing a given (sorted) vertex pair
  std::map<std::pair<size_t, size_t>, std::vector<size_t> > edgeFaces;
  /// Faces incident to a given vertex
  std::vector<std::vector<size_t> > vertexFaces;

  /// Closest point on one polygonal face, together with the feature that carries it
  struct FaceClosest {
    double d2{0.0};
    vec3r point;
    int feature{SDF_FEATURE_INTERIOR};
    size_t ia{0};  ///< For SDF_FEATURE_EDGE and SDF_FEATURE_VERTEX: global vertex index
    size_t ib{0};  ///< For SDF_FEATURE_EDGE: the other global vertex index
  };

  /// One skeleton primitive (a vertex, an edge or a face) referenced by the BVH
  struct Prim {
    int type{0};      ///< SDF_PRIM_VERTEX, SDF_PRIM_EDGE or SDF_PRIM_FACE
    size_t index{0};  ///< Its index in the Shape (face index refers to orientedFace)
  };

  /// A node of the bounding-volume hierarchy over the primitives. A leaf spans
  /// [start, start+count) of the reordered id array; an internal node has two
  /// children and count == 0
  struct BVHNode {
    AABB box;
    int left{-1};
    int right{-1};
    int start{0};
    int count{0};
  };

  std::vector<Prim> prims;         ///< All primitives, indexed by the all-prim BVH
  std::vector<AABB> primBox;       ///< AABB of each primitive
  std::vector<vec3r> primCentroid; ///< Centroid of each primitive, for the split
  std::vector<BVHNode> bvhNodes;   ///< BVH over all primitives (for closest point)
  std::vector<int> bvhIds;         ///< Primitive ids, reordered by the BVH build

  std::vector<BVHNode> faceNodes;  ///< BVH over faces only (for the sign)
  std::vector<int> faceIds;        ///< Face ids, reordered by the face BVH build

  void buildAdjacency();
  void orientFacesOutward();
  void buildBVH();

  /// Recursive median split; returns the index of the created node
  int buildNode(std::vector<int>& ids, std::vector<BVHNode>& nodes, int start, int count) const;

  /// Squared distance from a point to an AABB (0 inside)
  static double sqDistToBox(const vec3r& x, const AABB& b);

  FaceClosest closestOnFace(const vec3r& x, size_t f) const;

  /// Nearest face to x (its FaceClosest and index), found through the face BVH.
  /// Returns false when there is no usable face
  bool nearestFace(const vec3r& x, FaceClosest& out, size_t& faceOut) const;

  /// Angle-weighted pseudonormal at the closest point of the face set.
  /// Bearentzen & Aanaes, "Signed distance computation using the angle weighted
  /// pseudonormal", IEEE TVCG 2005.
  vec3r pseudoNormal(const FaceClosest& fc, size_t f) const;

  /// -1 strictly inside the solid, +1 outside. Uses the exact angle-weighted
  /// pseudonormal for a clean manifold, ray casting for an imperfect face soup
  /// (same test as Shape::inside), and +1 when the sign is undefined (surface,
  /// or no faces)
  double signAt(const vec3r& x) const;

  /// Point-in-polyhedron by ray casting along +x, counting face crossings.
  /// Matches the polyhedron test of Shape::inside(), and tolerates a face soup
  /// that is not a clean 2-manifold
  bool rayCastInside(const vec3r& x) const;

  double angleAtVertex(size_t f, size_t v) const;
};

#endif /* end of include guard: SHAPESDF_HPP */
