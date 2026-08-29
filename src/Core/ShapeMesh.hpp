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

#ifndef SHAPEMESH_HPP
#define SHAPEMESH_HPP

#include <cstddef>
#include <istream>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "Shape.hpp"
#include "ShapeSDF.hpp"
#include "vec3.hpp"

/// A triangle of the skin mesh, three indices into ShapeMesh::P
struct Tri {
  size_t a{0}, b{0}, c{0};
};

/**
   Triangle surface mesh of the r-shape skin (the level set f = 0 of ShapeSDF).

   Vertices carry an exact analytic normal (the SDF gradient at the vertex) and
   the skeleton primitive that generated them, so the mesh keeps the link back
   to the vertices, edges and faces of the Shape.

   The mesh lives in the body frame, exactly like the Shape it comes from. Since
   Shape::homothety() scales both the vertices and the Minkowski radius, two
   homothetic r-shapes are similar: one mesh plus a scale factor describes both.
*/
struct ShapeMesh {
  std::vector<vec3r> P;    ///< Vertex positions, on the surface f = 0
  std::vector<vec3r> N;    ///< Exact unit outward normal at each vertex
  std::vector<int> primType;    ///< SDF_PRIM_* of the closest skeleton primitive
  std::vector<size_t> primIndex;  ///< Its index in the Shape
  std::vector<Tri> tri;    ///< Triangles

  double epsilon{0.0};     ///< Geometric tolerance the mesh was built with
  double volume{0.0};      ///< Volume enclosed by the mesh (divergence theorem)
  double area{0.0};        ///< Surface area

  size_t nbVertices() const { return P.size(); }
  size_t nbTriangles() const { return tri.size(); }

  /// Volume enclosed by the closed triangle mesh, by the divergence theorem
  double computeVolume() const;
  /// Total surface area
  double computeArea() const;

  /// Write the mesh as a Wavefront OBJ (positions + exact vertex normals)
  void writeOBJ(std::ostream& os) const;
  /// Write the mesh as an ASCII PLY (positions + exact vertex normals)
  void writePLY(std::ostream& os) const;

  /// Write one named block of the ".rmsh" companion format (positions +
  /// exact normals + triangles), the cache read back by the viewer.
  void writeRmsh(std::ostream& os, const std::string& name) const;
};

/// Load a ".rmsh" companion file: the skin mesh of each shape, keyed by name.
/// Returns an empty map when the file is absent or unreadable.
std::map<std::string, ShapeMesh> loadRmsh(const std::string& path);

/// Parameters driving the mesher
struct ShapeMeshOptions {
  double epsilon{0.0};       ///< Target sag (max chordal error). 0 asks for the default below
  double epsilonRel{0.01};   ///< Default epsilon = epsilonRel * radius when epsilon == 0
  int newtonIters{2};        ///< Newton steps projecting each vertex onto f = 0
  int marginCells{2};        ///< Extra grid cells of padding around the shape AABB
  size_t maxCells{20000000};  ///< Budget on grid nodes; the step is coarsened to fit. 0 disables
};

/**
   Build the skin mesh of a r-shape by marching cubes on the ShapeSDF, then
   projecting every marching-cubes vertex onto f = 0 by a couple of Newton
   iterations along the analytic gradient. After projection the vertex position
   error falls to rounding, the normals are exact, and the only remaining error
   is the chordal error between vertices, set by epsilon.

   The grid step is chosen from epsilon through the sphere sag relation
   h ~ sqrt(8 epsilon R), and additionally capped at radius/2 so that a thin
   shell (thickness 2R) is never stepped over.
*/
ShapeMesh buildShapeMesh(const Shape& shape, const ShapeMeshOptions& opt = ShapeMeshOptions());

/// Parameters for the coplanar-merge simplification
struct SimplifyOptions {
  double angleToleranceDeg{0.25};  ///< Max normal deviation to consider two triangles coplanar
  /// Max distance of a vertex to the region plane, as a fraction of the mesh sag
  /// (ShapeMesh::epsilon). Dropping an interior vertex can move the surface by at
  /// most this much, so the merge error scales with the mesh's own tolerance and
  /// vanishes as the mesh is refined. A tolerance tied to the bounding box instead
  /// would stay fixed while the mesh error shrinks, and would end up dominating it.
  double planeToleranceSag{0.05};
  /// Fallback used only when the mesh carries no sag (epsilon == 0), relative to
  /// the bounding-box diagonal.
  double planeToleranceRel{1.0e-4};
};

/**
   Merge adjacent coplanar triangles and re-triangulate the resulting planar
   polygons. Flat regions (the offset faces of a r-shape are exact planes) are
   grown into planar patches, then re-triangulated from their boundary only, so
   the interior vertices vanish. Curved regions (rounded edges and corners) are
   left untouched. Only interior vertices are removed and boundary edges are
   preserved verbatim, so the mesh stays watertight.

   The surface is unchanged wherever the merged region is exactly planar. The
   tolerances above bound how much of the adjacent rounding may pass as flat, and
   hence how far the surface can move; both are tied to the mesh sag, so the
   residual shrinks with the mesh's own discretisation error.
*/
ShapeMesh simplifyCoplanar(const ShapeMesh& in, const SimplifyOptions& opt = SimplifyOptions());

#endif /* end of include guard: SHAPEMESH_HPP */
