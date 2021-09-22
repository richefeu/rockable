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

#ifndef SHAPE_HPP
#define SHAPE_HPP

#include <set>
#include <vector>

#include "AABB.hpp"
#include "OBB.hpp"
#include "OBBtree.hpp"
#include "quat.hpp"

struct subBox {
  size_t isub;
  int nbPoints;
};

/// Shape of a sphero-polyhedron
class Shape {
 public:
  std::vector<vec3r> vertex;                     ///< Vertex positions expressed in the body frame
                                                 ///< (Inertia principal axis) in the real world
  std::vector<std::pair<size_t, size_t> > edge;  ///< List of edges, each of which is defined by the
                                                 ///< two vertex local-number
  std::vector<std::vector<size_t> > face;        ///< List of faces, each of which is
                                                 ///< defined by a list of vertex local-number

  OBBtree<subBox> tree;  ///< The OBB tree of the model
  int OBBtreeLevel;      ///< Deepness of the OBB tree
  bool treeComputed;     ///< Flag to know wether the tree is built

  vec3r position;      ///< A saving of OG (after the computation of mass center)
  quat orientation;    ///< A saving of the rotation of eigen vectors resulting
                       ///< from the computation of inertia matrix
  std::string name;    ///< Name of the shape
  char preCompDone;    ///< Can be 'n' or 'y' to allow or not (resp.) the
                       ///< numerical computation of mass center, volume and I/m
  bool isSurface;      ///< If true, the shape is a surface (not a volume formed by a polyedron)
  double radius;       ///< Minskowski radius
  double volume;       ///< Volume
  vec3r inertia_mass;  ///< Diagonal terms of the matrix of inertia/mass (in the principal framework)
  OBB obb;             ///< Oriented Bounding Box
  size_t MCnstep;      ///< Number of steps for Monte Carlo integration
  double Rmax;         ///< Radius of the circumscribed sphere
  double Rswing;       ///< Maximum distance between a vertexe and the center

  // Ctor
  Shape();

  void rotate(const quat& Q);
  void homothety(const double H);

  void read(std::istream& is);
  void write(std::ostream& os);

  void fitObb();
  void getAABB(AABB& aabb);
  bool inside(const vec3r& point);
  void massProperties();
  void clean();

  void buildOBBtree();
};

#endif /* end of include guard: SHAPE_HPP */
