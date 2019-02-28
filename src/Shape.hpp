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

#ifndef SHAPE_HPP_EEEA5A98
#define SHAPE_HPP_EEEA5A98

#include <vector>
#include <set>

#include "AABB.hpp"
#include "OBB.hpp"
#include "quat.hpp"

struct OBBnode {
  std::vector<size_t> faceID;  // Face contained in the OBB of the node
  OBB obb;                     // The sub-OBB of the node
  int level;                   // level in the tree (0 is root)
  int parent;                  // index position of the parent node in OBBtree.nodes
  std::vector<int> children;   // index position of the ??? in OBBtree.nodes
  OBBnode() : level(0), parent(-1) {}
};

struct OBBtree {
  std::vector<OBBnode> nodes;
  // Storage solution:
  // id:      0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 ...
  // level:   0 | 1 | 1 | 2 | 2 | 2 | 2 | 3 ...
  // childID: 1 | 3 | 5 | ...
  //          2 | 4 | 6 | ...
  // >>>>>>> FIXME: this solution is not the one used. RE-DO <<<< ======= ************ !!!!!!!!!!!

  void add(OBBnode& root) { nodes.push_back(root); }

  void add(int parent, OBBnode& node) {
    node.parent = parent;
    nodes[parent].children.push_back((int)nodes.size());
    nodes.push_back(node);
  }

  void clear() { nodes.clear(); }
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

  std::vector<size_t> vgroup;
  std::vector<size_t> egroup;
  std::vector<size_t> fgroup;
  std::vector<std::string> groupDefinitionSolutions;

  std::vector<std::vector<size_t> > edgesConnectedWithVertex;  ///< iedge = ___[ivertex][i]
  std::vector<std::vector<size_t> > facesConnectedWithVertex;  ///< iface = ___[ivertex][i]
  bool ConnectionDefined;                                      ///< Flag to know wether the connections have been
                                                               ///< defined

  OBBtree tree;       ///< The OBB tree of the model
  int OBBtreeLevel;   ///< Deepness of the OBB tree
  bool treeComputed;  ///< Flag to know wether the tree is built

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

  void defineVertexConnectivity();
  void buildOBBtree();

 private:
  void updateObb(OBBnode& node);
};

#endif /* end of include guard: SHAPE_HPP_EEEA5A98 */
