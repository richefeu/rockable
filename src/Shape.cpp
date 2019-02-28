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

// common hearders
#include "AABB.hpp"
#include "kwParser.hpp"
#include "message.hpp"
#include "polyhTool.hpp"

// local hearder
#include "Shape.hpp"

Shape::Shape()
    : ConnectionDefined(false),
      OBBtreeLevel(0),
      treeComputed(false),
      preCompDone('n'),
      isSurface(false),
      radius(0.0),
      volume(0.0),
      MCnstep(10000),
      Rmax(0.0),
      Rswing(0.0) {}

void Shape::rotate(const quat& Q) {
  for (size_t i = 0; i < vertex.size(); i++) {
    vertex[i] = Q * vertex[i];
  }
}

void Shape::homothety(const double H) {
  for (size_t i = 0; i < vertex.size(); i++) {
    vertex[i] = H * vertex[i];
  }
  radius *= H;
}

void Shape::read(std::istream& is) {
  kwParser parser;
  parser.breakStr = ">";
  parser.kwMap["name"] = __GET__(is, name);
  parser.kwMap["radius"] = __GET__(is, radius);
  parser.kwMap["preCompDone"] = __GET__(is, preCompDone);
  parser.kwMap["volume"] = __GET__(is, volume);
  parser.kwMap["I/m"] = __GET__(is, inertia_mass);
  parser.kwMap["obb.extent"] = __GET__(is, obb.extent);
  parser.kwMap["obb.center"] = __GET__(is, obb.center);
  parser.kwMap["obb.e1"] = __GET__(is, obb.e[0]);
  parser.kwMap["obb.e2"] = __GET__(is, obb.e[1]);
  parser.kwMap["obb.e3"] = __GET__(is, obb.e[2]);
  parser.kwMap["OBBtreeLevel"] = __GET__(is, OBBtreeLevel);
  parser.kwMap["position"] = __GET__(is, position);
  parser.kwMap["orientation"] = __GET__(is, orientation);
  parser.kwMap["isSurface"] = __DO__() { isSurface = true; };
  parser.kwMap["MCnstep"] = __GET__(is, MCnstep);
  parser.kwMap["nv"] = __DO__(is) {
    size_t nv;
    is >> nv;
    vec3r pos;
    for (size_t v = 0; v < nv; ++v) {
      is >> pos;
      vertex.push_back(pos);
    }
  };
  parser.kwMap["ne"] = __DO__(is) {
    size_t ne;
    is >> ne;
    size_t from, to;
    for (size_t e = 0; e < ne; ++e) {
      is >> from >> to;
      edge.push_back(std::pair<size_t, size_t>(from, to));
    }
  };
  parser.kwMap["nf"] = __DO__(is) {
    size_t nf;
    is >> nf;
    std::vector<size_t> faceNodes;
    for (uint f = 0; f < nf; ++f) {
      size_t nb_node_face;
      is >> nb_node_face;
      std::vector<size_t> ids;
      for (size_t n = 0; n < nb_node_face; ++n) {
        size_t id;
        is >> id;
        ids.push_back(id);
      }
      face.push_back(ids);
    }
  };
  /*
  // The groups must to defined after the elements (vertices, edges and faces)
  have been defined parser.kwMap["defaultGroup"] = __DO__(is) {
          // TODO
  };

  // groupInPolygonXZ
  // groupInBox
  // egroup, vgroup, fgroup


  // egroup
  // 0 2
  // 1 2

  // egroupRange 1 10 2

  */

  parser.parse(is);

  // Limit MCnstep value within a reasonnable range
  if (MCnstep > 10000000) MCnstep = 10000000;
  if (MCnstep < 1000) MCnstep = 1000;
}

void Shape::write(std::ostream& os) {
  msg::bestPrecision(os);
  os << "<" << std::endl;
  os << "name " << name << std::endl;
  os << "radius " << radius << std::endl;
  os << "preCompDone " << preCompDone << std::endl;
  os << "volume " << volume << std::endl;
  os << "I/m " << inertia_mass << std::endl;
  os << "obb.extent " << obb.extent << std::endl;
  os << "obb.center " << obb.center << std::endl;
  os << "obb.e1 " << obb.e[0] << std::endl;
  os << "obb.e2 " << obb.e[1] << std::endl;
  os << "obb.e3 " << obb.e[2] << std::endl;
  os << "OBBtreeLevel " << OBBtreeLevel << std::endl;
  os << "position " << position << std::endl;
  os << "orientation " << orientation << std::endl;
  if (isSurface == true) os << "isSurface" << std::endl;
  os << "MCnstep " << MCnstep << std::endl;

  os << "nv " << vertex.size() << std::endl;
  for (size_t v = 0; v < vertex.size(); v++) {
    os << vertex[v] << std::endl;
  }
  msg::normalPrecision(os);

  os << "ne " << edge.size() << std::endl;
  for (size_t e = 0; e < edge.size(); e++) {
    os << edge[e].first << " " << edge[e].second << std::endl;
  }

  os << "nf " << face.size() << std::endl;
  for (size_t f = 0; f < face.size(); f++) {
    os << face[f].size();
    for (size_t n = 0; n < face[f].size(); n++) {
      os << " " << face[f][n];
    }
    os << std::endl;
  }

  os << ">" << std::endl;
}

// See pdf document "Minimum-Area Rectangle Containing a Convex Polygon" (@see
// http://www.geometrictools.com/)
void Shape::fitObb() {
  if (face.empty()) {
    std::cerr << "No face in this shape. Cannot define the Oriented Bounding Box" << std::endl;
    return;
  }

  double ext_test;
  vec3r e0, e1, e2;
  vec3r P0, Vtest;
  double minArea = 1.0e20;
  for (size_t f = 0; f < face.size(); ++f) {
    size_t nvert_in_face = face[f].size();

    P0 = vertex[face[f][0]];
    e0 = (vertex[face[f][1]] - P0).normalized();
    e1 = vertex[face[f][nvert_in_face - 1]] - P0;
    e2 = cross(e0, e1).normalized();
    e1 = cross(e2, e0);

    double ext0_min, ext1_min, ext2_min;
    double ext0_max, ext1_max, ext2_max;
    ext0_min = ext1_min = ext2_min = 1.0e20;
    ext0_max = ext1_max = ext2_max = -1.0e20;
    for (uint v = 0; v < vertex.size(); ++v) {
      Vtest = vertex[v] - P0;
      ext_test = Vtest * e0;
      if (ext_test < ext0_min)
        ext0_min = ext_test;
      else if (ext_test > ext0_max)
        ext0_max = ext_test;
      ext_test = Vtest * e1;
      if (ext_test < ext1_min)
        ext1_min = ext_test;
      else if (ext_test > ext1_max)
        ext1_max = ext_test;
      ext_test = Vtest * e2;
      if (ext_test < ext2_min)
        ext2_min = ext_test;
      else if (ext_test > ext2_max)
        ext2_max = ext_test;
    }
    double area = (ext0_max - ext0_min) * (ext1_max - ext1_min) * (ext2_max - ext2_min);
    if (area < minArea) {
      minArea = area;
      obb.e[0] = e0;
      obb.e[1] = e1;
      obb.e[2] = e2;
      obb.center = P0 + 0.5 * ((ext0_max + ext0_min) * e0 + (ext1_max + ext1_min) * e1 + (ext2_max + ext2_min) * e2);
      obb.extent.set(0.5 * (ext0_max - ext0_min), 0.5 * (ext1_max - ext1_min), 0.5 * (ext2_max - ext2_min));
    }
  }

  obb.enlarge(radius);  // Add the Minskowski radius
}

// Say whether a point is inside the shape
bool Shape::inside(const vec3r& point) {
  // === inside POLYHEDRON ===
  if (!isSurface) {
    vec3r v1, v2, v3;
    size_t nb_intersect = 0;
    for (size_t f = 0; f < face.size(); ++f) {
      v1 = vertex[face[f][0]];
      for (size_t v = 1; v < face[f].size() - 1; ++v) {
        v2 = vertex[face[f][v]];
        v3 = vertex[face[f][v + 1]];
        if (geoTool::intersectTriangle(point, vec3r::unit_x(), v1, v2, v3) > 0) {
          ++(nb_intersect);
          break;
        }
      }
    }

    if (nb_intersect % 2 != 0) return true;  // Should be optimized by the compilator
  }

  // === inside NODES or EDGES ===
  vec3r E, V;
  double d2;
  for (size_t e = 0; e < edge.size(); ++e) {
    size_t i0 = edge[e].first;
    size_t i1 = edge[e].second;
    E = vertex[i1] - vertex[i0];
    V = point - vertex[i0];
    double r = (V * E) / (E * E);
    if (r <= 0.0) {  // First extremity (a sphere)
      d2 = norm2(point - vertex[i0]);
      if (d2 < radius * radius) return true;
    } else if (r >= 1.0) {  // Second extremity (a sphere)
      d2 = norm2(point - vertex[i1]);
      if (d2 < radius * radius) return true;
    } else {  // The tube part
      d2 = norm2(point - (vertex[i0] + r * E));
      if (d2 < radius * radius) return true;
    }
  }

  // === inside FACE thickness ===
  for (size_t f = 0; f < face.size(); ++f) {

    // First, we project the node position onto the face plane.
    vec3r P;
    size_t nb_vertices = face[f].size();
    vec3r v = point - vertex[face[f][0]];
    vec3r v1 = vertex[face[f][1]] - vertex[face[f][0]];
    v1.normalize();
    vec3r v2 = vertex[face[f][nb_vertices - 1]] - vertex[face[f][0]];
    vec3r n = cross(v1, v2);
    n.normalize();
    double dn = v * n;
    if (dn < 0.0) {
      n = -n;
      dn = -dn;
    }
    P = point - dn * n;

    // Then, we check whether the projected point is inside the face (a
    // polygon). We use the crossing number algorithm (also known as even-odd
    // rule algorithm)
    size_t ODD = 0;
    v2 = cross(n, v1);  // both n and v1 are normalized
    double ori1 = P * v1;
    double ori2 = P * v2;
    for (size_t iva = 0; iva < nb_vertices; ++iva) {
      size_t ivb = iva + 1;
      if (ivb == nb_vertices) ivb = 0;
      double pa1 = vertex[face[f][iva]] * v1;
      double pb1 = vertex[face[f][ivb]] * v1;
      double pa2 = vertex[face[f][iva]] * v2;
      double pb2 = vertex[face[f][ivb]] * v2;

      // @see http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
      // @see http://alienryderflex.com/polygon/
      if ((pa2 < ori2 && pb2 >= ori2) || (pb2 < ori2 && pa2 >= ori2)) {
        if (pa1 + (ori2 - pa2) / (pb2 - pa2) * (pb1 - pa1) < ori1) {
          ODD = 1 - ODD;
        }
      }
    }

    if (ODD) {
      dn -= radius;
      if (dn < 0.0) return true;
    }
  }

  // If we arrive here, the point is not inside
  return false;
}

// This method aims to remove duplicate edges
void Shape::clean() {
  std::set<std::pair<size_t, size_t> > edgeSet;
  for (size_t e = 0; e < edge.size(); e++) {
    std::pair<size_t, size_t> P;
    if (edge[e].first < edge[e].second) {
      P.first = edge[e].first;
      P.second = edge[e].second;
    }
    else {
      P.first = edge[e].second;
      P.second = edge[e].first;
    }
    edgeSet.insert(P); // That way, duplication is not allowed 
  }
  
  size_t nbePrev = edge.size();
  size_t nbeNew = edgeSet.size();
  
  edge.clear();
  std::copy(edgeSet.begin(), edgeSet.end(), std::back_inserter(edge));
  
  std::cout << name << ":\n";
  std::cout << nbePrev - nbeNew << " duplicated edges have been removed\n";
}

// Precompute by means of Monte Carlo integration
// the volume, mass-center, inertia and body-frame (so that the inertia matrix
// is diagonal)
void Shape::massProperties() {
  std::cout << std::endl;
  std::cout << "Computation of mass properties (volume, mass-center, inertia "
               "and body-frame) for shape "
            << name << std::endl
            << std::flush;

  // 1- Get the bounding volume
  AABB box(vertex);
  box.enlarge(radius);
  double Vbox = (box.max.x - box.min.x) * (box.max.y - box.min.y) * (box.max.z - box.min.z);

  // 2- Monte Carlo integration to compute the volume and the mass-center
  double int_vol = 0.0, vol_err = 0.0;
  vec3r OG;
  std::vector<double> vv(3);
  polyhTool::sobolSequence(-3, vv);  // Initialize the Sobol sequence
  vec3r pt3;
  for (size_t i = 0; i < MCnstep; ++i) {
    polyhTool::sobolSequence(3, vv);
    pt3.set(box.min.x + vv[0] * (box.max.x - box.min.x), box.min.y + vv[1] * (box.max.y - box.min.y),
            box.min.z + vv[2] * (box.max.z - box.min.z));
    if (inside(pt3)) {
      OG += pt3;
      int_vol += 1.0;
      vol_err += 1.0;
    }
  }

  if (int_vol == 0.0) {
    std::cout << "No point inside!!\n";  // should be FATAL (?)
    return;
  }

  double inv_nstep = 1.0 / (double)MCnstep;
  volume = Vbox * int_vol * inv_nstep;

  vol_err = Vbox * sqrt((vol_err * inv_nstep - (int_vol * inv_nstep) * (int_vol * inv_nstep)) * inv_nstep);

  OG = inv_nstep * Vbox * OG;
  if (volume > 1.0e-20) OG = (1.0 / volume) * OG;

  // 3- Monte Carlo integration to compute (1/m) I(G)
  double I11_m = 0.0;
  double I12_m = 0.0;
  double I13_m = 0.0;
  double I22_m = 0.0;
  double I23_m = 0.0;
  double I33_m = 0.0;

  double lx, ly, lz;
  for (size_t i = 0; i < MCnstep; ++i) {
    polyhTool::sobolSequence(3, vv);
    pt3.set(box.min.x + vv[0] * (box.max.x - box.min.x), box.min.y + vv[1] * (box.max.y - box.min.y),
            box.min.z + vv[2] * (box.max.z - box.min.z));
    if (inside(pt3)) {
      lx = pt3.x - OG.x;
      ly = pt3.y - OG.y;
      lz = pt3.z - OG.z;

      I11_m += (ly * ly + lz * lz);
      I12_m -= (lx * ly);
      I13_m -= (lx * lz);
      I22_m += (lx * lx + lz * lz);
      I23_m -= (ly * lz);
      I33_m += (lx * lx + ly * ly);
    }
  }

  double fact = (Vbox / volume) * inv_nstep;
  I11_m *= fact;
  I12_m *= fact;
  I13_m *= fact;
  I22_m *= fact;
  I23_m *= fact;
  I33_m *= fact;

  mat9r I_m(I11_m, I12_m, I13_m, I12_m, I22_m, I23_m, I13_m, I23_m, I33_m);

  mat9r VP;  // Eigen vectors
  vec3r D;   // Eigen values
  I_m.sym_eigen(VP, D);

  // 4- set the precomputed properties
  orientation.set_rot_matrix(VP.c_mtx());
  orientation.normalize();

  quat Qinv = orientation.get_conjugated();
  for (size_t i = 0; i < vertex.size(); ++i) {
    vertex[i] = Qinv * (vertex[i] - OG);
  }
  position = OG;
  inertia_mass = D;
  fitObb();

  // 5- Display some information
  std::cout << "Number of steps in the Monte Carlo integration: " << MCnstep << std::endl;
  std::cout << "      Estimated error for the volume (err/vol): " << vol_err / volume << std::endl;
  std::cout << "                                        Volume: " << volume << std::endl;
  std::cout << "                                   Mass center: " << position << std::endl;
  std::cout << "                                  inertia/mass: " << inertia_mass << std::endl;
  std::cout << "                 Angular position (quaternion): " << orientation << std::endl;
  std::cout << "                                      Position: " << position << std::endl;
  std::cout << std::endl;
}

void Shape::getAABB(AABB& aabb) {
  aabb.set_single(vertex[0]);
  for (size_t i = 1; i < vertex.size(); ++i) {
    aabb.add(vertex[i]);
  }
  aabb.enlarge(radius);
}

void Shape::defineVertexConnectivity() {
  if (ConnectionDefined == true) {
    edgesConnectedWithVertex.clear();
    facesConnectedWithVertex.clear();
  }

  // Edges connected to Vertices
  std::vector<std::set<size_t> > edgeConSet(vertex.size());
  for (size_t e = 0; e < edge.size(); ++e) {
    edgeConSet[edge[e].first].insert(e);
    edgeConSet[edge[e].second].insert(e);
  }
  edgesConnectedWithVertex.resize(vertex.size());
  for (size_t i = 0; i < vertex.size(); ++i) {
    std::copy(edgeConSet[i].begin(), edgeConSet[i].end(), std::back_inserter(edgesConnectedWithVertex[i]));
  }

  // Faces connected to Vertices
  std::vector<std::set<size_t> > faceConSet(vertex.size());
  for (size_t f = 0; f < face.size(); ++f) {
    for (size_t i = 0; i < face[f].size(); ++i) {
      faceConSet[face[f][i]].insert(f);
    }
  }
  facesConnectedWithVertex.resize(vertex.size());
  for (size_t i = 0; i < vertex.size(); ++i) {
    std::copy(faceConSet[i].begin(), faceConSet[i].end(), std::back_inserter(facesConnectedWithVertex[i]));
  }

  ConnectionDefined = true;
}

void Shape::buildOBBtree() {
  tree.clear();

  OBBnode root;
  root.level = 0;
  root.obb = obb;
  for (size_t f = 0; f < face.size(); ++f) root.faceID.push_back(f);
  tree.add(root);

  // int p_beg = 0;
  int p_end = 0;
  for (int level = 1; level <= OBBtreeLevel; level++) {
    int p_beg = p_end;
    p_end = tree.nodes.size();
    int nAdd = 0;
    for (int parent = p_beg; parent < p_end; parent++) {

      int ie[3] = {2, 1, 0};
      if (tree.nodes[parent].obb.extent[ie[0]] < tree.nodes[parent].obb.extent[ie[1]]) {
        ie[0] = 1;
        ie[1] = 0;
      }
      if (tree.nodes[parent].obb.extent[ie[0]] < tree.nodes[parent].obb.extent[ie[2]]) {
        int tmp = ie[0];
        ie[0] = ie[2];
        ie[2] = tmp;
      }
      if (tree.nodes[parent].obb.extent[ie[1]] < tree.nodes[parent].obb.extent[ie[2]]) {
        int tmp = ie[1];
        ie[1] = ie[2];
        ie[2] = tmp;
      }

      for (int ii = 0; ii < 3; ii++) {
        vec3r dir = tree.nodes[parent].obb.e[ie[ii]];
        OBBnode child1, child2;
        child1.level = level;
        child2.level = level;
        for (size_t f = 0; f < tree.nodes[parent].faceID.size(); f++) {
          size_t fid = tree.nodes[parent].faceID[f];
          vec3r p;
          for (size_t v = 0; v < face[fid].size(); v++) p += vertex[face[fid][v]];
          p /= (double)(face[fid].size());
          p -= tree.nodes[parent].obb.center;
          double s = p * dir;
          if (s >= 0.0)
            child2.faceID.push_back(fid);
          else
            child1.faceID.push_back(fid);
        }

        if (!(child1.faceID.empty()) && !(child2.faceID.empty())) {
          updateObb(child1);
          updateObb(child2);
          tree.add(parent, child1);
          tree.add(parent, child2);
          nAdd++;
          break;
        } else {
          child1.faceID.clear();
          child2.faceID.clear();
        }
      }
    }  // recursive-loop over parents

    if (nAdd == 0) break;

  }  // loop over levels
}

/// @brief Define a sub-OBB
/// @see   http://jamesgregson.blogspot.fr/2011/03/latex-test.html
void Shape::updateObb(OBBnode& node) {
  if (node.faceID.empty()) {
    std::cerr << "@updateObb, No face in this OBBnode. Cannot define the "
                 "Oriented Bounding Box"
              << std::endl;
    return;
  }

  std::set<size_t> vertexSet;
  for (size_t f = 0; f < node.faceID.size(); ++f) {
    size_t fid = node.faceID[f];
    for (size_t v = 0; v < face[fid].size(); v++) {
      vertexSet.insert(face[fid][v]);
    }
  }
  std::vector<size_t> vertexID;
  std::copy(vertexSet.begin(), vertexSet.end(), std::back_inserter(vertexID));

  // ==== Build the covariance matrix
  vec3r mu;
  mat9r C;

  // loop over the points to find the mean point
  // location
  for (size_t id = 0; id < vertexID.size(); id++) {
    size_t i = vertexID[id];
    mu += vertex[i];
  }
  mu /= (double)vertexID.size();

  // loop over the points again to build the
  // covariance matrix.  Note that we only have
  // to build terms for the upper trianglular
  // portion since the matrix is symmetric
  double cxx = 0.0, cxy = 0.0, cxz = 0.0, cyy = 0.0, cyz = 0.0, czz = 0.0;
  for (size_t id = 0; id < vertexID.size(); id++) {
    size_t i = vertexID[id];

    vec3r p = vertex[i];
    cxx += p.x * p.x - mu.x * mu.x;
    cxy += p.x * p.y - mu.x * mu.y;
    cxz += p.x * p.z - mu.x * mu.z;
    cyy += p.y * p.y - mu.y * mu.y;
    cyz += p.y * p.z - mu.y * mu.z;
    czz += p.z * p.z - mu.z * mu.z;
  }

  // now build the covariance matrix
  C.xx = cxx;
  C.xy = cxy;
  C.xz = cxz;
  C.yx = cxy;
  C.yy = cyy;
  C.yz = cyz;
  C.zx = cxz;
  C.zy = cyz;
  C.zz = czz;

  // ==== set the OBB parameters from the covariance matrix
  // extract the eigenvalues and eigenvectors from C
  mat9r eigvec;
  vec3r eigval;
  C.sym_eigen(eigvec, eigval);

  // find the right, up and forward vectors from the eigenvectors
  vec3r r(eigvec.xx, eigvec.yx, eigvec.zx);
  vec3r u(eigvec.xy, eigvec.yy, eigvec.zy);
  vec3r f(eigvec.xz, eigvec.yz, eigvec.zz);
  r.normalize();
  u.normalize(), f.normalize();

  // now build the bounding box extents in the rotated frame
  vec3r minim(1e20, 1e20, 1e20), maxim(-1e20, -1e20, -1e20);
  for (size_t id = 0; id < vertexID.size(); id++) {
    size_t i = vertexID[id];
    vec3r p_prime(r * vertex[i], u * vertex[i], f * vertex[i]);
    if (minim.x > p_prime.x) minim.x = p_prime.x;
    if (minim.y > p_prime.y) minim.y = p_prime.y;
    if (minim.z > p_prime.z) minim.z = p_prime.z;
    if (maxim.x < p_prime.x) maxim.x = p_prime.x;
    if (maxim.y < p_prime.y) maxim.y = p_prime.y;
    if (maxim.z < p_prime.z) maxim.z = p_prime.z;
  }

  // set the center of the OBB to be the average of the
  // minimum and maximum, and the extents be half of the
  // difference between the minimum and maximum
  node.obb.center = eigvec * (0.5 * (maxim + minim));
  node.obb.e[0] = r;
  node.obb.e[1] = u;
  node.obb.e[2] = f;
  node.obb.extent = 0.5 * (maxim - minim);

  node.obb.enlarge(radius);  // Add the Minskowski radius
}
