// STATUS: [ ] STABLE  [ ] EXPERIMENTAL  [x] DRAFT

// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

// --------------------------------------------------------------------------------------
// COPY OF toofus/polyhTool.hpp FOR stl2blocks
//
// Modifications:
//   - the class is in the namespace 'stl2blocks', so that it does not clash with
//     the polyhTool of toofus used by Rockable_core (the include guard is changed too)
//   - functions for polyhedra that are NOT CONVEX (closed meshes): cut_mesh (cut by a
//     plane, with caps that can have several loops, holes, and be concave),
//     split_components, offset_inward, orient_faces, close_holes, volume, ...
// The function cut_poly (convex polyhedra) is unchanged.
// --------------------------------------------------------------------------------------

#ifndef STL2BLOCKS_POLYHTOOL_HPP
#define STL2BLOCKS_POLYHTOOL_HPP

#include "geoTool.hpp"
#include "mat9.hpp"
#include "vec2.hpp"
#include "vec3.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

namespace stl2blocks {

class polyhTool {
public:
  // =======================================================================================
  // Basic structures
  // =======================================================================================
  struct plan {
    vec3r pos, normal;
  };

  // A node in a polyhedron
  struct pnode {
    vec3r pos;
    pnode() : pos() {}
  };

  // An edge in a polyhedron
  struct pedge {
    int node0, node1; // index of vector<pnode> in poly
  };

  // A face in a polyhedron
  struct pface {
    std::list<int> nodes; // Index of vector<pnode> in poly
  };

  // polyhedron
  struct poly {
    std::vector<pnode> nodes;
    std::vector<pedge> edges;
    std::vector<pface> faces;
  };
  // =======================================================================================
  // =======================================================================================

  /// @brief Generate a sobol sequence
  /// @param[in] n number of data generated
  /// @param[out] x a vector of data of length n
  static void sobolSequence(const int n, std::vector<double> &x) {
    const int MAXBIT = 30, MAXDIM = 6;
    int j, k, l;
    unsigned int i, im, ipp;
    static int mdeg[MAXDIM] = {1, 2, 3, 3, 4, 4};
    static unsigned int in;
    static std::vector<unsigned int> ix(MAXDIM);
    static std::vector<unsigned int *> iu(MAXBIT);
    static unsigned int ip[MAXDIM]          = {0, 1, 1, 2, 1, 4};
    static unsigned int iv[MAXDIM * MAXBIT] = {1, 1, 1, 1, 1, 1, 3,  1,  3, 3,  1,  1,
                                               5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9};
    static double fac;

    if (n < 0) {
      for (k = 0; k < MAXDIM; k++) { ix[k] = 0; }
      in = 0;
      if (iv[0] != 1) { return; }
      fac = 1.0 / (1 << MAXBIT);
      for (j = 0, k = 0; j < MAXBIT; j++, k += MAXDIM) { iu[j] = &iv[k]; }
      for (k = 0; k < MAXDIM; k++) {
        for (j = 0; j < mdeg[k]; j++) { iu[j][k] <<= (MAXBIT - 1 - j); }
        for (j = mdeg[k]; j < MAXBIT; j++) {
          ipp = ip[k];
          i   = iu[j - mdeg[k]][k];
          i ^= (i >> mdeg[k]);
          for (l = mdeg[k] - 1; l >= 1; l--) {
            if (ipp & 1) { i ^= iu[j - l][k]; }
            ipp >>= 1;
          }
          iu[j][k] = i;
        }
      }
    } else {
      im = in++;
      for (j = 0; j < MAXBIT; j++) {
        if (!(im & 1)) { break; }
        im >>= 1;
      }
      if (j >= MAXBIT) { std::cerr << "MAXBIT too small in sobseq" << std::endl; }
      im       = j * MAXDIM;
      int kmax = (n < MAXDIM) ? n : MAXDIM;
      for (k = 0; k < kmax; k++) {
        ix[k] ^= iv[im + k];
        x[k] = ix[k] * fac;
      }
    }
  }

  // Build a sub_block corresponding to the part of the block
  // which is in the same side than the normal to the plan.
  // Return 1 if the block has been cutted.
  static int cut_poly(const poly &block, const plan &pl, poly &sub_block) {
    // At the beginning sub_block IS block, so we copy it
    sub_block.nodes.clear();
    sub_block.edges.clear();
    sub_block.faces.clear();

    sub_block.nodes = block.nodes;
    sub_block.edges = block.edges;
    sub_block.faces = block.faces;

    // Loop over each edge to intersect with plan
    // then add a node and split the edge in two edges
    std::map<std::pair<int, int>, int> intersected_edges; // associate a pair of nodes with a new node number
    std::pair<int, int> node_pair;
    int node0, node1;
    double alpha;
    size_t initial_node_number =
        sub_block.nodes.size(); // The number of nodes before we add those from the intersections
    size_t new_node_number = initial_node_number;
    std::set<int> new_nodes;
    vec3r A, B;
    pnode new_node;
    pedge new_edge;
    int tmp;
    int nb_cut = 0;

    size_t ne = sub_block.edges.size();
    for (size_t e = 0; e < ne; e++) {
      node0 = sub_block.edges[e].node0;
      node1 = sub_block.edges[e].node1;

      A = sub_block.nodes[node0].pos;
      B = sub_block.nodes[node1].pos;
      geoTool::intersect_edge(A, B, pl.pos, pl.normal, alpha);

      if (alpha > 0.0 && alpha < 1.0) {
        node_pair.first  = node0;
        node_pair.second = node1;
        if (node_pair.second < node_pair.first) { // Take care that first number is lower than the second one
          tmp              = node_pair.first;
          node_pair.first  = node_pair.second;
          node_pair.second = tmp;
        }
        intersected_edges[node_pair] = new_node_number;
        new_nodes.insert(new_node_number);
        new_node.pos = (1.0 - alpha) * A + alpha * B;
        sub_block.nodes.push_back(new_node);

        sub_block.edges[e].node1 = new_node_number;
        new_edge.node0           = new_node_number;
        new_edge.node1           = node1;
        sub_block.edges.push_back(new_edge);

        new_node_number++;
        nb_cut++;
      }
    }

    if (nb_cut < 3) { return 0; } // Because, we need 3 points (at least) to create the face that close the shape

    // Insert the new nodes in the faces
    std::list<int>::iterator it0, it1;
    std::map<std::pair<int, int>, int>::iterator mp;
    size_t nf = sub_block.faces.size();
    for (size_t f = 0; f < nf; f++) {
      for (it0 = sub_block.faces[f].nodes.begin(); it0 != sub_block.faces[f].nodes.end(); ++it0) {
        it1 = it0;
        ++it1;
        if (it1 == sub_block.faces[f].nodes.end()) it1 = sub_block.faces[f].nodes.begin();
        node0 = *it0;
        node1 = *it1;
        if (node1 < node0) {
          tmp   = node0;
          node0 = node1;
          node1 = tmp;
        }
        node_pair.first  = node0;
        node_pair.second = node1;

        mp = intersected_edges.find(node_pair);
        if (mp != intersected_edges.end()) { // If the pair has been found
          sub_block.faces[f].nodes.insert(it1, mp->second);
        }
      }
    }

    // Identify the nodes to be removed
    std::set<int> nodes_to_remove;
    std::set<int>::iterator inode;
    for (size_t n = 0; n < initial_node_number; n++) {
      if ((sub_block.nodes[n].pos - pl.pos) * pl.normal < 0.0) { nodes_to_remove.insert(n); }
    }

    // Remove node "to be removed" in each face
    for (size_t f = 0; f < sub_block.faces.size(); f++) {
      std::list<int>::iterator it = sub_block.faces[f].nodes.begin();
      size_t nb = sub_block.faces[f].nodes.size(), count = 0;
      while (count <= nb && it != sub_block.faces[f].nodes.end()) {
        inode = nodes_to_remove.find(*it);
        if (inode != nodes_to_remove.end()) // If the node needs to be removed
        {
          it = sub_block.faces[f].nodes.erase(it);
        } else {
          ++it;
        }
        ++count;
      }
    }

    // Add a face that close the sub-poly
    // together with the corresponding edges
    std::list<pedge> edge_loop0;
    int num_search = 0;
    pedge E;
    E.node0 = E.node1 = -1;
    for (size_t f = 0; f < sub_block.faces.size(); f++) {
      num_search = 0;
      E.node0 = E.node1 = -1;
      size_t count      = 0;
      for (std::list<int>::iterator it = sub_block.faces[f].nodes.begin(); it != sub_block.faces[f].nodes.end(); ++it) {
        if (*it >= (int)initial_node_number) { // It means that this is a node that has been added
          if (num_search == 0) {
            E.node0    = *it;
            num_search = 1;
          } else {
            E.node1 = *it;
            edge_loop0.push_back(E);
            E.node0 = E.node1;
          }
          ++count;
        }
      }
      if (count == 1) { std::cout << "1 seul nouveau noeud dans une face\n"; }
      if (count > 2) { std::cout << "plus de 2 nouveaux noeuds dans une face\n"; }
    }

    if (edge_loop0.size() != (size_t)nb_cut) {
      std::cout << "edge_loop0.size() = " << edge_loop0.size() << std::endl;
      std::cout << "nb_cut = " << nb_cut << std::endl;
    }

    for (std::list<pedge>::iterator ite = edge_loop0.begin(); ite != edge_loop0.end(); ++ite) {
      sub_block.edges.push_back(*ite);
    }

    // Assemble new edges to make a loop
    std::list<pedge>::iterator ite = edge_loop0.begin();
    std::vector<pedge> edge_loop1;
    pedge current_edge = *ite;
    edge_loop1.push_back(current_edge);
    edge_loop0.erase(ite);
    ite = edge_loop0.begin();

    int swap;
    while (!edge_loop0.empty()) {
      if (ite == edge_loop0.end()) { break; }
      for (ite = edge_loop0.begin(); ite != edge_loop0.end(); ++ite) {
        if (current_edge.node1 == ite->node1) {
          swap         = ite->node0;
          ite->node0   = ite->node1;
          ite->node1   = swap;
          current_edge = *ite;
          edge_loop1.push_back(current_edge);
          edge_loop0.erase(ite);
          ite = edge_loop0.begin();
          break;
        } else if (current_edge.node1 == ite->node0) {
          current_edge = *ite;
          edge_loop1.push_back(current_edge);
          edge_loop0.erase(ite);
          ite = edge_loop0.begin();
          break;
        }
      }
    }

    /*
    if (edge_loop1.size() != (size_t)nb_cut) {
    std::cout << "edge_loop1.size() = " << edge_loop1.size() << std::endl;
    std::cout << "nb_cut = " << nb_cut << std::endl;
    }
    */

    pface F;
    F.nodes.clear();
    for (size_t e = 0; e < edge_loop1.size(); e++) { F.nodes.push_back(edge_loop1[e].node0); }
    sub_block.faces.push_back(F);

    // Remove faces with less than 3 nodes
    std::vector<pface> faces_copy;
    for (size_t f = 0; f < sub_block.faces.size(); f++) {
      if (sub_block.faces[f].nodes.size() >= 3) { faces_copy.push_back(sub_block.faces[f]); }
    }
    sub_block.faces.clear();
    sub_block.faces = faces_copy;

    // Remove edges with less than 2 nodes
    std::vector<pedge> edges_copy;
    for (size_t e = 0; e < sub_block.edges.size(); e++) {
      if (nodes_to_remove.find(sub_block.edges[e].node0) == nodes_to_remove.end() &&
          nodes_to_remove.find(sub_block.edges[e].node1) == nodes_to_remove.end()) {
        edges_copy.push_back(sub_block.edges[e]);
      }
    }
    sub_block.edges.clear();
    sub_block.edges = edges_copy;

    // ====== Remove and re-number the nodes
    std::vector<int> new_number(sub_block.nodes.size());
    std::vector<pnode> nodes_copy;
    int NN = 0;
    for (size_t n = 0; n < sub_block.nodes.size(); n++) {
      inode = nodes_to_remove.find(n);
      if (inode == nodes_to_remove.end()) { // ... has not been removed
        new_number[n] = NN++;
        nodes_copy.push_back(sub_block.nodes[n]);
      }
    }
    sub_block.nodes.clear();
    sub_block.nodes = nodes_copy;

    // Re-numbering the nodes in faces
    for (size_t f = 0; f < sub_block.faces.size(); f++) {
      for (std::list<int>::iterator it = sub_block.faces[f].nodes.begin(); it != sub_block.faces[f].nodes.end(); ++it) {
        (*it) = new_number[*it];
      }
    }

    // Re-numbering the nodes in edges
    for (size_t e = 0; e < sub_block.edges.size(); e++) {
      sub_block.edges[e].node0 = new_number[sub_block.edges[e].node0];
      sub_block.edges[e].node1 = new_number[sub_block.edges[e].node1];
    }

    return 1;
  }

  // =======================================================================================
  // Polyhedra that are NOT CONVEX (closed meshes)
  // =======================================================================================
  // The faces are convex polygons (triangles, mostly), with their nodes in counterclockwise
  // order seen from outside the polyhedron.

  /// @brief Unnormalized normal of a face (its length is twice the area), Newell's method
  static vec3r face_area_vector(const poly &P, const pface &F) {
    vec3r N;
    for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
      auto jt = std::next(it);
      if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
      N += cross(P.nodes[*it].pos, P.nodes[*jt].pos);
    }
    return N;
  }

  /// @brief Volume enclosed by a closed polyhedron (positive when the faces are oriented outwards)
  static double volume(const poly &P) {
    double V = 0.0;
    for (const pface &F : P.faces) {
      if (F.nodes.size() < 3) { continue; }
      auto it      = F.nodes.begin();
      const vec3r &a = P.nodes[*it].pos;
      ++it;
      auto jt = std::next(it);
      for (; jt != F.nodes.end(); ++it, ++jt) { V += a * cross(P.nodes[*it].pos, P.nodes[*jt].pos); }
    }
    return V / 6.0;
  }

  /// @brief Size of the bounding box (its diagonal)
  static double size(const poly &P) {
    if (P.nodes.empty()) { return 0.0; }
    vec3r mn = P.nodes[0].pos, mx = P.nodes[0].pos;
    for (const pnode &n : P.nodes) {
      for (int c = 0; c < 3; c++) {
        mn[c] = std::min(mn[c], n.pos[c]);
        mx[c] = std::max(mx[c], n.pos[c]);
      }
    }
    return (mx - mn).length();
  }

  /// @brief Rebuild the edges from the faces (each side of a face, once)
  static void rebuild_edges(poly &P) {
    std::set<std::pair<int, int>> S;
    for (const pface &F : P.faces) {
      for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
        auto jt = std::next(it);
        if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
        int a = *it, b = *jt;
        if (a > b) { std::swap(a, b); }
        if (a != b) { S.insert(std::make_pair(a, b)); }
      }
    }
    P.edges.clear();
    for (const auto &p : S) {
      pedge E;
      E.node0 = p.first;
      E.node1 = p.second;
      P.edges.push_back(E);
    }
  }

  /// @brief Remove the nodes that no face uses, renumber the others, and rebuild the edges
  static void remove_unused_nodes(poly &P) {
    std::vector<int> num(P.nodes.size(), -1);
    for (const pface &F : P.faces) {
      for (int n : F.nodes) { num[n] = 0; }
    }
    std::vector<pnode> nodes;
    for (size_t n = 0; n < P.nodes.size(); n++) {
      if (num[n] == 0) {
        num[n] = (int)nodes.size();
        nodes.push_back(P.nodes[n]);
      }
    }
    P.nodes = nodes;
    for (pface &F : P.faces) {
      for (int &n : F.nodes) { n = num[n]; }
    }
    rebuild_edges(P);
  }

  /// @brief Make the orientation of the faces consistent (a face and its neighbour run along
  ///        their common side in opposite directions), then outwards (positive volume).
  /// @return the number of faces that have been reversed
  static int orient_faces(poly &P) {
    std::map<std::pair<int, int>, std::vector<int>> sideFaces; // unordered side -> faces
    for (size_t f = 0; f < P.faces.size(); f++) {
      const pface &F = P.faces[f];
      for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
        auto jt = std::next(it);
        if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
        sideFaces[std::minmax(*it, *jt)].push_back((int)f);
      }
    }
    // Does face f run from a to b?
    auto runs = [&P](int f, int a, int b) {
      const pface &F = P.faces[f];
      for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
        auto jt = std::next(it);
        if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
        if (*it == a && *jt == b) { return true; }
      }
      return false;
    };

    int nbReversed = 0;
    std::vector<int> done(P.faces.size(), 0);
    for (size_t f0 = 0; f0 < P.faces.size(); f0++) {
      if (done[f0]) { continue; }
      done[f0] = 1;
      std::vector<int> stack(1, (int)f0);
      while (!stack.empty()) {
        int f = stack.back();
        stack.pop_back();
        std::vector<int> nodes(P.faces[f].nodes.begin(), P.faces[f].nodes.end());
        for (size_t k = 0; k < nodes.size(); k++) {
          int a = nodes[k], b = nodes[(k + 1) % nodes.size()];
          for (int g : sideFaces[std::minmax(a, b)]) {
            if (g == f || done[g]) { continue; }
            if (runs(g, a, b)) { // same direction: g must be reversed
              P.faces[g].nodes.reverse();
              nbReversed++;
            }
            done[g] = 1;
            stack.push_back(g);
          }
        }
      }
    }

    if (volume(P) < 0.0) {
      for (pface &F : P.faces) { F.nodes.reverse(); }
      nbReversed = (int)P.faces.size() - nbReversed;
    }
    return nbReversed;
  }

  /// @brief The number of sides of faces that are not run the other way by another face
  ///        (0 for a closed polyhedron with consistently oriented faces)
  static int open_sides(const poly &P) {
    std::map<std::pair<int, int>, int> count;
    for (const pface &F : P.faces) {
      for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
        auto jt = std::next(it);
        if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
        count[std::make_pair(*it, *jt)]++;
      }
    }
    int nb = 0;
    for (const auto &c : count) {
      auto r = count.find(std::make_pair(c.first.second, c.first.first));
      if (c.second != 1 || r == count.end() || r->second != 1) { nb++; }
    }
    return nb;
  }

  /// @brief Remove the dangling faces: a face that has a free side (a side of this face only)
  ///        and a side shared by more than two faces is a fin stuck on the mesh.
  /// @return the number of faces removed
  static int remove_dangling_faces(poly &P) {
    int nbRemoved = 0;
    bool removed  = true;
    while (removed) {
      removed = false;
      std::map<std::pair<int, int>, int> count;
      for (const pface &F : P.faces) {
        for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
          auto jt = std::next(it);
          if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
          count[std::minmax(*it, *jt)]++;
        }
      }
      std::vector<pface> kept;
      for (const pface &F : P.faces) {
        bool freeSide = false, sharedSide = false;
        for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
          auto jt = std::next(it);
          if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
          int c = count[std::minmax(*it, *jt)];
          if (c == 1) { freeSide = true; }
          if (c > 2) { sharedSide = true; }
        }
        if (freeSide && sharedSide) {
          nbRemoved++;
          removed = true;
        } else {
          kept.push_back(F);
        }
      }
      P.faces = kept;
    }
    if (nbRemoved > 0) { remove_unused_nodes(P); }
    return nbRemoved;
  }

  /// @brief Close the holes of a mesh: each loop of free sides (sides of a single face) is
  ///        closed by a fan of triangles around a new node, at the centre of the loop.
  ///        The faces must be consistently oriented (orient_faces).
  /// @param[out] nbNonManifold the number of sides shared by more than two faces
  /// @return the number of holes closed
  static int close_holes(poly &P, int &nbNonManifold) {
    std::map<std::pair<int, int>, int> count;
    std::map<int, int> next; // free side a->b of a face: the hole runs from b to a
    for (const pface &F : P.faces) {
      for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
        auto jt = std::next(it);
        if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
        count[std::minmax(*it, *jt)]++;
      }
    }
    nbNonManifold = 0;
    for (const auto &c : count) {
      if (c.second > 2) { nbNonManifold++; }
    }
    for (const pface &F : P.faces) {
      for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
        auto jt = std::next(it);
        if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
        if (count[std::minmax(*it, *jt)] == 1) { next[*jt] = *it; }
      }
    }

    int nbHoles = 0;
    while (!next.empty()) {
      std::vector<int> loop;
      int start = next.begin()->first;
      int cur   = start;
      while (next.count(cur)) {
        loop.push_back(cur);
        int nxt = next[cur];
        next.erase(cur);
        cur = nxt;
        if (cur == start) { break; }
      }
      if (cur != start || loop.size() < 3) { continue; } // not a closed loop
      pnode C;
      for (int n : loop) { C.pos += P.nodes[n].pos; }
      C.pos /= (double)loop.size();
      int c = (int)P.nodes.size();
      P.nodes.push_back(C);
      for (size_t k = 0; k < loop.size(); k++) {
        pface F;
        F.nodes = {loop[k], loop[(k + 1) % loop.size()], c};
        P.faces.push_back(F);
      }
      nbHoles++;
    }
    rebuild_edges(P);
    return nbHoles;
  }

  /// @brief Move each node inwards so that the planes of the faces around it move by R
  ///        (least squares: each node is moved by d with d.n = -R for the normals n of the
  ///        faces around it, weighted by the angle of the face at the node).
  /// @param[out] maxError the largest distance between a node and the moved plane of one of
  ///             its faces (relative to R)
  /// @param[out] nbFlipped the number of faces that are still turned over (parts thinner than 2R)
  /// @param[out] nbClamped the number of nodes whose move was limited to 4R
  /// @param[out] nbReduced the number of nodes moved less, so that no face is turned over
  static void offset_inward(poly &P, double R, double &maxError, int &nbFlipped, int &nbClamped, int &nbReduced) {
    std::vector<vec3r> unitNormal(P.faces.size());
    for (size_t f = 0; f < P.faces.size(); f++) { unitNormal[f] = face_area_vector(P, P.faces[f]).normalized(); }

    std::vector<mat9r> A(P.nodes.size());
    std::vector<vec3r> b(P.nodes.size());
    for (size_t f = 0; f < P.faces.size(); f++) {
      std::vector<int> v(P.faces[f].nodes.begin(), P.faces[f].nodes.end());
      const vec3r &n = unitNormal[f];
      for (size_t k = 0; k < v.size(); k++) {
        vec3r e1 = P.nodes[v[(k + 1) % v.size()]].pos - P.nodes[v[k]].pos;
        vec3r e2 = P.nodes[v[(k + v.size() - 1) % v.size()]].pos - P.nodes[v[k]].pos;
        double l = e1.length() * e2.length();
        if (l <= 0.0) { continue; }
        double w = acos(std::max(-1.0, std::min(1.0, (e1 * e2) / l))); // angle at the node
        mat9r &M = A[v[k]];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) { M[3 * i + j] += w * n[i] * n[j]; }
        }
        b[v[k]] -= w * R * n;
      }
    }

    std::vector<vec3r> d(P.nodes.size());
    nbClamped = 0;
    for (size_t i = 0; i < P.nodes.size(); i++) {
      mat9r V;
      vec3r D;
      A[i].sym_eigen(V, D);
      double Dmax = std::max(std::fabs(D[0]), std::max(std::fabs(D[1]), std::fabs(D[2])));
      for (int k = 0; k < 3; k++) {
        if (std::fabs(D[k]) <= 1e-4 * Dmax) { continue; } // pseudo-inverse
        vec3r e(V[k], V[3 + k], V[6 + k]);
        d[i] += ((e * b[i]) / D[k]) * e;
      }
      double l = d[i].length();
      if (l > 4.0 * R) {
        d[i] *= 4.0 * R / l;
        nbClamped++;
      }
    }

    // The moves of the nodes of a face that would be turned over (a face smaller than R,
    // or a part thinner than 2R) are reduced until no face is turned over
    std::vector<pnode> initial = P.nodes;
    nbReduced = 0;
    std::set<int> reduced;
    for (int iter = 0; iter < 100; iter++) {
      for (size_t i = 0; i < P.nodes.size(); i++) { P.nodes[i].pos = initial[i].pos + d[i]; }
      std::set<int> toReduce;
      for (size_t f = 0; f < P.faces.size(); f++) {
        if (face_area_vector(P, P.faces[f]) * unitNormal[f] <= 0.0) {
          for (int n : P.faces[f].nodes) { toReduce.insert(n); }
        }
      }
      if (toReduce.empty()) { break; }
      for (int n : toReduce) {
        d[n] *= 0.8;
        reduced.insert(n);
      }
    }
    nbReduced = (int)reduced.size();

    maxError = 0.0;
    for (size_t f = 0; f < P.faces.size(); f++) {
      for (int n : P.faces[f].nodes) {
        maxError = std::max(maxError, std::fabs(d[n] * unitNormal[f] + R) / R);
      }
    }
    nbFlipped = 0;
    for (size_t f = 0; f < P.faces.size(); f++) {
      if (face_area_vector(P, P.faces[f]) * unitNormal[f] <= 0.0) { nbFlipped++; }
    }
  }

  /// @brief Simplify the faces of a closed polyhedron, once the cuts are done: two
  ///        neighbouring faces in the same plane are merged when the merged face is convex,
  ///        and the nodes in the middle of a straight side of a face are removed from it (a
  ///        node that no face uses any more is removed). Each face is then a convex polygon
  ///        without three aligned nodes, as needed by the tests of Rockable.
  static void simplify_faces(poly &P) {
    double L       = size(P);
    double epsDist = 1e-9 * L;

    // Merge coplanar neighbours
    bool merged = true;
    while (merged) {
      merged = false;
      std::map<std::pair<int, int>, int> faceOfSide; // directed side -> face
      for (size_t f = 0; f < P.faces.size(); f++) {
        const pface &F = P.faces[f];
        for (auto it = F.nodes.begin(); it != F.nodes.end(); ++it) {
          auto jt = std::next(it);
          if (jt == F.nodes.end()) { jt = F.nodes.begin(); }
          faceOfSide[std::make_pair(*it, *jt)] = (int)f;
        }
      }
      for (size_t f = 0; f < P.faces.size() && !merged; f++) {
        vec3r nf = face_area_vector(P, P.faces[f]);
        nf.normalize();
        std::vector<int> a(P.faces[f].nodes.begin(), P.faces[f].nodes.end());
        for (size_t k = 0; k < a.size() && !merged; k++) {
          int n0 = a[k], n1 = a[(k + 1) % a.size()];
          auto g = faceOfSide.find(std::make_pair(n1, n0));
          if (g == faceOfSide.end() || g->second == (int)f) { continue; }
          const pface &G = P.faces[g->second];
          vec3r ng       = face_area_vector(P, G);
          ng.normalize();
          if (nf * ng < 1.0 - 1e-10) { continue; }
          bool inPlane = true;
          for (int n : G.nodes) {
            if (std::fabs((P.nodes[n].pos - P.nodes[n0].pos) * nf) > epsDist) { inPlane = false; }
          }
          if (!inPlane) { continue; }
          std::vector<int> b(G.nodes.begin(), G.nodes.end());
          int shared = 0;
          std::set<int> sa(a.begin(), a.end());
          for (size_t j = 0; j < b.size(); j++) {
            if (sa.count(b[j]) && sa.count(b[(j + 1) % b.size()])) { shared++; }
          }
          if (shared != 1) { continue; }
          // Merged loop: a from n1 round to n0, then b from n0 round to n1 (both excluded)
          std::vector<int> m;
          for (size_t j = 0; j < a.size(); j++) { m.push_back(a[(k + 1 + j) % a.size()]); }
          size_t jb = std::find(b.begin(), b.end(), n0) - b.begin();
          for (size_t j = 1; j + 1 < b.size(); j++) { m.push_back(b[(jb + j) % b.size()]); }
          bool convex = true;
          for (size_t j = 0; j < m.size() && convex; j++) {
            const vec3r &p0 = P.nodes[m[(j + m.size() - 1) % m.size()]].pos;
            const vec3r &p1 = P.nodes[m[j]].pos;
            const vec3r &p2 = P.nodes[m[(j + 1) % m.size()]].pos;
            if (cross(p1 - p0, p2 - p1) * nf < -epsDist * (p1 - p0).length()) { convex = false; }
          }
          if (!convex) { continue; }
          P.faces[f].nodes.assign(m.begin(), m.end());
          P.faces.erase(P.faces.begin() + g->second);
          merged = true;
        }
      }
    }

    // Remove the nodes in the middle of a straight side
    for (pface &F : P.faces) {
      bool removed = true;
      while (removed && F.nodes.size() > 3) {
        removed = false;
        std::vector<int> v(F.nodes.begin(), F.nodes.end());
        for (size_t k = 0; k < v.size(); k++) {
          const vec3r &p0 = P.nodes[v[(k + v.size() - 1) % v.size()]].pos;
          const vec3r &p1 = P.nodes[v[k]].pos;
          const vec3r &p2 = P.nodes[v[(k + 1) % v.size()]].pos;
          vec3r u = p1 - p0, w = p2 - p1;
          if (cross(u, w).length() <= 1e-9 * u.length() * w.length() && u * w > 0.0) {
            v.erase(v.begin() + k);
            F.nodes.assign(v.begin(), v.end());
            removed = true;
            break;
          }
        }
      }
    }
    remove_unused_nodes(P);
  }

  /// @brief Closest point to p on the triangle abc (Ericson, Real-Time Collision Detection)
  static vec3r closest_on_triangle(const vec3r &p, const vec3r &a, const vec3r &b, const vec3r &c) {
    vec3r ab = b - a, ac = c - a, ap = p - a;
    double d1 = ab * ap, d2 = ac * ap;
    if (d1 <= 0.0 && d2 <= 0.0) { return a; }
    vec3r bp  = p - b;
    double d3 = ab * bp, d4 = ac * bp;
    if (d3 >= 0.0 && d4 <= d3) { return b; }
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) { return a + (d1 / (d1 - d3)) * ab; }
    vec3r cp  = p - c;
    double d5 = ab * cp, d6 = ac * cp;
    if (d6 >= 0.0 && d5 <= d6) { return c; }
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) { return a + (d2 / (d2 - d6)) * ac; }
    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) { return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b); }
    double denom = 1.0 / (va + vb + vc);
    return a + (vb * denom) * ab + (vc * denom) * ac;
  }

  /// @brief Shortest distance between the segments [p1, q1] and [p2, q2]
  static double segment_distance(const vec3r &p1, const vec3r &q1, const vec3r &p2, const vec3r &q2) {
    vec3r d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    double a = d1 * d1, e = d2 * d2, f = d2 * r;
    double s = 0.0, t = 0.0;
    if (a <= 1e-300 && e <= 1e-300) { return r.length(); }
    if (a <= 1e-300) {
      t = std::max(0.0, std::min(1.0, f / e));
    } else {
      double c = d1 * r;
      if (e <= 1e-300) {
        s = std::max(0.0, std::min(1.0, -c / a));
      } else {
        double b     = d1 * d2;
        double denom = a * e - b * b;
        s            = (denom > 0.0) ? std::max(0.0, std::min(1.0, (b * f - c * e) / denom)) : 0.0;
        t            = (b * s + f) / e;
        if (t < 0.0) {
          t = 0.0;
          s = std::max(0.0, std::min(1.0, -c / a));
        } else if (t > 1.0) {
          t = 1.0;
          s = std::max(0.0, std::min(1.0, (b - c) / a));
        }
      }
    }
    return ((p1 + s * d1) - (p2 + t * d2)).length();
  }

  /// @brief Shortest distance between the surfaces of two polyhedra (that do not cross each
  ///        other): the nodes of each one to the faces of the other, and the edges to the edges
  static double distance(const poly &A, const poly &B) {
    double dmin = 1e300;
    for (int pass = 0; pass < 2; pass++) {
      const poly &P = (pass == 0) ? A : B;
      const poly &Q = (pass == 0) ? B : A;
      for (const pface &F : Q.faces) {
        std::vector<int> v(F.nodes.begin(), F.nodes.end());
        for (size_t k = 1; k + 1 < v.size(); k++) {
          const vec3r &a = Q.nodes[v[0]].pos, &b = Q.nodes[v[k]].pos, &c = Q.nodes[v[k + 1]].pos;
          for (const pnode &N : P.nodes) { dmin = std::min(dmin, (N.pos - closest_on_triangle(N.pos, a, b, c)).length()); }
        }
      }
    }
    for (const pedge &E : A.edges) {
      for (const pedge &G : B.edges) {
        dmin = std::min(dmin, segment_distance(A.nodes[E.node0].pos, A.nodes[E.node1].pos, B.nodes[G.node0].pos,
                                               B.nodes[G.node1].pos));
      }
    }
    return dmin;
  }

  /// @brief Put two polyhedra together in one (that is not connected)
  static void append(poly &A, const poly &B) {
    int offset = (int)A.nodes.size();
    A.nodes.insert(A.nodes.end(), B.nodes.begin(), B.nodes.end());
    for (pedge E : B.edges) {
      E.node0 += offset;
      E.node1 += offset;
      A.edges.push_back(E);
    }
    for (pface F : B.faces) {
      for (int &n : F.nodes) { n += offset; }
      A.faces.push_back(F);
    }
  }

  /// @brief Separate the parts of a polyhedron that are not connected (by their nodes)
  static void split_components(const poly &P, std::vector<poly> &parts) {
    std::vector<int> root(P.nodes.size());
    std::iota(root.begin(), root.end(), 0);
    std::function<int(int)> find = [&](int i) { return root[i] == i ? i : (root[i] = find(root[i])); };
    for (const pface &F : P.faces) {
      int r0 = find(F.nodes.front());
      for (int n : F.nodes) { root[find(n)] = r0; }
    }
    std::map<int, size_t> partOfRoot;
    parts.clear();
    for (const pface &F : P.faces) {
      int r = find(F.nodes.front());
      auto it = partOfRoot.find(r);
      if (it == partOfRoot.end()) {
        it = partOfRoot.insert(std::make_pair(r, parts.size())).first;
        parts.push_back(poly());
        parts.back().nodes = P.nodes;
      }
      parts[it->second].faces.push_back(F);
    }
    for (poly &part : parts) { remove_unused_nodes(part); }
  }

  /// @brief Twice the signed area of a triangle in 2D
  static double cross2(const vec2r &a, const vec2r &b, const vec2r &c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  }

  /// @brief Is p inside the polygon (crossing number)?
  static bool inside2(const vec2r &p, const std::vector<vec2r> &poly2) {
    bool odd = false;
    for (size_t i = 0, j = poly2.size() - 1; i < poly2.size(); j = i++) {
      const vec2r &a = poly2[i], &b = poly2[j];
      if ((a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x)) { odd = !odd; }
    }
    return odd;
  }

  /// @brief Triangulate a polygon (counterclockwise), possibly with holes (clockwise), by
  ///        ear clipping. The holes are first joined to the outer loop by bridges.
  ///        No node is skipped: a node in the middle of a straight side is kept in the
  ///        triangles, so that the triangles match the faces around the polygon.
  /// @param[in] pos2 the 2D positions of the nodes (indexed by node number)
  /// @param[out] triangles the triangles, as node numbers
  /// @return false if the polygon is not simple: an ear had to be forced, or a triangle is
  ///         turned over (it is kept anyway, so that the polyhedron stays closed)
  static bool triangulate(const std::vector<int> &outer, std::vector<std::vector<int>> holes,
                          const std::map<int, vec2r> &pos2, std::vector<std::array<int, 3>> &triangles,
                          double eps) {
    std::vector<int> P = outer;
    auto X               = [&pos2](int n) -> const vec2r & { return pos2.at(n); };

    // The holes with the rightmost node first
    auto maxX = [&X](const std::vector<int> &h) {
      double m = -1e300;
      for (int n : h) { m = std::max(m, X(n).x); }
      return m;
    };
    std::sort(holes.begin(), holes.end(),
              [&maxX](const std::vector<int> &a, const std::vector<int> &b) { return maxX(a) > maxX(b); });

    for (const std::vector<int> &H : holes) {
      // M, the rightmost node of the hole
      size_t im = 0;
      for (size_t k = 1; k < H.size(); k++) {
        if (X(H[k]).x > X(H[im]).x) { im = k; }
      }
      vec2r M = X(H[im]);
      // The closest side of P crossed by the ray from M towards +x
      double bestX = 1e300;
      int bestK    = -1;
      for (size_t k = 0; k < P.size(); k++) {
        const vec2r &a = X(P[k]), &b = X(P[(k + 1) % P.size()]);
        if ((a.y - M.y) * (b.y - M.y) > 0.0 || a.y == b.y) { continue; }
        double x = a.x + (M.y - a.y) / (b.y - a.y) * (b.x - a.x);
        if (x >= M.x && x < bestX) {
          bestX = x;
          bestK = (int)k;
        }
      }
      if (bestK < 0) { return false; }
      vec2r I(bestX, M.y);
      // The visible node: the end of that side with the largest x...
      size_t kp = (X(P[bestK]).x > X(P[(bestK + 1) % P.size()]).x) ? bestK : (bestK + 1) % P.size();
      // ... unless a reflex node of P is inside the triangle M, I, P[kp]: take the one closest
      // in angle to the ray
      vec2r Pv     = X(P[kp]);
      double best  = 1e300;
      size_t kbest = kp;
      for (size_t k = 0; k < P.size(); k++) {
        if (k == kp) { continue; }
        const vec2r &q = X(P[k]);
        const vec2r &a = X(P[(k + P.size() - 1) % P.size()]), &c = X(P[(k + 1) % P.size()]);
        if (cross2(a, q, c) > 0.0) { continue; } // not reflex
        bool in = false;
        double s1 = cross2(M, I, q), s2 = cross2(I, Pv, q), s3 = cross2(Pv, M, q);
        if ((s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0)) { in = true; }
        if (!in || q.x < M.x) { continue; }
        double ang = std::fabs(q.y - M.y) / std::max(q.x - M.x, 1e-300);
        if (ang < best) {
          best  = ang;
          kbest = k;
        }
      }
      if (best < 1e300) { kp = kbest; }
      // Splice: ... P[kp], M, (the hole from M, round to M), M, P[kp], ...
      std::vector<int> Q(P.begin(), P.begin() + kp + 1);
      for (size_t k = 0; k <= H.size(); k++) { Q.push_back(H[(im + k) % H.size()]); }
      Q.push_back(P[kp]);
      Q.insert(Q.end(), P.begin() + kp + 1, P.end());
      P = Q;
    }

    // Ear clipping
    bool forced = false;
    while (P.size() > 3) {
      size_t m     = P.size();
      int earK     = -1;
      double bestC = -1e300;
      int bestK    = 0;
      for (size_t k = 0; k < m; k++) {
        int ia = P[(k + m - 1) % m], ib = P[k], ic = P[(k + 1) % m];
        const vec2r &a = X(ia), &b = X(ib), &c = X(ic);
        double cr = cross2(a, b, c);
        if (cr > bestC) {
          bestC = cr;
          bestK = (int)k;
        }
        if (cr <= eps) { continue; } // reflex or flat
        bool ear = true;
        for (size_t j = 0; j < m && ear; j++) {
          int n = P[j];
          if (n == ia || n == ib || n == ic) { continue; }
          const vec2r &q = X(n);
          if (cross2(a, b, q) >= -eps && cross2(b, c, q) >= -eps && cross2(c, a, q) >= -eps) { ear = false; }
        }
        if (ear) {
          earK = (int)k;
          break;
        }
      }
      if (earK < 0) {
        earK   = bestK;
        forced = true;
      }
      int ia = P[(earK + m - 1) % m], ib = P[earK], ic = P[(earK + 1) % m];
      double cr = cross2(X(ia), X(ib), X(ic));
      if (cr != 0.0) { triangles.push_back({ia, ib, ic}); }
      if (cr < 0.0) { forced = true; }
      P.erase(P.begin() + earK);
    }
    if (P.size() == 3) {
      double cr = cross2(X(P[0]), X(P[1]), X(P[2]));
      if (cr != 0.0) { triangles.push_back({P[0], P[1], P[2]}); }
      if (cr < 0.0) { forced = true; }
    }
    return !forced;
  }

  /// @brief Keep the part of a closed polyhedron (convex or not) that is on the side of the
  ///        normal of the plane, and close it with faces in the plane (the caps). A cap can
  ///        have several loops (several parts), holes, and be concave; a cap that is convex
  ///        and has no hole is a single face, the others are triangulated.
  /// @param[out] nbBad the number of cap loops that could not be closed, or triangulated
  /// @return 0 if the plane does not cut the polyhedron (all is kept), 1 if it is cut, and 2
  ///         if nothing is kept
  static int cut_mesh(const poly &block, const plan &pl, poly &sub_block, int &nbBad) {
    nbBad   = 0;
    vec3r n = pl.normal;
    n.normalize();
    double L   = size(block);
    double tol = 1e-10 * L;

    size_t nn = block.nodes.size();
    std::vector<double> d(nn);
    size_t nbIn = 0;
    for (size_t i = 0; i < nn; i++) {
      d[i] = (block.nodes[i].pos - pl.pos) * n;
      if (std::fabs(d[i]) < tol) { d[i] = 0.0; } // a node in the plane is kept
      if (d[i] >= 0.0) { nbIn++; }
    }
    if (nbIn == nn) {
      sub_block = block;
      return 0;
    }
    sub_block.nodes.clear();
    sub_block.edges.clear();
    sub_block.faces.clear();
    if (nbIn == 0) { return 2; }

    sub_block.nodes = block.nodes;
    std::map<std::pair<int, int>, int> cutNode;
    // The node where the side from a (kept) to b (removed) crosses the plane
    auto crossing = [&](int a, int b) -> int {
      if (d[a] == 0.0) { return a; }
      auto key = std::minmax(a, b);
      auto it  = cutNode.find(key);
      if (it != cutNode.end()) { return it->second; }
      double alpha = d[a] / (d[a] - d[b]);
      pnode N;
      N.pos = (1.0 - alpha) * block.nodes[a].pos + alpha * block.nodes[b].pos;
      int c = (int)sub_block.nodes.size();
      sub_block.nodes.push_back(N);
      cutNode[key] = c;
      return c;
    };

    // Clip each face (Sutherland-Hodgman); its side in the plane, from the node where it
    // leaves the kept part to the node where it comes back, is run the other way by the cap
    std::multimap<int, int> capSide;
    for (const pface &F : block.faces) {
      std::vector<int> v(F.nodes.begin(), F.nodes.end());
      size_t m  = v.size();
      bool all  = true;
      bool none = true;
      for (int k : v) {
        if (d[k] >= 0.0) {
          none = false;
        } else {
          all = false;
        }
      }
      if (all) {
        sub_block.faces.push_back(F);
        continue;
      }
      if (none) { continue; }
      std::vector<int> out;
      int exitNode = -1, entryNode = -1;
      for (size_t k = 0; k < m; k++) {
        int a = v[k], b = v[(k + 1) % m];
        bool ia = (d[a] >= 0.0), ib = (d[b] >= 0.0);
        if (ia) { out.push_back(a); }
        if (ia && !ib) {
          exitNode = crossing(a, b);
          out.push_back(exitNode);
        } else if (!ia && ib) {
          entryNode = crossing(b, a);
          out.push_back(entryNode);
        }
      }
      std::vector<int> clean;
      for (int k : out) {
        if (clean.empty() || clean.back() != k) { clean.push_back(k); }
      }
      while (clean.size() > 1 && clean.back() == clean.front()) { clean.pop_back(); }
      if (clean.size() >= 3) {
        pface G;
        G.nodes.assign(clean.begin(), clean.end());
        sub_block.faces.push_back(G);
      }
      if (exitNode >= 0 && entryNode >= 0 && exitNode != entryNode) {
        capSide.insert(std::make_pair(entryNode, exitNode));
      }
    }

    // The caps, in a 2D frame of the plane where their outward normal (-n) points to the
    // viewer: the outer loops are counterclockwise, the holes clockwise
    vec3r e1 = (std::fabs(n.x) < 0.9) ? cross(n, vec3r(1, 0, 0)) : cross(n, vec3r(0, 1, 0));
    e1.normalize();
    vec3r e2 = cross(-n, e1);
    std::map<int, vec2r> pos2;
    for (const auto &s : capSide) {
      for (int k : {s.first, s.second}) { pos2[k] = vec2r(sub_block.nodes[k].pos * e1, sub_block.nodes[k].pos * e2); }
    }

    // Chain the sides of the cap into loops. Where several sides leave a node (two loops
    // touch there), the one that turns the most to the right is taken: the loop follows
    // the border of the cap, which stays on its left.
    std::vector<std::vector<int>> loops;
    while (!capSide.empty()) {
      auto it   = capSide.begin();
      int start = it->first;
      std::vector<int> loop(1, start);
      int prev = start;
      int cur  = it->second;
      capSide.erase(it);
      while (cur != start) {
        auto range = capSide.equal_range(cur);
        if (range.first == range.second) { break; }
        auto jt = range.first;
        if (std::next(jt) != range.second) {
          vec2r back      = pos2[prev] - pos2[cur];
          double bestTurn = 1e300;
          for (auto kt = range.first; kt != range.second; ++kt) {
            vec2r out   = pos2[kt->second] - pos2[cur];
            double turn = std::atan2(-(back.x * out.y - back.y * out.x), back.x * out.x + back.y * out.y);
            if (turn <= 0.0) { turn += 2.0 * M_PI; } // clockwise angle from 'back' to 'out'
            if (turn < bestTurn) {
              bestTurn = turn;
              jt       = kt;
            }
          }
        }
        loop.push_back(cur);
        prev = cur;
        cur  = jt->second;
        capSide.erase(jt);
      }
      if (cur != start) {
        nbBad++;
        continue;
      }
      if (loop.size() >= 3) { loops.push_back(loop); }
    }

    std::vector<double> area(loops.size(), 0.0);
    std::vector<std::vector<vec2r>> poly2(loops.size());
    for (size_t l = 0; l < loops.size(); l++) {
      for (int k : loops[l]) { poly2[l].push_back(pos2[k]); }
      for (size_t i = 0; i < poly2[l].size(); i++) {
        const vec2r &a = poly2[l][i], &b = poly2[l][(i + 1) % poly2[l].size()];
        area[l] += 0.5 * (a.x * b.y - b.x * a.y);
      }
    }
    double eps = 1e-14 * L * L;
    std::vector<int> outerOfHole(loops.size(), -1);
    for (size_t h = 0; h < loops.size(); h++) {
      if (area[h] >= 0.0) { continue; }
      double best = 1e300;
      for (size_t o = 0; o < loops.size(); o++) {
        if (area[o] <= 0.0) { continue; }
        if (area[o] < best && inside2(poly2[h][0], poly2[o])) {
          best           = area[o];
          outerOfHole[h] = (int)o;
        }
      }
      if (outerOfHole[h] < 0 && -area[h] > eps) { nbBad++; }
    }
    for (size_t o = 0; o < loops.size(); o++) {
      if (area[o] <= eps) { continue; }
      std::vector<std::vector<int>> holes;
      for (size_t h = 0; h < loops.size(); h++) {
        if (outerOfHole[h] == (int)o) { holes.push_back(loops[h]); }
      }
      bool convex = holes.empty();
      size_t m    = loops[o].size();
      for (size_t k = 0; k < m && convex; k++) {
        if (cross2(poly2[o][(k + m - 1) % m], poly2[o][k], poly2[o][(k + 1) % m]) < -eps) { convex = false; }
      }
      if (convex) {
        pface F;
        F.nodes.assign(loops[o].begin(), loops[o].end());
        sub_block.faces.push_back(F);
      } else {
        std::vector<std::array<int, 3>> tri;
        if (!triangulate(loops[o], holes, pos2, tri, eps)) { nbBad++; }
        for (const auto &t : tri) {
          pface F;
          F.nodes = {t[0], t[1], t[2]};
          sub_block.faces.push_back(F);
        }
      }
    }

    remove_unused_nodes(sub_block);
    if (sub_block.faces.empty()) { return 2; }
    return 1;
  }

}; // end class polyhTool

} // namespace stl2blocks

namespace std {
template <> struct less<stl2blocks::polyhTool::pedge> {
  bool operator()(const stl2blocks::polyhTool::pedge &lhs, const stl2blocks::polyhTool::pedge &rhs) const {
    if (lhs.node0 < rhs.node0) { return true; }
    if ((lhs.node0 == rhs.node0) && (lhs.node1 < rhs.node1)) { return true; }
    return false;
  }
};
} // namespace std

#endif /* end of include guard: POLYHTOOL_HPP */
