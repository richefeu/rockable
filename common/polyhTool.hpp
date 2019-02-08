// Copyright (C) polyhTool <vincent.richefeu@3sr-grenoble.fr>
// 
// This file is part of mbox.
// 
// polyhTool can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note 
// Without a license, the code is copyrighted by default. 
// People can read the code, but they have no legal right to use it. 
// To use the code, you must contact the author directly and ask permission.

#ifndef POLYHTOOL_HPP_F6351922
#define POLYHTOOL_HPP_F6351922

#include "vec2.hpp"
#include "vec3.hpp"
#include "mat9.hpp"
#include "geoTool.hpp"

#include <vector>
#include <set>
#include <map>
#include <list>
#include <utility>

class polyhTool
{
public:
	
	// =======================================================================================
	// Basic structures	
	// =======================================================================================
	struct plan
	{
		vec3r pos, normal;
	};

	// A node in a polyhedron
	struct pnode
	{
		vec3r pos;
		pnode():pos() { }
	};

	// An edge in a polyhedron
	struct pedge
	{
		int node0, node1; // index of vector<pnode> in poly
	};

	// A face in a polyhedron
	struct pface
	{
		std::list<int> nodes; // Index of vector<pnode> in poly
	};

	// polyhedron
	struct poly
	{
		std::vector<pnode> nodes;
		std::vector<pedge> edges;
		std::vector<pface> faces;
	};
	// =======================================================================================
	// =======================================================================================

	
	/// @brief Generate a sobol sequence
	/// @param[in] n number of data generated
	/// @param[out] x a vector of data of length n
	static void sobolSequence(const int n, std::vector<double> & x)
	{
		const int MAXBIT = 30, MAXDIM = 6;
		int j, k, l;
		unsigned int i, im, ipp;
		static int mdeg[MAXDIM] = {1, 2, 3, 3, 4, 4};
		static unsigned int in;
		static std::vector<unsigned int> ix(MAXDIM);
		static std::vector<unsigned int *> iu(MAXBIT);
		static unsigned int ip[MAXDIM] = {0, 1, 1, 2, 1, 4};
		static unsigned int iv[MAXDIM * MAXBIT] =
		{1, 1, 1, 1, 1, 1, 3, 1, 3, 3, 1, 1, 5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9};
		static double fac;

		if (n < 0) {
			for (k = 0; k < MAXDIM; k++) ix[k] = 0;
			in = 0;
			if (iv[0] != 1) return;
			fac = 1.0 / (1 << MAXBIT);
			for (j = 0, k = 0; j < MAXBIT; j++, k += MAXDIM) iu[j] = &iv[k];
			for (k = 0; k < MAXDIM; k++) {
				for (j = 0; j < mdeg[k]; j++) iu[j][k] <<= (MAXBIT - 1 - j);
				for (j = mdeg[k]; j < MAXBIT; j++) {
					ipp = ip[k];
					i = iu[j - mdeg[k]][k];
					i ^= (i >> mdeg[k]);
					for (l = mdeg[k] - 1; l >= 1; l--) {
						if (ipp & 1) i ^= iu[j - l][k];
						ipp >>= 1;
					}
					iu[j][k] = i;
				}
			}
		}
		else {
			im = in++;
			for (j = 0; j < MAXBIT; j++) {
				if (!(im & 1)) break;
				im >>= 1;
			}
			if (j >= MAXBIT) std::cerr << "MAXBIT too small in sobseq" << std::endl;
			im = j * MAXDIM;
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
	static int cut_poly(const poly & block, const plan & pl, poly & sub_block)
	{
		// At the beginning sub_block IS block, so we copy it
		sub_block.nodes.clear();
		sub_block.edges.clear();
		sub_block.faces.clear();

		sub_block.nodes = block.nodes;
		sub_block.edges = block.edges;
		sub_block.faces = block.faces;

		// Loop over each edge to intersect with plan
		// then add a node and split the edge in two edges
		std::map < std::pair<int,int>, int > intersected_edges; // associate a pair of nodes with a new node number
		std::pair<int,int> node_pair;
		int node0,node1;
		double alpha;
		size_t initial_node_number = sub_block.nodes.size(); // The number of nodes before we add those from the intersections
		size_t new_node_number = initial_node_number;
		std::set<int> new_nodes;
		vec3r A,B;
		pnode new_node; 
		pedge new_edge;
		int tmp;
		int nb_cut = 0;

		size_t ne = sub_block.edges.size();
		for (size_t e = 0 ; e < ne ; e++) {
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
				new_edge.node0 = new_node_number;
				new_edge.node1 = node1;
				sub_block.edges.push_back(new_edge);
		
				new_node_number++;
				nb_cut++;
			}
		}

		if (nb_cut < 3) { return 0; } // Because, we need 3 points (at least) to create the face that close the shape

		// Insert the new nodes in the faces
		std::list<int>::iterator it0, it1;
		std::map< std::pair<int,int>, int >::iterator mp;
		size_t nf = sub_block.faces.size();
		for (size_t f = 0 ; f < nf ; f++) {
			for (it0 = sub_block.faces[f].nodes.begin() ; it0 != sub_block.faces[f].nodes.end() ; ++it0) {
				it1 = it0; ++it1;
				if (it1 == sub_block.faces[f].nodes.end()) it1 = sub_block.faces[f].nodes.begin();
				node0 = *it0;
				node1 = *it1;
				if (node1 < node0) { tmp = node0; node0 = node1; node1 = tmp; }
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
		for (size_t n = 0 ; n < initial_node_number ; n++) {
			if ( (sub_block.nodes[n].pos - pl.pos) * pl.normal < 0.0 ) nodes_to_remove.insert(n);
		}


		// Remove node "to be removed" in each face
		for (size_t f = 0 ; f < sub_block.faces.size() ; f++) {
			std::list<int>::iterator it = sub_block.faces[f].nodes.begin();
			size_t nb = sub_block.faces[f].nodes.size(), count = 0;	
			while (count <= nb && it != sub_block.faces[f].nodes.end()) {
				inode = nodes_to_remove.find(*it);
				if (inode != nodes_to_remove.end()) // If the node needs to be removed 
					it = sub_block.faces[f].nodes.erase(it);
				else ++it;
				++count;
			}
		}

		// Add a face that close the sub-poly
		// together with the corresponding edges
		std::list<pedge> edge_loop0;
		int num_search = 0;
		pedge E; E.node0 = E.node1 = -1;
		for (size_t f = 0 ; f < sub_block.faces.size() ; f++) {
			num_search = 0; E.node0 = E.node1 = -1; 
			size_t count = 0;
			for (std::list<int>::iterator it = sub_block.faces[f].nodes.begin() ; it != sub_block.faces[f].nodes.end() ; ++it) {
				if ( *it >= (int)initial_node_number ) { // It means that this is a node that has been added 
					if (num_search == 0) {
						E.node0 = *it;
						num_search = 1;
					}
					else {
						E.node1 = *it;
						edge_loop0.push_back(E);
						E.node0 = E.node1;
					}
					++count;
				}
			}
			if (count == 1) std::cout << "1 seul nouveau noeud dans une face\n";
			if (count > 2)  std::cout << "plus de 2 nouveaux noeuds dans une face\n";
		}


		if (edge_loop0.size() != (size_t)nb_cut) {
			std::cout << "edge_loop0.size() = " << edge_loop0.size() << std::endl;
			std::cout << "nb_cut = " << nb_cut << std::endl;
		}

		for (std::list<pedge>::iterator ite = edge_loop0.begin() ; ite != edge_loop0.end() ; ++ite) {
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
			if (ite == edge_loop0.end()) break;
			for (ite = edge_loop0.begin() ; ite != edge_loop0.end() ; ++ite) {
				if (current_edge.node1 == ite->node1) {
					swap = ite->node0;
					ite->node0 = ite->node1;
					ite->node1 = swap;
					current_edge = *ite;
					edge_loop1.push_back(current_edge);
					edge_loop0.erase(ite);
					ite = edge_loop0.begin();
					break;
				}
				else if (current_edge.node1 == ite->node0) {
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

		pface F; F.nodes.clear();
		for (size_t e = 0 ; e <  edge_loop1.size() ; e++) {
			F.nodes.push_back(edge_loop1[e].node0);
		}
		sub_block.faces.push_back(F);

		// Remove faces with less than 3 nodes
		std::vector<pface> faces_copy;
		for (size_t f = 0 ; f < sub_block.faces.size() ; f++) {
			if (sub_block.faces[f].nodes.size() >= 3) {
				faces_copy.push_back(sub_block.faces[f]);
			} 
		}
		sub_block.faces.clear();
		sub_block.faces = faces_copy;

		// Remove edges with less than 2 nodes
		std::vector<pedge> edges_copy;
		for (size_t e = 0 ; e < sub_block.edges.size() ; e++) {
			if (nodes_to_remove.find( sub_block.edges[e].node0 ) == nodes_to_remove.end() 
				&& nodes_to_remove.find( sub_block.edges[e].node1 ) == nodes_to_remove.end())
					edges_copy.push_back(sub_block.edges[e]);
		}
		sub_block.edges.clear();
		sub_block.edges = edges_copy;

		// ====== Remove and re-number the nodes
		std::vector<int> new_number(sub_block.nodes.size());
		std::vector<pnode> nodes_copy;
		int NN = 0;
		for (size_t n = 0 ; n < sub_block.nodes.size() ; n++) {
			inode = nodes_to_remove.find(n);
			if (inode == nodes_to_remove.end()) { // ... has not been removed
				new_number[n] = NN++;
				nodes_copy.push_back(sub_block.nodes[n]);
			}	
		}
		sub_block.nodes.clear();
		sub_block.nodes = nodes_copy;

		// Re-numbering the nodes in faces
		for (size_t f = 0 ; f < sub_block.faces.size() ; f++) {
			for (std::list<int>::iterator it = sub_block.faces[f].nodes.begin() ; it != sub_block.faces[f].nodes.end() ; ++it) {
				(*it) = new_number[*it];
			}
		}

		// Re-numbering the nodes in edges
		for (size_t e = 0 ; e < sub_block.edges.size() ; e++) {
			sub_block.edges[e].node0 = new_number[ sub_block.edges[e].node0 ];
			sub_block.edges[e].node1 = new_number[ sub_block.edges[e].node1 ];
		}

		return 1;
	}

}; // end class polyhTool

namespace std {
	template<>
	struct less<polyhTool::pedge> {
		bool operator()(const polyhTool::pedge& lhs, const polyhTool::pedge& rhs) const {
			if (lhs.node0 < rhs.node0) return true;
			if ((lhs.node0 == rhs.node0) && (lhs.node1 < rhs.node1)) return true;
			return false;
		}	
	};
}


#endif /* end of include guard: POLYHTOOL_HPP_F6351922 */
