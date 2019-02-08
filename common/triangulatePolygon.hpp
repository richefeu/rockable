// COTD Entry submitted by John W. Ratcliff [jratcliff@verant.com]

//  THIS IS A CODE SNIPPET WHICH WILL EFFICIEINTLY TRIANGULATE ANY
//  POLYGON/CONTOUR (without holes) AS A STATIC CLASS.  THIS SNIPPET
//  IS COMPRISED OF 3 FILES, TRIANGULATE.H, THE HEADER FILE FOR THE
//  TRIANGULATE BASE CLASS, TRIANGULATE.CPP, THE IMPLEMENTATION OF
//  THE TRIANGULATE BASE CLASS, AND TEST.CPP, A SMALL TEST PROGRAM
//  DEMONSTRATING THE USAGE OF THE TRIANGULATOR.  THE TRIANGULATE
//  BASE CLASS ALSO PROVIDES TWO USEFUL HELPER METHODS, ONE WHICH
//  COMPUTES THE AREA OF A POLYGON, AND ANOTHER WHICH DOES AN EFFICENT
//  POINT IN A TRIANGLE TEST.
//  SUBMITTED BY JOHN W. RATCLIFF (jratcliff@verant.com) July 22, 2000

// =========================================================================
// THIS VERSION HAS BEEN ADAPTED FOR MBOX (vincent.richefeu@3sr-grenoble.fr)
// =========================================================================

#ifndef TRIANGULATE_POLYGON_H
#define TRIANGULATE_POLYGON_H

// Static class to triangulate any contour/polygon efficiently
// You should replace Vector2d with whatever your own Vector
// class might be.  Does not support polygons with holes.
// Uses STL vectors to represent a dynamic array of vertices.
// This code snippet was submitted to FlipCode.com by
// John W. Ratcliff (jratcliff@verant.com) on July 22, 2000
// I did not write the original code/algorithm for this
// this triangulator, in fact, I can't even remember where I
// found it in the first place.  However, I did rework it into
// the following black-box static class so you can make easy
// use of it in your own code.  Simply replace Vector2d with
// whatever your own Vector implementation might be.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>  // Include STL vector class.
#include "vec2.hpp"
#include <iostream>
//using namespace std;

class TriangulatePolygon
{
public:

	// triangulate a contour/polygon, places results in STL vector
	// as series of triangles.						
	static bool Process(const std::vector<vec2r> &contour, std::vector<int> &result) {
		// allocate and initialize list of Vertices in polygon

		int n = contour.size();
		if ( n < 3 ) return false;

		int *V = new int[n];

		// we want a counter-clockwise polygon in V

		if ( 0.0f < Area(contour) )
			for (int v = 0; v < n; v++) V[v] = v;
		else
			for (int v = 0; v < n; v++) V[v] = (n - 1) - v;

		int nv = n;

		//  remove nv-2 Vertices, creating 1 triangle every time
		int count = 2 * nv; // error detection

		for (int m = 0, v = nv - 1; nv > 2; ) {
			// if we loop, it is probably a non-simple polygon
			if (0 >= (count--)) {
				// TriangulatePolygon: ERROR - probable bad polygon!
				return false;
			}

			// three consecutive vertices in current polygon, <u,v,w>
			int u = v  ;
			if (nv <= u) u = 0; // previous
			v = u + 1;
			if (nv <= v) v = 0;     // new v
			int w = v + 1;
			if (nv <= w) w = 0; // next

			if ( Snip(contour, u, v, w, nv, V) ) {
				int a, b, c, s, t;

				// true names of the vertices
				a = V[u];
				b = V[v];
				c = V[w];

				// output Triangle
				result.push_back(a);
				result.push_back(b);
				result.push_back(c);
				m++;

				// remove v from remaining polygon
				for (s = v, t = v + 1; t < nv; s++, t++) V[s] = V[t];
				nv--;

				// resest error detection counter
				count = 2 * nv;
			}
		}

		delete V;

		return true;
	}

	// compute area of a contour/polygon
	static double Area(const std::vector<vec2r> &contour) {
		int n = contour.size();

		double A = 0.0f;

		for (int p = n - 1, q = 0; q < n; p = q++) {
			A += contour[p].x * contour[q].y - contour[q].x * contour[p].y;
		}
		return A * 0.5f;
	}
	
	// InsideTriangle decides if a point P is Inside of the triangle
	// defined by A, B, C.
	static bool InsideTriangle(double Ax, double Ay,
	                           double Bx, double By,
	                           double Cx, double Cy,
	                           double Px, double Py) {
		double ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
		double cCROSSap, bCROSScp, aCROSSbp;

		ax = Cx - Bx;
		ay = Cy - By;
		bx = Ax - Cx;
		by = Ay - Cy;
		cx = Bx - Ax;
		cy = By - Ay;
		apx = Px - Ax;
		apy = Py - Ay;
		bpx = Px - Bx;
		bpy = Py - By;
		cpx = Px - Cx;
		cpy = Py - Cy;

		aCROSSbp = ax * bpy - ay * bpx;
		cCROSSap = cx * apy - cy * apx;
		bCROSScp = bx * cpy - by * cpx;

		return ((aCROSSbp >= 0.0f) && (bCROSScp >= 0.0f) && (cCROSSap >= 0.0f));
	}
	
private:
	
	static bool Snip(const std::vector<vec2r> &contour, int u, int v, int w, int n, int *V) {
		int p;
		double Ax, Ay, Bx, By, Cx, Cy, Px, Py;
		const double EPSILON = 0.0000000001f;

		Ax = contour[V[u]].x;
		Ay = contour[V[u]].y;

		Bx = contour[V[v]].x;
		By = contour[V[v]].y;

		Cx = contour[V[w]].x;
		Cy = contour[V[w]].y;

		if ( EPSILON > (((Bx - Ax) * (Cy - Ay)) - ((By - Ay) * (Cx - Ax))) ) return false;

		for (p = 0; p < n; p++) {
			if ( (p == u) || (p == v) || (p == w) ) continue;
			Px = contour[V[p]].x;
			Py = contour[V[p]].y;
			if (InsideTriangle(Ax, Ay, Bx, By, Cx, Cy, Px, Py)) return false;
		}

		return true;
	}
};

#endif // TRIANGULATE_POLYGON_H
