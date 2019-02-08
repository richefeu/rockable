#ifndef GEOTOOL_HPP_D1916D8F
#define GEOTOOL_HPP_D1916D8F

#include <vector>
#include <set>
#include <map>
#include <list>
#include <utility>

#include "vec2.hpp"
#include "vec3.hpp"
#include "mat9.hpp"

class geoTool
{
public:
	
	//https://www.gamedev.net/topic/338987-aabb---line-segment-intersection-test/
	static bool isIntersectedSegmentAABB(const vec3r& p1, const vec3r& p2, const vec3r& min, const vec3r& max)
	{
#define EPSILON 1e-12
		vec3r d = (p2 - p1) * 0.5f;
		vec3r e = (max - min) * 0.5f;
		vec3r c = p1 + d - (min + max) * 0.5f;
		vec3r ad = component_abs(d); // Returns same vector with all components positive

		if (fabs(c[0]) > e[0] + ad[0]) return false;
		if (fabs(c[1]) > e[1] + ad[1]) return false;
		if (fabs(c[2]) > e[2] + ad[2]) return false;

		if (fabs(d[1] * c[2] - d[2] * c[1]) > e[1] * ad[2] + e[2] * ad[1] + EPSILON) return false;
		if (fabs(d[2] * c[0] - d[0] * c[2]) > e[2] * ad[0] + e[0] * ad[2] + EPSILON) return false;
		if (fabs(d[0] * c[1] - d[1] * c[0]) > e[0] * ad[1] + e[1] * ad[0] + EPSILON) return false;
		
		return true;
#undef EPSILON
	}
	
	// Adapted from @see http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#Segment-Triangle
	static int intersectTriangle(const vec3r & orig, const vec3r & dir, const vec3r & vert0, const vec3r & vert1, const vec3r & vert2)
	{
		vec3r u, v, n;
		vec3r w0, w;
		double r, a, b;
		vec3r I;

		u = vert1 - vert0;
		v = vert2 - vert0;
		n = cross(u, v);
		// Here we suppose that the triangle is not degenerated
		w0 = orig - vert0;
		a = -(n * w0);
		b = n * dir;
		if (fabs(b) < 1.0e-15) {
			if (a == 0) return 2;
			else        return 0;
		}
		r = a / b;
		if (r < 0.0) return 0;
		I = orig + r * dir; // This is the intersection point (not returned by the function)
		double uu, uv, vv, wu, wv, D;
		uu = u * u;
		uv = u * v;
		vv = v * v;
		w = I - vert0;
		wu = w * u;
		wv = w * v;
		D = 1.0 / (uv * uv - uu * vv);
		float s, t;
		s = (uv * wv - vv * wu) * D;
		if (s < 0.0 || s > 1.0) return 0;
		t = (uv * wu - uu * wv) * D;
		if (t < 0.0 || (s + t) > 1.0) return 0;
		return 1;
	}


	// alpha = 0.0 means intersect in A
	// alpha = 1.0 means intersect in B
	// alpha in [0.0; 1.1] -> intersect segment AB
	//static void intersect_edge(const vec3r & A, const vec3r & B, const plan & pl, double & alpha)
	static void intersect_edge(const vec3r & A, const vec3r & B, const vec3r & pl_pos, const vec3r & pl_normal, double & alpha)	
	{
		vec3r PA = A - pl_pos;
		vec3r u = B - A; // do not normalize!
		double un = u * pl_normal;
		double PAn = PA * pl_normal;
		alpha = - PAn / un;
	}

	// 
	static vec3r rotatePoint(vec3r const &p, vec3r const &center, vec3r const &axis, double theta)
	{
		double const c = cos(theta), s = sin(theta);
		double const C = 1.0 - c;
		vec3r tmp = p - center;
		return center +
		       vec3r(tmp[0] * (axis[0] * axis[0] * C + c) +
		             tmp[1] * (axis[0] * axis[1] * C - axis[2] * s) +
		             tmp[2] * (axis[0] * axis[2] * C + axis[1] * s),
		             tmp[0] * (axis[1] * axis[0] * C + axis[2] * s) +
		             tmp[1] * (axis[1] * axis[1] * C + c) +
		             tmp[2] * (axis[1] * axis[2] * C - axis[0] * s),
		             tmp[0] * (axis[2] * axis[0] * C - axis[1] * s) +
		             tmp[1] * (axis[2] * axis[1] * C + axis[0] * s) +
		             tmp[2] * (axis[2] * axis[2] * C + c));
	}
	
	/// Return the volume (positive value) of a tetrahedron (four vextices)
	static double volumeTetrahedron (const vec3r &v1, const vec3r &v2, const vec3r &v3, const vec3r &v4)
	{
		double x31 = v3.x - v1.x;
		double x41 = v4.x - v1.x;
		double y31 = v3.y - v1.y;
		double y41 = v4.y - v1.y;
		double z31 = v3.z - v1.z;
		double z41 = v4.z - v1.z;
		double V =
			  (v2.x - v1.x) * (y31 * z41 - z31 * y41)
			- (v2.y - v1.y) * (x31 * z41 - z31 * x41)
			+ (v2.z - v1.z) * (x31 * y41 - y31 * x41);
		return (fabs(V) / 6.0);
	}


	/// Compute the matrix I/m of a tetrahedron.
	/// Parameters are the vertex-positions minus the point where the matrix is computed
	static mat9r inertiaDivMassTetrahedron (const vec3r &v1, const vec3r &v2, const vec3r &v3, const vec3r &v4)
	{
		double Ixx =
			 0.1 * (  v1.y * v1.y + v1.y * v2.y + v2.y * v2.y + v1.y * v3.y + v2.y * v3.y
			+ v3.y * v3.y + v1.y * v4.y + v2.y * v4.y + v3.y * v4.y + v4.y * v4.y
			+ v1.z * v1.z + v1.z * v2.z + v2.z * v2.z + v1.z * v3.z + v2.z * v3.z
			+ v3.z * v3.z + v1.z * v4.z + v2.z * v4.z + v3.z * v4.z + v4.z * v4.z );
		double Iyy =
			 0.1 * (  v1.x * v1.x + v1.x * v2.x + v2.x * v2.x + v1.x * v3.x + v2.x * v3.x
			+ v3.x * v3.x + v1.x * v4.x + v2.x * v4.x + v3.x * v4.x + v4.x * v4.x
			+ v1.z * v1.z + v1.z * v2.z + v2.z * v2.z + v1.z * v3.z + v2.z * v3.z
			+ v3.z * v3.z + v1.z * v4.z + v2.z * v4.z + v3.z * v4.z + v4.z * v4.z );
		double Izz =
			 0.1 * (  v1.x * v1.x + v1.x * v2.x + v2.x * v2.x + v1.x * v3.x + v2.x * v3.x
			+ v3.x * v3.x + v1.x * v4.x + v2.x * v4.x + v3.x * v4.x + v4.x * v4.x
			+ v1.y * v1.y + v1.y * v2.y + v2.y * v2.y + v1.y * v3.y + v2.y * v3.y
			+ v3.y * v3.y + v1.y * v4.y + v2.y * v4.y + v3.y * v4.y + v4.y * v4.y );
		double Ixy =
			 0.05 * ( 2.0 * v1.x * v1.y + v2.x * v1.y + v3.x * v1.y + v4.x * v1.y + v1.x * v2.y
			+ 2.0 * v2.x * v2.y + v3.x * v2.y + v4.x * v2.y + v1.x * v3.y + v2.x * v3.y
			+ 2.0 * v3.x * v3.y + v4.x * v3.y + v1.x * v4.y + v2.x * v4.y + v3.x * v4.y + 2.0 * v4.x * v4.y );
		double Ixz = 
			 0.05 * ( 2.0 * v1.x * v1.z + v2.x * v1.z + v3.x * v1.z + v4.x * v1.z + v1.x * v2.z
			+ 2.0 * v2.x * v2.z + v3.x * v2.z + v4.x * v2.z + v1.x * v3.z + v2.x * v3.z
			+ 2.0 * v3.x * v3.z + v4.x * v3.z + v1.x * v4.z + v2.x * v4.z + v3.x * v4.z + 2.0 * v4.x * v4.z );
		double Iyz = 0.05 *
			( 2.0 * v1.y * v1.z + v2.y * v1.z + v3.y * v1.z + v4.y * v1.z + v1.y * v2.z
			+ 2.0 * v2.y * v2.z + v3.y * v2.z + v4.y * v2.z + v1.y * v3.z + v2.y * v3.z
			+ 2.0 * v3.y * v3.z + v4.y * v3.z + v1.y * v4.z + v2.y * v4.z + v3.y * v4.z + 2.0 * v4.y * v4.z );
		return mat9r (
			Ixx, Ixy, Ixz,
			Ixy, Iyy, Iyz,
			Ixz, Iyz, Izz
		);
	}

}; // end class geoTool


#endif /* end of include guard: GEOTOOL_HPP_D1916D8F */
