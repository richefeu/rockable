#ifndef LINKCELLS_HPP_DA634EBB
#define LINKCELLS_HPP_DA634EBB

#include <vector>

#include "vec3.hpp"
#include "AABB.hpp"

// An axis aligned box used for neighbor tracking
class AABB_Cell {
public:
	std::vector<size_t>   bodies; ///< Holded bodies
	std::vector<AABB_Cell*>  pcells; ///< Surroundind cells (+ current cell)

	AABB_Cell()  { }
	~AABB_Cell() { }
};


class linkCells {
private:
  vec3r  factor;    ///< Used to convert position to index
  
public:
	AABB   box; ///< Overall surrounding box
	vec3r  minSizes; ///< Wanted minimum size (along x, y and z) of the cells
	std::vector<std::vector<std::vector<AABB_Cell> > > cells; ///< Cells that hold only free bodies
	AABB_Cell oversized_bodies; ///< A particular cell that hold only 'too big' bodies
	vec3ui N; ///< Number of cells in each direction (x, y, and z)

	linkCells(AABB & Box, vec3r & CellMinSizes): box(Box), minSizes(CellMinSizes) { init(); }
	~linkCells() { }

	void init() {
		// Partition
		N.x = (size_t)floor((box.max.x - box.min.x) / minSizes.x);
		N.y = (size_t)floor((box.max.y - box.min.y) / minSizes.y);
		N.z = (size_t)floor((box.max.z - box.min.z) / minSizes.z);
		if (N.x < 1) N.x = 1;
		if (N.y < 1) N.y = 1;
		if (N.z < 1) N.z = 1;

		real dx, dy, dz;
		dx = (box.max.x - box.min.x) / (real)N.x;
		dy = (box.max.y - box.min.y) / (real)N.y;
		dz = (box.max.z - box.min.z) / (real)N.z;

		factor.set(
			(dx > 0.0) ? 1.0 / dx : 0.0,
			(dy > 0.0) ? 1.0 / dy : 0.0,
			(dz > 0.0) ? 1.0 / dz : 0.0
		);

		// Reserve memory
		AABB_Cell A;
		std::vector<AABB_Cell> kvec;
		for (size_t k = 0 ; k < N.z ; ++k) { kvec.push_back(A); }
		std::vector<std::vector<AABB_Cell> > jvec;
		for (size_t j = 0 ; j < N.y ; ++j) jvec.push_back(kvec);
		for (size_t i = 0 ; i < N.x ; ++i) cells.push_back(jvec);

		// Link the cells
		size_t ix0, ix1;
		size_t iy0, iy1;
		size_t iz0, iz1;
		for (size_t ix = 0 ; ix < N.x ; ++ix) {
			for (size_t iy = 0 ; iy < N.y ; ++iy) {
				for (size_t iz = 0 ; iz < N.z ; ++iz) {

					ix0 = (ix > 0)       ? ix - 1 : ix;
					ix1 = (ix < N.x - 1) ? ix + 1 : ix;
					iy0 = (iy > 0)       ? iy - 1 : iy;
					iy1 = (iy < N.y - 1) ? iy + 1 : iy;
					iz0 = (iz > 0)       ? iz - 1 : iz;
					iz1 = (iz < N.z - 1) ? iz + 1 : iz;

					for (size_t iix = ix0 ; iix <= ix1 ; ++iix) {
						for (size_t iiy = iy0 ; iiy <= iy1 ; ++iiy) {
							for (size_t iiz = iz0 ; iiz <= iz1 ; ++iiz) {
								cells[ix][iy][iz].pcells.push_back(&cells[iix][iiy][iiz]);
							}
						}
					}

					cells[ix][iy][iz].pcells.push_back(&oversized_bodies); // any cell 'see' the oversized_bodies
				}
			}
		}
	}

	void clear() {
		for (size_t ix = 0 ; ix < N.x ; ++ix) {
			for (size_t iy = 0 ; iy < N.y ; ++iy) {
				for (size_t iz = 0 ; iz < N.z ; ++iz) {
					cells[ix][iy][iz].bodies.clear();
				}
			}
		}
		oversized_bodies.bodies.clear();
	}

	void add_body(size_t B, vec3r & pos, vec3r diag) {
		int ix, iy, iz;
	    
		if (diag.x > minSizes.x || diag.y > minSizes.y || diag.z > minSizes.z) oversized_bodies.bodies.push_back(B);
		else {
			ix = (int)trunc((pos.x - box.min.x) * factor.x);
			iy = (int)trunc((pos.y - box.min.y) * factor.y);
			iz = (int)trunc((pos.z - box.min.z) * factor.z);

			if (ix < 0 || ix >= (int)N.x) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
			if (iy < 0 || iy >= (int)N.y) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
			if (iz < 0 || iz >= (int)N.z) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
			cells[ix][iy][iz].bodies.push_back(B);
		}
	}

	void add_body(size_t B, vec3r & pos, AABB & aabb) {
		int ix, iy, iz;
	  vec3r diag = aabb.max - aabb.min;
    
		if (diag.x > minSizes.x || diag.y > minSizes.y || diag.z > minSizes.z) oversized_bodies.bodies.push_back(B);
		else {
			ix = (int)trunc((pos.x - box.min.x) * factor.x);
			iy = (int)trunc((pos.y - box.min.y) * factor.y);
			iz = (int)trunc((pos.z - box.min.z) * factor.z);

			if (ix < 0 || ix >= (int)N.x) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
			if (iy < 0 || iy >= (int)N.y) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
			if (iz < 0 || iz >= (int)N.z) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
			cells[ix][iy][iz].bodies.push_back(B);
		}
	}
  
};

#endif /* end of include guard: LINKCELLS_HPP_DA634EBB */
