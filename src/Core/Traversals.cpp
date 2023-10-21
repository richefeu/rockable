#include <Core/Traversals.hpp>

linkCells build_linkCells(AABB & Box, vec3r & CellMinSizes)
{
	linkCells res(Box, CellMinSizes);
	res.clear();
	res.init();
	return res;
}
