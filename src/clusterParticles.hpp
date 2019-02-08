//        Rockable, 3D-DEM with sphero-polyhedra
//        Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>
//        
//        This program is free software: you can redistribute it and/or modify
//        it under the terms of the GNU General Public License as published by
//        the Free Software Foundation, either version 3 of the License, or
//        (at your option) any later version.
//        
//        This program is distributed in the hope that it will be useful,
//        but WITHOUT ANY WARRANTY; without even the implied warranty of
//        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//        GNU General Public License for more details.
//        
//        You should have received a copy of the GNU General Public License
//        along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef CLUSTERPARTICLES_HPP_2A8381E6
#define CLUSTERPARTICLES_HPP_2A8381E6

/// @brief A container struct with the cluster identifier 
///        and the identifiers of the held particles 
struct clusterParticles {
	size_t clusterId;
	std::vector<size_t> particleId; ///< A particleId is the index of the particle in the vector Rockable::Particles 
};

namespace std {
	template <>
	struct less<clusterParticles> {
		bool operator() (const clusterParticles & lhs, const clusterParticles & rhs) const
		{
			if (lhs.clusterId < rhs.clusterId) return true;
			return false;
		}
	};
}

#endif /* end of include guard: CLUSTERPARTICLES_HPP_2A8381E6 */
