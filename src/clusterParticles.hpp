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
