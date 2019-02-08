// Copyright (C) Mth <vincent.richefeu@3sr-grenoble.fr>
// 
// This file is part of mbox.
// 
// Mth can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note 
// Without a license, the code is copyrighted by default. 
// People can read the code, but they have no legal right to use it. 
// To use the code, you must contact the author directly and ask permission.

#ifndef MTH_HPP_001266EB
#define MTH_HPP_001266EB

#include <cstddef>
#include <cmath>
#include <vector>

namespace Mth
{
	
	// All const are with external linkage by default
	const double pi      = 3.14159265358979323846;
	const double piSqr   = pi * pi;
	const double pi_2    = pi / 2.0;
	const double pi_4    = pi / 4.0;
	const double _2pi    = 2.0 * pi;
	const double _1_3    = 1.0 / 3.0;
	const double _4_3    = 1.0 / 3.0;
	const double e       = 2.71828182845904523536;
	const double deg2rad = pi / 180.0;
	const double rad2deg = 180.0 / pi;

	/// @brief Gives angle between 0 and 4, 
	/// while atan2 gives an angle between -PI and PI
	template <typename T>
	T DiamondAngle(T x, T y) {
		if (y >= 0.0)
			return ( x >= 0.0 ? y / (x + y) : 1.0 - x / (-x + y) ); 
		else
			return ( x < 0.0 ? 2.0 - y / (-x - y) : 3.0 + x / (x - y) ); 
	}

	/// @brief Return the sign of a value (1 is positive, 0 is negative)
	template <typename T>
	T sign(T value) { return copysign(1.0, value); }

	/// @brief Compute the mean value and the variance of some data
	template <typename T>
	void MeanAndVariance(std::vector<T> & data, double & mean, double & var) {
		mean = 0.0;
		var = 0.0;
		if (data.size() < 2) return;
		for (size_t i = 0 ; i < data.size() ; i++) mean += data[i];
		mean /= data.size();
		double s = 0.0;
		for (size_t i = 0 ; i < data.size() ; i++) {
			s = data[i] - mean;
			var += s * s;
		}
		var /= (data.size() - 1);
	}

	/// @brief Coefficient of variation (CV) or relative standard deviation (RSD)
	template <typename T>
	T RSD(std::vector<T> & data) {
		double mean = 0.0;
		double var = 0.0;
		MeanAndVariance(data, mean, var);
		double CV = 0.0;
		if (fabs(mean) > 1e-20) CV = sqrt(var) / mean; 
		return CV;
	}

} // End of namespace Mth


#endif /* end of include guard: MTH_HPP_001266EB */
