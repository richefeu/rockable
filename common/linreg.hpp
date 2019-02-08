#ifndef LINREG_HPP
#define LINREG_HPP

#include <vector>
#include <iostream>
#include <cmath>

namespace {
	
	/// @brief Make a linear regression of a set of (x, y) pairs
	/// @return Returns a, b so that the line ax+b is the linear regression of data x and y.
	///         r is the coefficient of correlation.
	void linreg(std::vector<double> &x, std::vector<double> &y, double & a, double & b, double &r)
	{
		if (x.size() != y.size()) {
			std::cerr << "@linreg, x and y must have the same size" << std::endl;
			return;
		}
		int N = static_cast<int>(x.size());
		if (N == 0) return;
		double invN = 1.0f / (double)N;

		double xmean = 0.0;
		double ymean = 0.0;
		for (int i = 0 ; i < N ; i++) {
			xmean += x[i];
			ymean += y[i];
		}
		xmean *= invN;
		ymean *= invN;

		double Sx = 0.0;
		double Sy = 0.0;
		double Sxy = 0.0;
		for (int i = 0 ; i < N ; i++) {
			double dx = x[i] - xmean;
			double dy = y[i] - ymean;
			Sx += dx * dx;
			Sy += dy * dy;
			Sxy += dx * dy;
		}

		a = Sxy / Sx; // suppose that Sx cannot be zero
		b = ymean - a * xmean;
		if (Sy == 0.0) {
			r = 1.0;
		}
		else {
			r = Sxy / sqrt(Sx * Sy);
		}
	}

} // end namespace

#endif /* end of include guard: LINREG_HPP */
