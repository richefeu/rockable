#ifndef HISTO_HPP
#define HISTO_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include "cubicSpline.hpp"

/**
@file histo.hpp
EXAMPLE:
@code{.cpp}
int main (int argc, char const *argv[])
{	
	std::vector<double> v;
	//for (size_t i = 0 ; i < 1000 ; i++) v.push_back(rand()/(double)RAND_MAX);
	for (size_t i = 0 ; i < 1000 ; i++) v.push_back(i/1000.);
	
	//histo H = histo::histoNumBins(v, 20, false);
	//histo H = histo::pdfNumBins(v, 4);
	histo H = histo::pdfMaxPerBin(v, 200);
	for (size_t i = 0 ; i < H.data.size() ; i++) std::cout << H.data[i].X << " " << H.data[i].ProbDensity << " " << H.data[i].Width << std::endl;
	
	return 0;
}
@endcode
*/

// +-----+  ^
// |     |  |
// |<-W->|  | P
// |     |  |
// +--+--+  v
//    X
struct bar {
	double X;           // mean position of the bar
	double ProbDensity; // height of the bar
	double Width;       // width of the bar
};

class histo {
public:
	double min;
	double max;
	std::vector<bar> data;

	// Ctors
	histo(size_t nbins) { data.resize(nbins); }
	histo() { }
	void clear() { data.clear(); }
	
	double entropyShannon() {
		double S = 0.0;
		for (size_t i = 0 ; i < data.size() ; i++) {
			if (data[i].ProbDensity > 0.0)
				S += data[i].ProbDensity * log10(data[i].ProbDensity);
		}
		return -S;
	}
	
	double entropyFisher() {
		double S = 0.0;
		for (size_t i = 0 ; i < data.size() ; i++) {
			// TODO
		}
		return S;
	}
	
	static histo histoNumBins(std::vector<double> & value, int nbins, bool normalize = true)
	{
		histo H(nbins);
		size_t nb = value.size();
		std::sort(value.begin(), value.end());
		double fact;
		if (normalize == true && nb > 0) fact = 1.0 / nb;
		else fact = 1.0;
		H.min = value[0];
		H.max = value[nb - 1];
		double binWidth = (H.max - H.min) / (double)nbins;
		double binAbscise;
		size_t amount;
		double threshold;
		size_t count = 0;
		for (int b = 0; b < nbins; b++) {
			binAbscise = H.min + (b + 0.5) * binWidth;
			amount = 0;
			threshold = binAbscise + 0.5 * binWidth;
			while (count < nb && value[count] <= threshold) { ++amount; ++count; }
			H.data[b].Width = binWidth;
			H.data[b].X = binAbscise;
			H.data[b].ProbDensity = amount * fact;
		}
		return H;
	}

	static histo pdfNumBins(std::vector<double> & value, int nbins)
	{
		histo H(nbins);
		std::sort(value.begin(), value.end());
		H.min = value[0];
		H.max = value[value.size() - 1];
		double binWidth = (H.max - H.min) / (double)nbins;
		double binAbscise;
		size_t amount;
		double threshold;
		size_t count = 0;
		for (int b = 0; b < nbins; b++) {
			binAbscise = H.min + (b + 0.5) * binWidth;
			amount = 0;
			threshold = binAbscise + 0.5 * binWidth;
			while (count < value.size() && value[count] <= threshold) { ++amount; ++count; }
			H.data[b].Width = binWidth;
			H.data[b].X = binAbscise;
			H.data[b].ProbDensity = amount / binWidth; // a density (not yet normalized)
		}
	
		// normalization: \int P dx = 1
		double sum = 0.0;
		for (int b = 0 ; b < nbins ; b++) sum += H.data[b].ProbDensity * H.data[b].Width;
		double invSum = 1.0;
		if (sum > 0.0) invSum = 1.0 / sum;
		for (int b = 0 ; b < nbins ; b++) H.data[b].ProbDensity *= invSum;

		return H;
	}
	
	
	
	// remarque : une meilleur solution serait de fixer nbins à une valeur pas trop grande pour limiter les "vagues"
	// et de faire une optimisation des positions des quelques points de la spline, 
  // puis on dérive une fois pour obtenir le pdf 
	static void pdfNumBinsSpline(std::vector<double> & value, std::vector<double> & xs, std::vector<double> & ys, int nbins)
	{
		std::sort(value.begin(), value.end());
		std::vector<double> x, y;
		
		size_t deltab = value.size() / nbins;
		for (size_t i = 0 ; i < value.size() ; i += nbins) {
			if (i % deltab == 0 || i == 0 || i == value.size() - 1) {
				x.push_back(value[i]);
				y.push_back((double)i / value.size());
			}
		}
		
		std::vector<SplineSet> cs = spline(x, y);
		std::vector<SplineSet> deriv;
		derivSpline(cs, deriv);
		getSlineCurve(deriv, xs, ys, 10);
	}
	
	static histo pdfNumBinsRange(std::vector<double> & value, int nbins, double min, double max)
	{
		histo H(nbins);
		std::sort(value.begin(), value.end());
		H.min = min;
		H.max = max;
		double binWidth = (H.max - H.min) / (double)nbins;
		double binAbscise;
		size_t amount;
		double threshold;
		size_t count = 0;
		for (int b = 0; b < nbins; b++) {
			binAbscise = H.min + (b + 0.5) * binWidth;
			amount = 0;
			threshold = binAbscise + 0.5 * binWidth;
			while (count < value.size() && value[count] <= threshold) { ++amount; ++count; }
			H.data[b].Width = binWidth;
			H.data[b].X = binAbscise;
			H.data[b].ProbDensity = amount / binWidth; // a density (not yet normalized)
		}
	
		// normalization: \int P dx = 1
		double sum = 0.0;
		for (int b = 0 ; b < nbins ; b++) sum += H.data[b].ProbDensity * H.data[b].Width;
		double invSum = 1.0;
		if (sum > 0.0) invSum = 1.0 / sum;
		for (int b = 0 ; b < nbins ; b++) H.data[b].ProbDensity *= invSum;

		return H;
	}

	static histo pdfMaxPerBin(std::vector<double> & value, int maxEltPerBin)
	{
		histo H;
		std::sort(value.begin(), value.end());
		H.min = value[0];
		H.max = value[value.size() - 1];
		double x0 = H.min;
		double x1 = H.min;
		size_t ivalue = 0;
		size_t prev_ivalue = 0, amount = 0;
		bar B;
		while (ivalue < value.size()) {
			ivalue = prev_ivalue + maxEltPerBin;
			if (ivalue >= value.size()) { ivalue = value.size() - 1; }
			x1 = value[ivalue];
			amount = ivalue - prev_ivalue;
			if (amount == 0) break;
			B.Width = x1 - x0;
			B.X = 0.5 * (x0 + x1);
			B.ProbDensity = (double)amount / B.Width; // a density (not yet normalized)
			H.data.push_back(B);
			x0 = x1;
			prev_ivalue = ivalue;
		}
		
		// normalization: \int P dx = 1
		double sum = 0.0;
		for (size_t b = 0 ; b < H.data.size() ; b++) sum += H.data[b].ProbDensity * H.data[b].Width;
		double invSum = 1.0;
		if (sum > 0.0) invSum = 1.0 / sum;
		for (size_t b = 0 ; b < H.data.size() ; b++) H.data[b].ProbDensity *= invSum;

		return H;
	}	
};

#endif /* end of include guard: HISTO_HPP */
