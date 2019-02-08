#ifndef CUBICSPLINE_HPP_4E96E93E
#define CUBICSPLINE_HPP_4E96E93E


#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>

namespace
{

struct SplineSet {
	double a;
	double b;
	double c;
	double d;
	double x;
};

std::vector<SplineSet> spline(std::vector<double> &x, std::vector<double> &y)
{
	int n = x.size() - 1;
	std::vector<double> a;
	a.insert(a.begin(), y.begin(), y.end());
	std::vector<double> b(n);
	std::vector<double> d(n);
	std::vector<double> h;

	for (int i = 0; i < n; ++i)
		h.push_back(x[i + 1] - x[i]);

	std::vector<double> alpha;
	for (int i = 0; i < n; ++i)
		alpha.push_back( 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1]  );

	std::vector<double> c(n + 1);
	std::vector<double> l(n + 1);
	std::vector<double> mu(n + 1);
	std::vector<double> z(n + 1);
	l[0] = 1;
	mu[0] = 0;
	z[0] = 0;

	for (int i = 1; i < n; ++i) {
		l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
		mu[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	l[n] = 1;
	z[n] = 0;
	c[n] = 0;

	for (int j = n - 1; j >= 0; --j) {
		c[j] = z [j] - mu[j] * c[j + 1];
		b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
		d[j] = (c[j + 1] - c[j]) / 3 / h[j];
	}

	std::vector<SplineSet> output_set(n+1);
	for (int i = 0; i < n; ++i) {
		output_set[i].a = a[i];
		output_set[i].b = b[i];
		output_set[i].c = c[i];
		output_set[i].d = d[i];
		output_set[i].x = x[i];
	}
	output_set[n].x = x[n]; // the last point is only used to save the last x
	return output_set;
}

void getSlineCurve(std::vector<SplineSet> & cs, std::vector<double> &xsv, std::vector<double> &ysv, int ndiv = 5)
{
	double dx = (cs[1].x - cs[0].x) / (double)ndiv;
	for (size_t i = 0; i < cs.size()-1; ++i) {	
		for (double xs = cs[i].x ; xs < cs[i + 1].x ; xs += dx) {
			double delta = xs - cs[i].x;
			double ys = cs[i].a + cs[i].b * delta + cs[i].c * delta * delta
				+ cs[i].d * delta * delta * delta;
			xsv.push_back(xs);
			ysv.push_back(ys);
		}
	}
	// last point
	int i2 = cs.size() - 1;
	int i1 = i2 - 1;
	double delta = cs[i2].x - cs[i1].x;
	double ys = cs[i1].a + cs[i1].b * delta + cs[i1].c * delta * delta
		+ cs[i1].d * delta * delta * delta;
	xsv.push_back(cs[i2].x);
	ysv.push_back(ys);
}

void derivSpline(std::vector<SplineSet> & cs, std::vector<SplineSet> & csd)
{
	SplineSet ss;
	
	for (size_t i = 0 ; i < cs.size() ; i++) {
		ss.a = cs[i].b;
		ss.b = 2.0 * cs[i].c;
		ss.c = 3.0 * cs[i].d;
		ss.d = 0.0;
		ss.x = cs[i].x;
		csd.push_back(ss);
	}
}


} // end of unnamed namespace

/**
@file cubicSpline.hpp
@see https://en.wikipedia.org/w/index.php?title=Spline_%28mathematics%29&oldid=288288033#Algorithm_for_computing_natural_cubic_splines
Example of usage:
@code{.cpp}
#include <iostream>
#include <cubicSpline.hpp>

int main()
{
	std::vector<double> x(11);
	std::vector<double> y(11);
	for(int i = 0; i < x.size(); ++i) {
		x[i] = i;
		y[i] = sin(i);
	}

	std::vector<SplineSet> cs = spline(x, y);
	for(int i = 0; i < cs.size(); ++i)
		std::cout << cs[i].d << "\t" << cs[i].c << "\t" << cs[i].b << "\t" << cs[i].a << std::endl;
}
@endcode
*/

#endif /* end of include guard: CUBICSPLINE_HPP_4E96E93E */
