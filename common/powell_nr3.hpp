#ifndef POWELL_NR3_HPP_C1B36365
#define POWELL_NR3_HPP_C1B36365

/**
@file  powell_nr3.hpp
@brief A c++ version of powell minimization in multidimensions.
       Adapted from Numerical Recipes Third Edition
Example of usage
@code{.cpp}
#include <iostream>


using namespace std;
#include "powell_nr3.hpp"

// functor
struct Func {
double operator() (vector<double> & X) {
	double xx = X[0] - 20.;
	double yy = X[1] - 100.;
	return xx*xx+yy*yy;
}
};

// classic function
double to_minimize (vector<double> & X) {
	double xx = X[0] - 20.;
	double yy = X[1] - 100.;
	return xx*xx+yy*yy;
}

int main()
{
	vector<double> X(2);
	X[0]=X[1]=0.0;
	
	// For functor
  //Func func;
	//Powell<Func> powell(func);

	// For classic functions
	Powell<double (vector<double> &)> powell(to_minimize);

	X = powell.minimize(X);
	cout << "X[0] = " << X[0] << endl;
	cout << "X[1] = " << X[1] << endl;
	cout << "fret = " << powell.fret << endl;

	return 0;
}
@endcode
*/


#include <vector>
#include <cmath>
#include <limits>
#include <iostream>


template<class T>
inline const T &MAX(const T &a, const T &b)
{
	return b > a ? (b) : (a);
}

template<class T>
inline void SWAP(T &a, T &b)
{
	T dum = a;
	a = b;
	b = dum;
}

template<class T>
inline T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

template<class T> inline T SQR(const T a)
{
	return a * a;
}


struct Bracketmethod {
	
	double ax, bx, cx, fa, fb, fc;
	template <class T>
	void bracket(const double a, const double b, T &func)
	{
		using namespace std;
		const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
		ax = a;
		bx = b;
		double fu;
		fa = func(ax);
		fb = func(bx);
		if (fb > fa) {
			SWAP(ax, bx);
			SWAP(fb, fa);
		}
		cx = bx + GOLD * (bx - ax);
		fc = func(cx);
		while (fb > fc) {
			double r = (bx - ax) * (fb - fc);
			double q = (bx - cx) * (fb - fa);
			double u = bx - ((bx - cx) * q - (bx - ax) * r) /
			           (2.0 * SIGN(MAX(abs(q - r), TINY), q - r));
			double ulim = bx + GLIMIT * (cx - bx);
			if ((bx - u) * (u - cx) > 0.0) {
				fu = func(u);
				if (fu < fc) {
					ax = bx;
					bx = u;
					fa = fb;
					fb = fu;
					return;
				}
				else if (fu > fb) {
					cx = u;
					fc = fu;
					return;
				}
				u = cx + GOLD * (cx - bx);
				fu = func(u);
			}
			else if ((cx - u) * (u - ulim) > 0.0) {
				fu = func(u);
				if (fu < fc) {
					shft3(bx, cx, u, u + GOLD * (u - cx));
					shft3(fb, fc, fu, func(u));
				}
			}
			else if ((u - ulim) * (ulim - cx) >= 0.0) {
				u = ulim;
				fu = func(u);
			}
			else {
				u = cx + GOLD * (cx - bx);
				fu = func(u);
			}
			shft3(ax, bx, cx, u);
			shft3(fa, fb, fc, fu);
		}
	}
	inline void shft2(double &a, double &b, const double c)
	{
		a = b;
		b = c;
	}
	inline void shft3(double &a, double &b, double &c, const double d)
	{
		a = b;
		b = c;
		c = d;
	}
	inline void mov3(double &a, double &b, double &c, const double d, const double e,
	                 const double f)
	{
		a = d;
		b = e;
		c = f;
	}
};


struct Brent : Bracketmethod {
	double xmin, fmin;
	const double tol;
	Brent(const double toll = 3.0e-8) : fmin(0.0), tol(toll) {}
	template <class T>
	double minimize(T &func)
	{
		using namespace std;
		const int ITMAX = 100;
		const double CGOLD = 0.3819660;
		const double ZEPS = numeric_limits<double>::epsilon() * 1.0e-3;
		double a, b, d = 0.0, etemp, fu, fv, fw, fx;
		double p, q, r, tol1, tol2, u, v, w, x, xm;
		double e = 0.0;

		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;
		fw = fv = fx = func(x);
		for (int iter = 0; iter < ITMAX; iter++) {
			xm = 0.5 * (a + b);
			tol2 = 2.0 * (tol1 = tol * abs(x) + ZEPS);
			if (abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
				fmin = fx;
				return xmin = x;
			}
			if (abs(e) > tol1) {
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if (q > 0.0) p = -p;
				q = abs(q);
				etemp = e;
				e = d;
				if (abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x)
				        || p >= q * (b - x))
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				else {
					d = p / q;
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = SIGN(tol1, xm - x);
				}
			}
			else {
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			}
			u = (abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
			fu = func(u);
			if (fu <= fx) {
				if (u >= x) a = x;
				else b = x;
				shft3(v, w, x, u);
				shft3(fv, fw, fx, fu);
			}
			else {
				if (u < x) a = u;
				else b = u;
				if (fu <= fw || w == x) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				}
				else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}
		}
		cerr << "Too many iterations in brent\n";
		return 0.0; // to shut down the compiler
	}
};



template <class T>
struct F1dim {
	const std::vector<double> &p;
	const std::vector<double> &xi;
	int n;
	T &func;
	std::vector<double> xt;
	F1dim(std::vector<double> &pp, std::vector<double> &xii, T &funcc) : p(pp),
		xi(xii), n(pp.size()), func(funcc), xt(n) {}
	double operator() (const double x)
	{
		for (int j = 0; j < n; j++)
			xt[j] = p[j] + x * xi[j];
		return func(xt);
	}
};

template <class T>
struct Linemethod {
	std::vector<double> p;
	std::vector<double> xi;
	T &func;
	int n;
	Linemethod(T &funcc) : func(funcc) {}
	double linmin()
	{
		double ax, xx, xmin;
		n = p.size();
		F1dim<T> f1dim(p, xi, func);
		ax = 0.0;
		xx = 1.0;
		Brent brent;
		brent.bracket(ax, xx, f1dim);
		xmin = brent.minimize(f1dim);
		for (int j = 0; j < n; j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return brent.fmin;
	}
};

template <class T>
struct Powell : Linemethod<T> {
	
	int iter;
	double fret;
	using Linemethod<T>::func;
	using Linemethod<T>::linmin;
	using Linemethod<T>::p;
	using Linemethod<T>::xi;
	const double ftol;
	Powell(T &func, const double ftoll = 3.0e-8) : Linemethod<T>(func),
		ftol(ftoll) {}
	std::vector<double> minimize(std::vector<double> &pp)
	{
		using namespace std;
		int n = pp.size();
		vector<vector<double> > ximat;
		ximat.resize(n);
		for (int i = 0; i < n; i++) ximat[i].resize(n);
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) ximat[i][j] = 0.0;
		for (int i = 0; i < n; i++) ximat[i][i] = 1.0;
		return minimize(pp, ximat);
	}
	// A version with possibility to initialize the diagonal of xi
	std::vector<double> minimize(std::vector<double> &pp, std::vector<double> &xivec)
	{
		using namespace std;
		int n = pp.size();
		vector<vector<double> > ximat;
		ximat.resize(n);
		for (int i = 0; i < n; i++) ximat[i].resize(n);
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) ximat[i][j] = 0.0;
		for (int i = 0; i < n; i++) ximat[i][i] = xivec[i];
		return minimize(pp, ximat);
	}
	std::vector<double> minimize(std::vector<double> &pp, std::vector<std::vector<double> > &ximat)
	{
		using namespace std;
		const int ITMAX = 200;
		const double TINY = 1.0e-25;
		double fptt;
		int n = pp.size();
		p = pp;
		vector<double> pt(n), ptt(n);
		xi.resize(n);
		fret = func(p);
		for (int j = 0; j < n; j++) pt[j] = p[j];
		for (iter = 0;; ++iter) {
			double fp = fret;
			int ibig = 0;
			double del = 0.0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) xi[j] = ximat[j][i];
				fptt = fret;
				fret = linmin();
				if (fptt - fret > del) {
					del = fptt - fret;
					ibig = i + 1;
				}
			}
			if (2.0 * (fp - fret) <= ftol * (abs(fp) + abs(fret)) + TINY) {
				return p;
			}
			if (iter == ITMAX) cerr << "powell exceeding maximum iterations.\n";
			for (int j = 0; j < n; j++) {
				ptt[j] = 2.0 * p[j] - pt[j];
				xi[j] = p[j] - pt[j];
				pt[j] = p[j];
			}
			fptt = func(ptt);
			if (fptt < fp) {
				double t = 2.0 * (fp - 2.0 * fret + fptt) * SQR(fp - fret - del) - del * SQR(fp - fptt);
				if (t < 0.0) {
					fret = linmin();
					for (int j = 0; j < n; j++) {
						ximat[j][ibig - 1] = ximat[j][n - 1];
						ximat[j][n - 1] = xi[j];
					}
				}
			}
		}
	}
};

#endif /* end of include guard: POWELL_NR3_HPP_C1B36365 */

