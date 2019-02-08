#ifndef SIGNALSTAT_HPP_7C22741C
#define SIGNALSTAT_HPP_7C22741C

#include <vector>
#include <utility>
#include <cmath>

class signalStat
{
public:

	/// Visibility graph according to Lacasa et al (PNAS)
	/// The basic algorithm
	static void visibilityGraphBasic(std::vector<double> & t, std::vector<double> & y, std::vector<std::pair<int,int>> & graph) {
		graph.clear();
		int itmax = t.size() - 1;
		for (int ia = 0 ; ia < itmax ; ia++) {
			graph.push_back(std::pair<int,int>(ia, ia + 1));
			for (int ib = ia + 2 ; ib <= itmax ; ib++) {
				bool connect = true;
				double fact = (y[ia] - y[ib]) / (t[ib] - t[ia]);
				for (int ic = ia + 1 ; ic < ib ; ic++) {
					double ycrit = y[ib] + fact * (t[ib] - t[ic]);
					if (y[ic] > ycrit) {
						connect = false;
						break;
					}
				}	
				if (connect == true) {
					graph.push_back(std::pair<int,int>(ia, ib));
				}
			}
		}
	}
	
	/// Visibility graph according to Lacasa et al (PNAS)
	/// A little bit faster algorithm
	static void visibilityGraphFast(std::vector<double> & t, std::vector<double> & y, std::vector<std::pair<int,int>> & graph) {
		graph.clear();
		int itmax = t.size() - 1;
		for (int ia = 0 ; ia < itmax ; ia++) {
			graph.push_back(std::pair<int,int>(ia, ia + 1));
			int im = ia + 1;
			for (int ib = ia + 2 ; ib <= itmax ; ib++) {
				bool connect = true;
				double fact = (y[ia] - y[ib]) / (t[ib] - t[ia]);
			
				double ycrit = y[ib] + fact * (t[ib] - t[im]);
				if (y[im] > ycrit) {
					connect = false;
				}
				else for (int ic = im ; ic < ib ; ic++) {		
					if (y[ic] > y[im]) im = ic;
					double ycrit = y[ib] + fact * (t[ib] - t[ic]);
					if (y[ic] > ycrit) {
						connect = false;
						break;
					}
				}	
				if (connect == true) {
					graph.push_back(std::pair<int,int>(ia, ib));
				}
			}
		}
	}
	
	// Adapted from the function weir of Allbens Attman
	// g is gamma, h is the Hurst, x0 and xf is the start and end of t, and nt is the number of points
	static void WeierstrassFunction(double g, double h, double x0, double xf, int nt, std::vector<double> & t, std::vector<double> & y) {

		int n, k, N = 10;
		double gn[115];
		double x, dx, c, xn, a, b, d1, d2;

		dx = (xf - x0) / nt;
		for (n = 0 ; n <= 2 * N ; ++n){
			xn = 1.0 * (n - N);
			gn[n] = pow(g, xn);
		}
		x = x0;
		while (x < xf) {
			c = 0.0;
			k = 0;
			for (n=0;n<=2*N;++n) c += (1.0 - cos(gn[n] * x)) / pow(gn[n], h);
			n = N;
			while(k == 0 && n < 100) {    
				++n;  
				a = pow(g, 1.0 * n);
				b = (1.0 - cos(a * x)) / pow(a, h);
				d1 = fabs(b / c);
				c += b;
				a = 1.0 / a;
				b = (1.0 - cos(a * x)) / pow(a, h);
				d2 = fabs(b / c);
				if (d1 < 1.0e-06 && d2 < 1.0e-06) k = 1;
				c += b;
			}

			//fprintf(stdout, "%12.9le %12.9le %5d %12.9le %12.9le \n",x,c,n,d1,d2);
			t.push_back(x);
			y.push_back(c);
			x = x + dx;
		}	
	}
};

#endif /* end of include guard: SIGNALSTAT_HPP_7C22741C */
