#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <histo.hpp>


int main (int argc, char const *argv[])
{
	
    std::random_device rd;
    std::mt19937 gen(rd());
 
    std::normal_distribution<> d(70,10);
	// plot "pdfspline.txt" w l, 1.0/(10*sqrt(2*pi))*exp(-0.5*((x-70)/10)**2)
 
 	std::vector<double> val;
    for(int n=0; n<10000; ++n) {
        val.push_back(d(gen));
    }
 
 	std::vector<double> xs, ys;
	histo::pdfNumBinsSpline(val, xs, ys, 200);
 	
	std::ofstream dat("pdfspline.txt");
	for (int i = 0; i < xs.size(); ++i) {	
		dat << xs[i] << ' ' << ys[i] << std::endl;
	} 
	
	return 0;
}

/*
#include <cubicSpline.hpp>

int main (int argc, char const *argv[])
{
	
	std::vector<double> x(11);
	std::vector<double> y(11);
	std::ofstream dat("sindat.txt");
	for(int i = 0; i < x.size(); ++i) {
		x[i] = i;
		y[i] = atan(2*(i-10))+1.52;
		dat << x[i] << ' ' << y[i] << std::endl;
	}
	std::cout << "x.size() = " << x.size() << std::endl;



	std::vector<SplineSet> cs = spline(x, y);
	std::cout << "cs.size() = " << cs.size() << std::endl;
	for(int i = 0; i < cs.size(); ++i)
		std::cout << cs[i].d << "\t" << cs[i].c << "\t" << cs[i].b << "\t" << cs[i].a << "\t" << cs[i].x << std::endl;
	
	std::vector<double> xs, ys;
	getSlineCurve(cs, xs, ys, 10);
	std::cout << "xs.size() = " << xs.size() << std::endl;
	std::ofstream dat2("sindat2.txt");
	for(int i = 0; i < xs.size(); ++i) {	
		dat2 << xs[i] << ' ' << ys[i] << std::endl;
	} 
	
	
	std::vector<SplineSet> deriv;
	derivSpline(cs,deriv);
	xs.clear();
	ys.clear();
	getSlineCurve(deriv, xs, ys, 10);
	std::ofstream dat3("sindat3.txt");
	for(int i = 0; i < xs.size(); ++i) {	
		dat3 << xs[i] << ' ' << ys[i] << std::endl;
	} 
	
	return 0;
}
*/