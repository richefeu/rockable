#ifndef POLYREG_HPP
#define POLYREG_HPP

#include <vector>
#include<iostream>
#include<iomanip>
#include<cmath>

namespace {

/// @fn      polyreg
/// @brief   Polynomial Fit
/// @return  'a' as a vector of n coefficients
///          n is the degree of Polynomial 
void polyreg(std::vector<double> &x, std::vector<double> &y, int n, std::vector<double> & a)
{
	if (x.size() != y.size()) {
		std::cerr << "@polyreg, vectors x and y are not the same size." << std::endl;
	}
	int N = (int)(x.size());
    a.resize(n + 1, 0.0); // values of the final coefficients
	
    double X[2 * n + 1]; // Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    for (int i = 0 ; i < 2 * n + 1 ; i++) {
        X[i] = 0.0;
        for (int j = 0 ; j < N ; j++) {
			// Consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
        	X[i] = X[i] + pow(x[j], i); 
        }
    }
    
	double B[n + 1][n + 2]; // B is the Normal matrix(augmented) that will store the equations
    for (int i = 0 ; i <= n ; i++) {
        for (int j = 0 ; j <= n ; j++) {
        	B[i][j] = X[i + j]; 
			// Build the Normal matrix by storing the corresponding coefficients 
			// at the right positions except the last column of the matrix
        }
    }
                    
    double Y[n+1]; // Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    for (int i = 0 ; i < n + 1 ; i++) {    
        Y[i] = 0.0;
        for (int j = 0 ; j < N ; j++) {
			// Consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
        	Y[i] = Y[i] + pow(x[j], i) * y[j]; 
        }
    }
	
	// Load the values of Y as the last column of B(Normal Matrix but augmented)
    for (int i = 0 ; i <= n ; i++) B[i][n + 1] = Y[i]; 
	
	// n is made n+1 because the Gaussian Elimination part below was for n equations, 
	// but here n is the degree of polynomial and for n degree we get n+1 equations
    n = n + 1; 
      
    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
	for (int i = 0 ; i < n ; i++) {
		for (int k = i + 1 ; k < n ; k++) {
			if (B[i][i] < B[k][i]) {
				for (int j = 0 ; j <= n ; j++) {
					std::swap(B[i][j],B[k][j]);
				}
			}
		}
	}                   
                
    //loop to perform the gauss elimination
	for (int i = 0 ; i < n - 1 ; i++) {           
		for (int k = i + 1 ; k < n ; k++) {
			double t = B[k][i] / B[i][i];
			// make the elements below the pivot elements equal to zero or elimnate the variables
			for (int j = 0 ; j <= n ; j++) B[k][j] = B[k][j] - t * B[i][j]; 
		}
	}
	
	// Back-substitution
    for (int i = n - 1 ; i >= 0 ; i--) { // x is an array whose values correspond to the values of x,y,z..
        a[i] = B[i][n]; // make the variable to be calculated equal to the rhs of the last equation
		// Then subtract all the lhs values except the coefficient of the variable whose value is being calculated
        for (int j = 0 ; j < n ; j++) if (j!=i) a[i] = a[i] - B[i][j] * a[j];
		//now finally divide the rhs by the coefficient of the variable to be calculated
        a[i] = a[i] / B[i][i];            
    }
}

} // end of unnamed namespace

#endif /* end of include guard: POLYREG_HPP */





