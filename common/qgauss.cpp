#include <iostream>
#include <cmath>
using namespace std;

// Tsallis q-exponetial
double eq(double q, double x) {
  if (q == 1.0) return exp(x);
  double minus = 1.0 - q;
  double n = 1.0 + minus * x;
  if (n <= 0.0) return 0.0; 
  return pow(n, 1.0 / minus);
}

// Normalisation coefficient for q-gaussian function
double C(double q) {
  const double sqrtPi = sqrt(M_PI);
  if (q < 1.0) {
    double onemq = 1.0 - q;
    double threemq = 3.0 - q;
    return 2.0*sqrtPi*tgamma(1/onemq)/(threemq*sqrt(onemq)*tgamma(threemq/(2.0*onemq)));
  }
  else if (q == 1.0) {
    return sqrtPi;
  }
  else if (q < 3.0) {
    double qm1 = q - 1.0;
    return sqrtPi*tgamma((3.0-q)/(2.0*qm1))/(sqrt(qm1)*tgamma(1.0/qm1));
  }
  return 0.0; 
}

double qgauss(double beta, double Cq, double q, double x) {
  return (sqrt(beta)/Cq)*eq(q,-beta*x*x);
}

// Function to evalute the value of integral from -infty to x
// n is the number of divisions between -10 and x
double qgauss_cdf(double beta, double Cq, double q, double x, double n) 
{ 
  double a = -5.0;
  while (qgauss(beta, Cq, q, a) > 1e-6) {
    a *= 2;
    n *= 2;
  }
  // Grid spacing 
  double h = (x - a) / (double)n;  
  // Computing sum of first and last terms 
  // in above formula 
  double s = qgauss(beta, Cq, q, a) + qgauss(beta, Cq, q, x); 
  // Adding middle terms in above formula 
  for (int i = 1; i < n; i++) 
    s += 2 * qgauss(beta, Cq, q, a + i * h); 
  // h/2 indicates (b-a)/2n. Multiplying h/2 with s. 
  return (h / 2) * s; 
} 


int main ()
{ 
  double q = 2.0;
  double beta = 3.0;
  double Cq = C(q);
  for (double x = -4.0; x < 4.0; x += 0.1) {
    cout << x << ' ' << qgauss(beta, Cq, q, x) << ' ' << qgauss_cdf(beta, Cq, q, x, 100) << endl;
  }
  
  
  return 0;
}
