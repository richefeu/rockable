#ifndef FITCDF_HPP
#define FITCDF_HPP

#include <cmath>
#include <vector>
#include <functional>
#include <iostream>
#include <map>
#include <utility>

#include "powell_nr3.hpp"
#include "KStest.hpp"

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
  //if (q<1.0 || q >= 3.0) return 0.0;
    
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

struct KS_Optim {
  std::function<double(std::vector<double> & p, double x)> func;
  std::vector<double> params;
  double *X_;
  double *Y_;
  size_t nb;
  
  double operator() (std::vector<double> & p) {
    double D = 0.0;
    for (size_t i = 0 ; i < nb ; i++) {
      double Dtest = fabs(func(p, X_[i]) - Y_[i]);
      if (Dtest > D) D = Dtest;
    }   
    return D;
  }
  
  void plugData(std::vector<double> & X, std::vector<double> & Y) {
    if (X.size() != Y.size()) {
      std::cerr << "@KS_Optim::plugData, the two vectors should have the same size\n";
      return;
    }
    X_ = &X[0];
    Y_ = &Y[0];
    nb = X.size();
    if (nb < 2) {
      std::cerr << "@KS_Optim::plugData, at least 2 values are required\n";
      return;
    }
  } 
};

template <typename T> 
void fit_cdf(std::vector<T> & data, std::vector<T> & Params, const char* cdf_func_name, T & D, T & PROB, T & pearson) {
  std::vector<double> X(data);
  std::sort(X.begin(), X.end());
  std::vector<double> Y(data.size());
  for (size_t i = 0 ; i < data.size() ; i++) {
    Y[i] = (double)i / (double)(data.size() - 1);
    //std::cout << X[i] << " " << Y[i] << '\n';
  }
  
  KS_Optim KS;
  KS.plugData(X, Y);
  std::string opt(cdf_func_name);
  if (opt == "gaussian") {
    KS.func = [](std::vector<double> & p, double x) -> double {
      if (p[1] == 0.0) return 1.0e20;
      double f = (x - p[0]) / (1.414213562373095 * p[1]);
      return 0.5 * (1.0 + erf(f));
    };
    KS.params.resize(2);
    Params.resize(2);
    KS.params[0] = 0.0;
    KS.params[1] = 1.0;
  }
  else if (opt == "q-gaussian") { // SEEMS NOT A GOOD IDEA (TOO LONG!!!!)
    KS.func = [](std::vector<double> & p, double x) -> double {
      double q = p[0];
      if (q>=3.0) p[0] = q = 2.9;
      if (q<=1.0) p[0] = q = 1.01;
      double beta = p[1];
      double Cq = C(q);
      return qgauss_cdf(beta, Cq, q, x, 50);
    };
    KS.params.resize(2);
    Params.resize(2);
    KS.params[0] = 1.0;
    KS.params[1] = 1.0;
  }
  else { // "uniform" is default
    KS.func = [](std::vector<double> & p, double x) -> double {
      return p[0] * x + p[1];
    };
    KS.params.resize(2);
    Params.resize(2);
    KS.params[0] = 1.0;
    KS.params[1] = 0.0;
  }
  
  Powell<KS_Optim> powell(KS);
  KS.params = powell.minimize(KS.params);
  for (size_t i = 0 ; i < KS.params.size() ; i++) {
    Params[i] = KS.params[i];
  }
  
  KSONE<double>(data, KS.params, KS.func, D, PROB);
  
  // pearson coefficient
  T YmeanData = 0.0;
  T YmeanFunc = 0.0;
  for (size_t i = 0 ; i < X.size() ; i++) {
    YmeanData += Y[i];
    YmeanFunc += KS.func(KS.params, X[i]);
  }
  YmeanData /= (T)X.size();
  YmeanFunc /= (T)X.size();
  T up = 0.0, downData = 0.0, downFunc = 0.0;
  for (size_t i = 0 ; i < X.size() ; i++) {
    T dData = Y[i] - YmeanData;
    T dFunc = KS.func(KS.params, X[i]) - YmeanFunc;
    up += dData * dFunc;
    downData += dData * dData;
    downFunc += dFunc * dFunc;
  }
  pearson = up / sqrt(downData * downFunc);
  
}

int main() {
  std::vector<double> dat = {
    9.3162352070e+0,  -3.2545862110e+1 ,
    4.4669927510e+1,   3.6645800380e+0 ,
    8.0808290830e+0,   2.3629311080e+0 ,
   -1.2733480010e+1,   2.0257893600e+0 ,
    5.4791947840e+0,   1.5409263380e+1 ,
   -2.3595309550e+0,  -1.0733516640e+1 ,
   -1.6899262950e+1,  -5.8189159160e+0 ,
    1.5080023840e+1,  -1.1555956730e+1 ,
    2.1044929450e+0,  -1.1076547770e+1 ,
    2.2346955510e+1,  -7.0385190690e+0 ,
   -1.8990626160e+1,   1.1074156640e+1 ,
   -2.7198891200e+1,   1.8810695100e+1 ,
    3.0499957430e+1,  -4.5409087100e+0 ,
   -1.1382574950e+1,  -9.5242045450e-1 ,
    1.1005613170e+1,   8.7685391390e-1 ,
   -1.5686441460e+1,   9.0814933260e+0 ,
   -1.0406197700e+1,  -1.7129972480e+0 ,
    2.3111045490e+1,  -2.0429355880e+1 ,
   -1.6082233770e+1,   8.3466353030e+0 ,
   -1.6294775720e+1,  -4.7211824900e+0 ,
    2.4399921570e+1,   2.4115930590e+1 ,
   -1.3390631540e+1,  -4.8785413440e+0 ,
    3.0683929750e+0,   2.3315040750e+1 ,
    1.2256443560e+1,   2.7768135560e+1 ,
    2.8230256120e+1,   8.2434614710e+0 ,
    1.1030154950e+1,   6.4625345040e+0 ,
   -1.1859315050e+1,   1.1542480980e+1 ,
   -9.6284078680e+0,   1.9036090280e+1 ,
   -6.6137205620e+0,  -8.1838489060e+0 ,
   -1.3594625300e+1,  -8.2707867490e+0 ,
    2.5436329120e+1,  -2.3590899130e+1 ,
    1.4705229290e+1,   8.7351468200e+0 ,
    1.2184040210e+1,   8.5847434510e+0 ,
   -1.2893224650e+1,  -1.6369065970e+1 ,
   -1.7599960830e+0,   2.1602172020e+1 ,
   -8.9880667070e+0,   2.2529256600e+1 ,
   -4.6563194460e+0,  -1.2762263370e+1 ,
   -2.7126736710e+1,   9.9306319010e+0 ,
   -4.6852083410e+0,   8.7273249230e+0 ,
   -1.0081846380e+1,  -9.3897553310e+0 ,
    3.7858619460e+0,   2.1386091710e+1 ,
   -1.4547575700e+1,   1.6305144730e+1 ,
    1.0759547450e+1,   3.1878885910e+1 ,
   -4.5300896260e+0,  -1.7108160430e+1 ,
    1.6384063890e+1,   2.9666981620e+0 ,
    6.0757141520e+0,  -1.6727809450e+1 ,
    6.1304043450e+0,   1.7511080530e+1 ,
   -1.7066986380e+1,  -2.0249726970e+1 ,
    9.4255720800e+0,   1.9916153580e+1 ,
   -2.3007143300e+1,  -1.4389331020e+1
  };
  
  double mu,sig;
  std::vector<double> P;
  double D, PROB, r;
  fit_cdf(dat, P, "gaussian", D, PROB, r);
  std::cout << "P[0] = " << P[0] << '\n';
  std::cout << "P[1] = " << P[1] << '\n';
  std::cout << "D = " << D << "\n";
  std::cout << "p-value = " << PROB << "\n";
  std::cout << "Pearson = " << r << '\n';
  return 0;
}

#endif /* end of include guard: FITCDF_HPP */
