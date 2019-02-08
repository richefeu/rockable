#ifndef FITTINGTOOLS_HPP
#define FITTINGTOOLS_HPP

#include <cmath>
#include <vector>
#include <functional>
#include <iostream>
#include <map>
#include <utility>

#include "powell_nr3.hpp"

typedef std::function<double(std::vector<double> & p, double x)> paramFunc;
typedef std::map<double, const char *> fitOrder;


struct LeastSquare {
  paramFunc func;
  std::vector<double> params;
  double *X_;
  double *Y_;
  size_t nb;
  
  double operator() (std::vector<double> & p) {
    double Sum = 0.0;
    for (size_t i = 0 ; i < nb ; i++) {
      double r = func(p, X_[i]) - Y_[i];
      Sum += r * r;
    }   
    return Sum;
  }
  
  void plugData(std::vector<double> & X, std::vector<double> & Y) {
    if (X.size() != Y.size()) {
      std::cerr << "@LeastSquare::plugData, the two vectors should have the same size\n";
      return;
    }
    X_ = &X[0];
    Y_ = &Y[0];
    nb = X.size();
    if (nb < 2) {
      std::cerr << "@LeastSquare::plugData, at least 2 values are required\n";
      return;
    }
  } 
};


struct LSfit {
  LeastSquare LS;
  static std::map<const char *, std::function<paramFunc(std::vector<double> & Params)> > FuncDataBase;

  LSfit() { }
  LSfit(std::vector<double> & X, std::vector<double> & Y, 
        std::vector<double> & Params, paramFunc Func)
  {
    run(X, Y, Params, Func); 
  }
  
  LSfit(std::vector<double> & X, std::vector<double> & Y, 
        std::vector<double> & Params, const char * FuncName) 
  {
    run(X, Y, Params, FuncName);
  }
  
  static void add_distributions();
  static void add_trigo();
  
  double run(std::vector<double> & X, std::vector<double> & Y, 
        std::vector<double> & Params, paramFunc Func) {
    LS.func = Func;
    for (size_t i = 0; i < Params.size(); i++)
      LS.params.push_back(Params[i]);
    LS.plugData(X, Y);
  
    Powell<LeastSquare> powell(LS);
    LS.params = powell.minimize(LS.params);
    
    for (size_t i = 0; i < Params.size(); i++)
      Params[i] = LS.params[i];
    return powell.fret; 
  }
  
  double run(std::vector<double> & X, std::vector<double> & Y, 
        std::vector<double> & Params, const char * FuncName) 
  {
    LS.func = LSfit::FuncDataBase[FuncName](Params);
      
    for (size_t i = 0; i < Params.size(); i++)
      LS.params.push_back(Params[i]);
    LS.plugData(X, Y);
  
    Powell<LeastSquare> powell(LS);
    LS.params = powell.minimize(LS.params);
    
    for (size_t i = 0; i < Params.size(); i++)
      Params[i] = LS.params[i];
    return powell.fret;  
  }
  
  fitOrder bestFit(std::vector<double> & X, std::vector<double> & Y) {
    fitOrder order;
    
    for (auto f : LSfit::FuncDataBase) {
      //std::string fname = f.first;
      std::vector<double> Params;
      LSfit myfit;
      double ret = myfit.run(X, Y, Params, f.first);
      order.insert(std::pair<double, const char*>(ret, f.first));
    }
    
    return order;
  }
};

// ===========================================================================================
std::map<const char *, std::function<paramFunc(std::vector<double> & Params)> > LSfit::FuncDataBase = {
  { "linear", 
    [](std::vector<double> & Params) -> paramFunc {
      if (Params.size() < 2) { 
        Params.resize(2);
        Params[0] = 1.0;
        Params[1] = 0.0;
      }
      return [](std::vector<double> & p, double x) -> double {
        return p[0] + p[1] * x;
      };
    } 
  },
  
  { "poly2", 
    [](std::vector<double> & Params) -> paramFunc {
      if (Params.size() < 3) { 
        Params.resize(3);
        Params[0] = 1.0;
        Params[1] = 0.0;
        Params[2] = 0.0;
      }
      return [](std::vector<double> & p, double x) -> double {
        return p[0] + p[1] * x + p[2] * x * x;
      };
    } 
  },
  
  { "poly3", 
    [](std::vector<double> & Params) -> paramFunc {
      if (Params.size() < 4) { 
        Params.resize(4);
        Params[0] = 1.0;
        Params[1] = 0.0;
        Params[2] = 0.0;
        Params[3] = 0.0;
      }
      return [](std::vector<double> & p, double x) -> double {
        double xx = x * x;
        return p[0] + p[1] * x + p[2] * xx + p[3] * xx * x;
      };
    } 
  },
  
  { "poly4", 
    [](std::vector<double> & Params) -> paramFunc {
      if (Params.size() < 5) { 
        Params.resize(5);
        Params[0] = 1.0;
        Params[1] = 0.0;
        Params[2] = 0.0;
        Params[3] = 0.0;
        Params[4] = 0.0;
      }
      return [](std::vector<double> & p, double x) -> double {
        double xx = x * x;
        return p[0] + p[1] * x + p[2] * xx + p[3] * xx * x + p[4] * xx * xx;
      };
    } 
  },
  
  { "power", 
    [](std::vector<double> & Params) -> paramFunc {
      if (Params.size() < 2) { 
        Params.resize(2);
        Params[0] = 1.0;
        Params[1] = 1.0;
      }
      return [](std::vector<double> & p, double x) -> double {
        return p[0] * pow(x, p[1]);
      };
    } 
  },
  
  { "exponential", 
    [](std::vector<double> & Params) -> paramFunc {
      if (Params.size() < 2) { 
        Params.resize(2);
        Params[0] = 1.0;
        Params[1] = 1.0;
      }
      return [](std::vector<double> & p, double x) -> double {
        return p[0] * exp(p[1] * x);
      };
    } 
  }    
};

// ===========================================================================================
void LSfit::add_trigo() {
  LSfit::FuncDataBase.insert (
    { "cos", 
      [](std::vector<double> & Params) -> paramFunc {
        if (Params.size() < 4) { 
          Params.resize(4);
          Params[0] = 0.0;
          Params[1] = 1.0;
          Params[2] = 1.0;
          Params[3] = 0.0;
        }
        return [](std::vector<double> & p, double x) -> double {
          return p[0] + p[1] * cos(p[2] * x + p[3]);
        };
      } 
    }
  );
}
// ===========================================================================================
void LSfit::add_distributions() {
  LSfit::FuncDataBase.insert (
    { "cauchy_cdf", 
      [](std::vector<double> & Params) -> paramFunc {
        if (Params.size() < 2) { 
          Params.resize(2);
          Params[0] = 0.0; // x0, location
          Params[1] = 1.0; // gamm, scale 
        }
        return [](std::vector<double> & p, double x) -> double {
          if (p[1] == 0.0) return 1.0e20;
          double f = (x - p[0]) / p[1];
          return 0.318309886183791 * atan(f) + 0.5;
        };
      } 
    }
  );
  
  LSfit::FuncDataBase.insert (
    { "cauchy_pdf", 
      [](std::vector<double> & Params) -> paramFunc {
        if (Params.size() < 2) { 
          Params.resize(2);
          Params[0] = 0.0; // x0, location
          Params[1] = 1.0; // gamma, scale
        }
        return [](std::vector<double> & p, double x) -> double {
          if (p[1] == 0.0) return 1.0e20;
          double f = (x - p[0]) / p[1];
          return 1.0 / (3.141592653589793 * p[1] * (1.0 + f * f));
        };
      } 
    }
  );
  
  LSfit::FuncDataBase.insert (
    { "gaussian_pdf", 
      [](std::vector<double> & Params) -> paramFunc {
        if (Params.size() < 2) { 
          Params.resize(2);
          Params[0] = 0.0; // mu, mean
          Params[1] = 1.0; // sigma, stddev 
        }
        return [](std::vector<double> & p, double x) -> double {
          if (p[1] == 0.0) return 1.0e20;
          double f = (x - p[0]) / p[1];
          return exp(-0.5 * f * f) / (p[1] * 2.506628274631001);
        };
      } 
    }
  );
  
  LSfit::FuncDataBase.insert (
    { "gaussian_cdf", 
      [](std::vector<double> & Params) -> paramFunc {
        if (Params.size() < 2) { 
          Params.resize(2);
          Params[0] = 0.0; // mu, mean
          Params[1] = 1.0; // sigma, stddev 
        }
        return [](std::vector<double> & p, double x) -> double {
          if (p[1] == 0.0) return 1.0e20;
          double f = (x - p[0]) / (1.414213562373095 * p[1]);
          return 0.5 * (1.0 + erf(f));
        };
      } 
    }
  );
  
}

#if 0
int main() {
  
  std::vector<double> x = {0., 1., 2., 3., 4., 5.};
  std::vector<double> y = {0, 2., 8, 18, 32., 50}; // y = 2 * x * x
  
  //LSfit::add_distributions();
  LSfit myfit;
  fitOrder order = myfit.bestFit(x, y);
  for (auto c : order) {
    std::cout << c.second << ", r^2 = " << c.first << '\n';
  }
  
  std::cout << "Fit with " << order.begin()->second << '\n';
  std::vector<double> Params;
  double ret = myfit.run(x, y, Params, order.begin()->second);
  std::cout << "ret = " << ret << '\n';
  for (size_t i = 0; i < Params.size(); i++) std::cout << Params[i] << '\n';
  /*
  std::vector<double> Params;
  //LSfit myfit(x, y, Params, [](std::vector<double> & p, double x) -> double {return p[0] * x * x;} );
  LSfit myfit;//(x, y, Params, "poly2" );
  double ret = myfit.run(x, y, Params, "poly2");
  std::cout << "ret = " << ret << '\n';
  for (size_t i = 0; i < Params.size(); i++) std::cout << Params[i] << '\n';
  */

  return 0;
}
#endif

#endif /* end of include guard: FITTINGTOOLS_HPP */