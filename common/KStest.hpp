#ifndef KSTEST_HPP
#define KSTEST_HPP

#include <algorithm>
#include <iostream>
#include <cmath>

// FROM: http://daviddeley.com/random/code.htm (This comes actually from numerical recipes)
// Transcripted from fortran

// The above functions use the following function for calculating the function QKS 
template <typename T>
T PROBKS(T ALAM) {
  const T EPS1 = 0.001;
  const T EPS2 = 1.0e-8;
  T A2 = -2.0 * ALAM * ALAM;
  T FAC = 2.0;
  T SUM = 0.0;
  T TERMBF = 0.0; // Previous term in sum.
  for (int J = 1; J <= 100; J++) {
    T TERM = FAC * exp(A2 * J * J);
    SUM += TERM;
    if (fabs(TERM) <= (EPS1 * TERMBF) || fabs(TERM) <= (EPS2 * SUM) ) return SUM;
    FAC = -FAC; // Alternating signs in sum.
    TERMBF = fabs(TERM);
  }
  return 1.0; // Get here only by failing to converge.
}


// Given a array of N values, DATA, and given a user-supplied function 
// of a single variable FUNC which is cumulative distribution function 
// ranging from 0 (for smallest values of argument) to 1 (for largest 
// values of its argument), this routine returns the K-S statistic D and 
// the significance level PROB. Small values of PROB show that the cumulative 
// distribution function of DATA is Significantly different from FUNC. 
// The array DATA is modified by being sorted into ascending order.
template <typename T>
void KSONE(const std::vector<T> & data, 
           std::vector<double> & Params,
           std::function<double(std::vector<double> & p, double x)> func, 
           T & D, T & PROB) {
  std::vector<T> DATA(data);
  size_t N = DATA.size();
  std::sort(DATA.begin(), DATA.end());
  
  T EN = (T)N;
  D = 0.0;
  T FO = 0.0;  // Data's c.d.f. before the next step.
  for (size_t J = 0; J < N ; J++) {     // Loop over the sorted data points.
    T FN = (T)(J + 1) / EN;     // Data's c.d.f. after this step.
    T FF = func(Params, DATA[J]);     // Compare to the user-supplied function.
    T DT = std::max(fabs(FO - FF), fabs(FN - FF));     // Maximum distance.
    if (DT > D) D = DT;
    FO = FN;
  }
  
  //PROB = PROBKS(sqrt(EN) * D); // Compute significance.
  EN = sqrt(EN);
  PROB = PROBKS((EN + 0.12 + 0.11 / EN) * D);
}

// Given an array DATA1 of N1 values, and an array DATA2 of N2 values, 
// this routine returns the K-S statistic D, and the significance level
// PROB for the null hypothesis that the data sets are drawn from the 
// same distribution. Small values of PROB show that the cumulative 
// distribution function of DATA1 is significantly different from that 
// of DATA2. The arrays DATA1 and DATA2 are modified by being sorted 
// into ascending order.
// 
template <typename T>
void KSTWO(const std::vector<T> & data1, const std::vector<T> & data2, T & D, T & PROB) {
  std::vector<T> DATA1(data1);
  std::vector<T> DATA2(data2);
  size_t N1 = DATA1.size();
  size_t N2 = DATA2.size();
  
  std::sort(DATA1.begin(), DATA1.end());
  std::sort(DATA2.begin(), DATA2.end());
  
  T EN1 = (T)N1;
  T EN2 = (T)N2;
  size_t J1 = 0; // Next value of DATA1 to be processed.
  size_t J2 = 0; // Ditto, DATA2.
  T FO1 = 0.0;   // value of c.d.f. before the next step.
  T FO2 = 0.0;   // Ditto, for DATA2.
  D = 0.0;
  while (J1 < N1 && J2 < N2) { // If we are not done...
    if (DATA1[J1] < DATA2[J2]) { // Next step is in DATA1.
      T FN1 = (T)(J1 + 1) / EN1;
      T DT = std::max(fabs(FN1 - FO2), fabs(FO1 - FO2));
      if (DT > D) D = DT;
      FO1 = FN1;
      J1++;
    }
    else { // Next step is in DATA2.
      T FN2 = (T)(J2 + 1) / EN2;
      T DT = std::max(fabs(FN2 - FO1), fabs(FO2 - FO1));
      if (DT > D) D = DT;
      FO2 = FN2;
      J2++;
    }
  }
  T EN = sqrt((EN1 * EN2) / (EN1 + EN2));
  PROB = PROBKS((EN + 0.12 + 0.11 / EN) * D);
}


/*
int main() {
  
  std::vector<double> d1 = {2.52, 3.02, 2.92, 2.84, 2.84, 2.87, 2.88, 2.18, 3.07, 2.95, 2.62, 2.95, 3.07, 3.05, 2.74, 2.74, 2.96, 3.14, 3.18, 2.92, 2.79, 2.71, 2.82, 2.89, 2.67, 3.11, 2.92, 3.10};
  std::vector<double> d2 = {2.72, 2.87, 2.90, 2.96, 3.20, 3.04, 3.20, 2.81, 3.16, 3.26, 3.33, 3.02, 3.09, 3.08, 2.89, 3.00, 2.81, 3.12, 2.93, 3.22, 2.96, 2.42, 2.93, 2.61, 3.09, 3.01, 3.10, 3.06, 3.21};
  
  double D, PROB;
  
  KSTWO<double>(d1,d2,D,PROB);
  std::cout << "D = " << D << "\n";
  std::cout << "p-value = " << PROB << "\n";
  
  return 0;
}
*/

#endif /* end of include guard: KSTEST_HPP */
