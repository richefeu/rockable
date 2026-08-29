#ifndef COLORIZER_HPP
#define COLORIZER_HPP

#include "ColorTable.hpp"

#include <algorithm>

// Base class for particle colorization strategies.
template <typename SimuClass>
class Colorizer {
 public:
  Colorizer(size_t N) { Values.resize(N); }

  virtual ~Colorizer() = default;

  virtual void init() {};
  
  virtual void calcRange(const SimuClass* box) {
    if (Values.empty()) {
      std::cout << "@Colorizer.calcBounds, colorLevels is empty" << std::endl;
    }

    // Use std::min_element and std::max_element to find min and max values
    auto minIt = std::min_element(Values.begin(), Values.end());
    auto maxIt = std::max_element(Values.begin(), Values.end());

    minVal = *minIt;
    maxVal = *maxIt;
    CT.setMinMax(minVal, maxVal);
  }
  
  void setRange(double t_minVal, double t_maxVal) {
    minVal = t_minVal;
    maxVal = t_maxVal;
    CT.setMinMax(minVal, maxVal);
  }
  
  virtual void computeValues(const SimuClass* box) const = 0;

  double minVal{0.0};
  double maxVal{0.0};
  std::vector<double> Values;
  ColorTable CT;
};

#endif  // COLORIZER_HPP
