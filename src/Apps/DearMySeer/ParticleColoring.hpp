#pragma once

#include <vector>
#include <string>
#include <GL/glew.h>

#include "ColorTable.hpp"


class Rockable;

struct ParticleColoring {
  ColorTable CT;
  std::vector<colorRGBA> pcolors;
  std::string selected{"none"};
  float xpos{0.01};
  float ypos{0.1};
  int W{20};
  int H{150};
  
  int rescale{1};
  int show_colorbar{0};
  double colorRangeMin{0.0};
  double colorRangeMax{1.0};
  
  ParticleColoring();
  void setColors(Rockable & box, const std::string & color_by);
  void draw_colorbar(Rockable & box, int width, int height);
};
