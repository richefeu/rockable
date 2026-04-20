#include "ParticleColoring.hpp"
#include "Core/Rockable.hpp"

#define GLTOOLS_NO_IMPLEMENTATION
#include "glTools.hpp"

ParticleColoring::ParticleColoring() {
  CT.setTableID(MATLAB_HOT);
  CT.setSize(128);
  CT.Rebuild();
  pcolors.clear();
}

void ParticleColoring::setColors(Rockable& box, const std::string& color_by) {

  if (color_by == "velocity magnitude") {
    selected = color_by;
    if (rescale == 1) {
      colorRangeMin = 0.0;
      colorRangeMax = 0.0;
      for (size_t i = box.nDriven; i < box.Particles.size(); ++i) {
        double v = norm(box.Particles[i].vel);
        if (v > colorRangeMax) colorRangeMax = v;
      }
    }

    CT.setTableID(MATLAB_HOT);
    CT.setSize(128);
    CT.setMinMax(colorRangeMin, colorRangeMax);
    CT.Rebuild();
    pcolors.clear();
    colorRGBA col;
    for (size_t i = 0; i < box.Particles.size(); ++i) {
      double v = norm(box.Particles[i].vel);
      CT.getRGB(v, &col);
      pcolors.push_back(col);
    }
  }

  else if (color_by == "particle shape name") {
    selected = color_by;
    colorRangeMin = 0.0;
    colorRangeMax = box.Shapes.size() - 1;

    CT.setTableID(RANDOM);
    CT.setSize(box.Shapes.size());
    CT.setMinMax(colorRangeMin, colorRangeMax);
    CT.Rebuild(23450);
    pcolors.clear();
    colorRGBA col;
    for (size_t i = 0; i < box.Particles.size(); ++i) {
      std::string pname = box.Particles[i].shape->name;
      size_t id = box.shapeId[pname];

      double v = static_cast<double>(id);
      CT.getRGB(v, &col);
      pcolors.push_back(col);
    }
  } else {
    pcolors.clear();
    show_colorbar = 0;
    selected = "none";
  }
}

void ParticleColoring::draw_colorbar(Rockable & box, int width, int height) {

  if (selected == "velocity magnitude") {
    glColorBar CB;
    CB.setPos(xpos*width, ypos*height);
    CB.setSize(W, H);
    CB.setTitle("Velocity");
    CB.show(width, height, CT);
  } else if (selected == "particle shape name") {
    glColorBar CB;
    CB.setPos(xpos*width, ypos*height);
    CB.setSize(W, H);
    CB.setTitle("Shape names");

    size_t step = 1;
    if (box.Shapes.size() > 8) {
      step = box.Shapes.size() / 8;
      int hh = 15 * box.Shapes.size();
      if (hh > 450) {
        hh = 450;
      }
      CB.setSize(20, hh);
    }

    for (size_t i = 0; i < box.Shapes.size(); i += step) {
      CB.addLabel(i, box.Shapes[i].name, CT);
    }
    CB.show(width, height, CT);
  }
}
