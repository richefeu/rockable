// ==========================================================
// 'shape' is the new version of shapeSurvey for Rockable.
// It is built with Dear ImGui. 
// ==========================================================

#include "toofus-gate/tclap/CmdLine.h"

#include "ImGuiSDLApplication.hpp"

int main(int argc, char* argv[]) {
  
  ImGuiSDLApplication app(800, 800, "Dear My Shape");

  std::string shapeFileName;

  try {
    TCLAP::CmdLine cmd("Dear My Shape", ' ', "0.1");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the shape-file", false, "shapes.txt", "shape-file");

    cmd.add(nameArg);
    cmd.parse(argc, argv);
    shapeFileName = nameArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

  if (app.scene3d.readShapeLib(shapeFileName.c_str()) == 0) {
    return 0;
  }
  
  
  return app.run();
}
