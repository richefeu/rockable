// ==========================================================
// 'seer' is a viewer for Rockable, built with Dear ImGui. 
// As a playful nod to Dear ImGui, it is named Dear My Seer
// ==========================================================

#include "toofus-gate/tclap/CmdLine.h"

#include "ImGuiSDLApplication.hpp"

int main(int argc, char* argv[]) {
  ImGuiSDLApplication app(1280, 720, "Dear My Seer");

  app.box().setInteractive(true);
  app.box().showBanner();

  std::string confFileName;
  int verboseLevel = 0;

  try {
    TCLAP::CmdLine cmd("Dear My Seer", ' ', "0.1");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the conf-file", false, "conf0", "conf-file");
    TCLAP::ValueArg<int> verboseArg("v", "verbose", "Verbose level", false, 0, "int");

    cmd.add(nameArg);
    cmd.add(verboseArg);

    cmd.parse(argc, argv);

    confFileName = nameArg.getValue();
    verboseLevel = verboseArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

  app.box().setVerboseLevel(verboseLevel);
  app.box().clearMemory();
  app.box().loadConf(confFileName.c_str());

  app.box().System.read();
  app.scene3d.confNum = app.box().iconf;

  if (app.box().Particles.empty()) {
    std::cerr << "No particles! Goodbye." << std::endl;
    return 1;
  }
  
  app.box().computeAABB(app.box().nDriven);
  app.scene3d.slice.updatePos(app.box().aabb.min, app.box().aabb.max);
  app.scene3d.slice.updateN();
  app.scene3d.slice.updateWidth(app.box().aabb.min, app.box().aabb.max);

  app.box().computeAABB();
  app.scene3d.fitView();
  return app.run();
}
