// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#include "Rockable.hpp"
#include "PostProcessor.hpp"

int firstConf = 0;
int lastConf = 0;
int stepConf = 1;

PostProcessor *processor;

PostProcessor* readPostproCommands(const char* fname) {
  std::ifstream file(fname);
  if (!file) {
    std::cout << "Problem with reading file " << fname << std::endl;
  }
  PostProcessor *processor;
  
  kwParser parser;
  parser.kwMap["firstConf"] = __GET__(file, firstConf);
  parser.kwMap["lastConf"] = __GET__(file, lastConf);
  parser.kwMap["stepConf"] = __GET__(file, stepConf);
  parser.kwMap["PostProcessor"] = [&](std::istream & file) {
    std::string postproName;
    file >> postproName;
    std::cout << "PostProcessor: " << postproName << std::endl;
    processor = Factory<PostProcessor>::Instance()->Create(postproName);
    processor->read(file);
  };
  parser.parse(file);
  return processor;
}

int main(int argc, char const* argv[]) {
  PostProcessor *processor;
  
  // Read postprocessing commands from input file
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " <postpro_commands.txt>" << std::endl;
    return 1;
  } else {
    if (!fileTool::fileExists(argv[1])) {
      std::cout << "File " << argv[1] << " has not been found" << std::endl;
      return 1;
    }
    processor = readPostproCommands(argv[1]);
  }

  Rockable box;
  box.setInteractive(true);
  processor->plug(&box);
  
  processor->init(); // Eventually open a file (it depends on the type of PostProcessor)
  for (int n = firstConf; n <= lastConf; n += stepConf) {
    char nameConf[256];
    sprintf(nameConf, "./conf%d", n);
    std::cout << "Processed file: " << nameConf << std::endl;
    box.loadConf(nameConf);
    processor->exec();
  }
  processor->end();
  
  return 0;
}
