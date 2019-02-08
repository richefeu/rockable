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

#include <tclap/CmdLine.h>

#include "Rockable.hpp"

int main(int argc, char const* argv[]) {
  std::string confFileName;
  int nbThreads = 1;

  try {
    TCLAP::CmdLine cmd("This is the main runner application for Rockable", ' ', "0.3");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the conf-file", true, "conf0", "conf-file");
    TCLAP::ValueArg<int> nbThreadsArg("j", "nbThreads", "Number of threads to be used", false, 1, "int");

    cmd.add(nameArg);
    cmd.add(nbThreadsArg);

    cmd.parse(argc, argv);

    confFileName = nameArg.getValue();
    nbThreads = nbThreadsArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

#ifdef _OPENMP
  omp_set_num_threads(nbThreads);
  std::cout << "OpenMP acceleration (Number of threads = " << nbThreads << ")\n";
#else
  std::cout << "No acceleration" << std::endl;
  if (nbThreads != 1) std::cout << "nbThreads = 1\n";
#endif

  Rockable box;
  box.showBanner();

  box.loadConf(confFileName.c_str());
  box.initOutputFiles();

  box.initialChecks();

  box.System.read(true);
  box.readDataExtractors();

  if (!box.dataExtractors.empty()) {
    std::ofstream docfile("extractedDataDoc.txt");
    for (size_t d = 0; d < box.dataExtractors.size(); d++) {
      box.dataExtractors[d]->generateHelp(docfile);
      // The initialisation of dataExtractors can be necessary after conf is loaded
      // and the System is read
      box.dataExtractors[d]->init();
    }
    docfile.close();
  }

  std::cout << std::endl << std::endl;
  std::cout << msg::info() << "COMPUTATION STARTS" << msg::normal() << std::endl;
  box.UpdateNL();
  box.integrate();
  std::cout << msg::info() << "COMPUTATION NORMALLY STOPPED" << msg::normal() << std::endl;

  return 0;
}
