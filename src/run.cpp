//        Rockable, 3D-DEM with sphero-polyhedra
//        Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>
//        
//        This program is free software: you can redistribute it and/or modify
//        it under the terms of the GNU General Public License as published by
//        the Free Software Foundation, either version 3 of the License, or
//        (at your option) any later version.
//        
//        This program is distributed in the hope that it will be useful,
//        but WITHOUT ANY WARRANTY; without even the implied warranty of
//        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//        GNU General Public License for more details.
//        
//        You should have received a copy of the GNU General Public License
//        along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
