//  Copyright or © or Copr. Rockable
//  
//  vincent.richefeu@3sr-grenoble.fr
//  
//  This software is a computer program whose purpose is 
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces. 
//  It is developed for an ACADEMIC USAGE
//  
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use, 
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info". 
//  
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability. 
//  
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or 
//  data to be ensured and,  more generally, to use and operate it in the 
//  same conditions as regards security. 
//  
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include <tclap/CmdLine.h>

#include "Rockable.hpp"
#include "stackTracer.hpp"

/**
 * @brief This is the command line interface (CLI) for using Rockable
 * 
 */
int main(int argc, char const* argv[]) {
  
  StackTracer::initSignals();
  
  std::string confFileName;
  int nbThreads = 1;
  int verboseLevel = 0;

  try {
    TCLAP::CmdLine cmd("This is the command line interface for Rockable", ' ', "0.3");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the conf-file", true, "conf0", "conf-file");
    TCLAP::ValueArg<int> nbThreadsArg("j", "nbThreads", "Number of threads to be used", false, 1, "int");
    TCLAP::ValueArg<int> verboseArg("v", "verbose", "Verbose level", false, 4, "int");

    cmd.add(nameArg);
    cmd.add(nbThreadsArg);
    cmd.add(verboseArg);

    cmd.parse(argc, argv);

    confFileName = nameArg.getValue();
    nbThreads = nbThreadsArg.getValue();
    verboseLevel = verboseArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

  Rockable box;
  box.setVerboseLevel(verboseLevel);
  box.showBanner();
  
#ifdef _OPENMP
  omp_set_num_threads(nbThreads);
  box.console->info("OpenMP acceleration (Number of threads = {})", nbThreads);
#else
  box.console->info("No multithreading");
#endif

  box.console_run(confFileName.c_str());
  
  return 0;
}