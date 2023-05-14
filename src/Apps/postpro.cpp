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

#include "Core/Rockable.hpp"
#include "PostProcessors/PostProcessor.hpp"

int firstConf = 0;
int lastConf = 0;
int stepConf = 1;

PostProcessor* processor;

PostProcessor* readPostproCommands(const char* fname) {
  std::ifstream file(fname);
  if (!file) {
    std::cout << "Problem with reading file " << fname << std::endl;
  }
  PostProcessor* processor;

  kwParser parser;
  parser.kwMap["firstConf"] = __GET__(file, firstConf);
  parser.kwMap["lastConf"] = __GET__(file, lastConf);
  parser.kwMap["stepConf"] = __GET__(file, stepConf);
  parser.kwMap["PostProcessor"] = [&](std::istream& file) {
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
	INIT_TIMERS();
	
  PostProcessor* processor;

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

  processor->init();  // Eventually open a file (it depends on the type of PostProcessor)
  for (int n = firstConf; n <= lastConf; n += stepConf) {
    char nameConf[256];
    snprintf(nameConf, 256, "./conf%d", n);
    std::cout << "Processed file: " << nameConf << std::endl;
    box.loadConf(nameConf);
    processor->exec();
  }
  processor->end();

  return 0;
}
