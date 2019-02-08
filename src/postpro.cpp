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
