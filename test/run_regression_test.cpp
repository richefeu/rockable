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
#include <filesystem>
#include <regex>
#include "tclap/CmdLine.h"
//#include <MATools.hxx>
#include "Core/Rockable.hpp"
#include "stackTracer.hpp"
#include <stdexcept>

bool compareConf(std::string a_newConf, std::string a_regConf)
{
	std::ifstream newConf(a_newConf);
	std::ifstream regConf(a_regConf);

  	if (!newConf.is_open()) {
    		std::cerr << msg::warn() << " Cannot read " << a_newConf << msg::normal() << std::endl;
    		return false;
  	}

  	if (!regConf.is_open()) {
    		std::cerr << msg::warn() << " Cannot read " << a_regConf << msg::normal() << std::endl;
    		return false;
  	}

	std::string lineNewConf = std::string();
	std::string lineRegConf = std::string();

	// catch the reference line that includes starter
	std::string starter="Particles";

	// get the first line
	std::getline(regConf,lineRegConf);
	// test each line until we get a good catch
	while(lineRegConf.find(starter) == std::string::npos) 
	{
		if(!(std::getline(regConf,lineRegConf))) 
			break;
	}
	// checking step
	if(lineRegConf.find(starter) == std::string::npos) 
	{
		std::cout << " The regression file doesn't contain Particles field" << std::endl;
		return false;
	}

	// get the first line
	std::getline(newConf,lineNewConf);
	// test each line until we get a good catch
	while(lineNewConf.find(starter) == std::string::npos) 
	{
		if(!(std::getline(newConf,lineNewConf))) 
			break;
	}
	// checking step 
	if(lineNewConf.find(starter) == std::string::npos) 
	{
		std::cout << " The new file doesn't contain Particles field" << std::endl;
		return false;
	}
	

		
	while(std::getline(regConf,lineRegConf))
	{
		std::getline(newConf,lineNewConf);
		if(lineRegConf != lineNewConf)
		{

			std::cout << " error " << lineRegConf << " != " <<   lineNewConf << std::endl;
			return false;
			return false;
		}
	}
	return true;
	
}

void cleanSimulationFolder() {
  std::vector<std::string> filesToDelete = {"kineticEnergy.txt", "perf.txt", "staticBalance.txt"};
  std::regex patternToDelete("conf.*");

  std::vector<std::filesystem::path> filesToRemove;
  size_t nbConfDeleted = 0;

  for (const auto& entry : std::filesystem::directory_iterator(".")) {
    const std::string& filename = entry.path().filename().string();

    // Collect files matching the pattern "conf*"
    if (std::filesystem::is_regular_file(entry) && std::regex_match(filename, patternToDelete)) {
      filesToRemove.push_back(entry.path());
      nbConfDeleted++;
    }
  }

  // Delete files matching the pattern "conf*"
  for (const auto& fileToRemove : filesToRemove) {
    std::filesystem::remove(fileToRemove);
  }
  std::cout << "Number of conf-files deleted: " << nbConfDeleted << std::endl;

  // Delete specific files
  for (const auto& fileToDelete : filesToDelete) {
    std::filesystem::path filePath(fileToDelete);
    if (std::filesystem::exists(filePath)) {
      std::filesystem::remove(filePath);
      std::cout << "File deleted: " << fileToDelete << std::endl;
    }
  }
}

int main(int argc, char const* argv[]) {
  std::string confFileName;
  int nbThreads = 1;
  int verboseLevel = 0;
  bool cleanAndLeave = false;
  bool printBannerAndLeave = false;

  std::string newconf;
  std::string correct;

  try {
    TCLAP::CmdLine cmd("This is the main runner application for Rockable", ' ', "0.3");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the conf-file", true, "conf0", "conf-file");
    TCLAP::ValueArg<int> nbThreadsArg("j", "nbThreads", "Number of threads to be used", false, 1, "int");
	TCLAP::ValueArg<int> verboseArg(
        "v", "verbose", "Verbose level (0='off', 1='critical', 2='err', 3='warn', 4='info', 5='debug', 6='trace')",
        false, 4, "int");
    TCLAP::SwitchArg cleanArg("c", "clean", "Clean files", false);
    TCLAP::SwitchArg bannerArg("b", "banner", "show banner", false);
    TCLAP::ValueArg<std::string> _regConf ("r","regressionFile", "archive conf", true, "error", "regression-conf-file");
    TCLAP::ValueArg<std::string> _newConf ("n","newFile", "New conf file to check", true, "error", "new-conf-file");

    cmd.add(nameArg);
    cmd.add(nbThreadsArg);
	cmd.add(verboseArg);
    cmd.add(cleanArg);
    cmd.add(bannerArg);
    cmd.add(_newConf);
    cmd.add(_regConf);

    cmd.parse(argc, argv);

    confFileName = nameArg.getValue();
    nbThreads = nbThreadsArg.getValue();
	verboseLevel = verboseArg.getValue();
    cleanAndLeave = cleanArg.getValue();
    printBannerAndLeave = bannerArg.getValue();
    newconf = _newConf.getValue();
    correct = _regConf.getValue();

  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }
  
  if (cleanAndLeave) {
    cleanSimulationFolder();
    return 0;
  }

  RockableProfiler::ProfilerManager prof;
  StackTracer::initSignals();

  Rockable box;
  
  box.showBanner();
  if (printBannerAndLeave) {
    return 0;
  }
  
  box.setVerboseLevel(verboseLevel);

  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  std::cout << " Compute test: " << newconf << std::endl;
  std::cout << " Regression file: " << correct << std::endl;
  std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;

#ifdef _OPENMP
  omp_set_num_threads(nbThreads);
  box.console->info("OpenMP acceleration (Number of threads = {})", nbThreads);
#else
  box.console->info("No multithreading (compiled without OpenMP)");
#endif


box.console_run(confFileName);

  //Rockable box;
  //box.showBanner();

  //box.loadConf(confFileName.c_str());
  //box.initOutputFiles();

  //box.initialChecks();

  //const bool warning = true;
  //box.updateDrivingSystem(warning);
  //box.readDataExtractors();

  //if (!box.dataExtractors.empty()) {
  //  std::ofstream docfile("extractedDataDoc.txt");
  //  for (size_t d = 0; d < box.dataExtractors.size(); d++) {
  //    box.dataExtractors[d]->generateHelp(docfile);
  //    // The initialisation of dataExtractors can be necessary after conf is loaded
  //    // and the System is read
  //    box.dataExtractors[d]->init();
  //  }
  //  docfile.close();
  //}

  //std::cout << std::endl << std::endl;
  //std::cout << msg::info() << "COMPUTATION STARTS" << msg::normal() << std::endl;
  //box.UpdateNL();
  //box.integrate();
  //std::cout << msg::info() << "COMPUTATION NORMALLY STOPPED" << msg::normal() << std::endl;

if (argc < 3) {
  std::cout << " not enough parameters " << std::endl;
  return -1; //std::abort();
}
else {
  bool succeed = compareConf(newconf,correct);
  if(succeed){
	std::cout << newconf << " and " << correct << " are the same" << std::endl;
  }
  else{
	std::cout << "Test not passed"<<std::endl;
    return -1; //std::abort();
  }	
}
std::cout << "Coucou"<<std::endl;
return -1;
}
