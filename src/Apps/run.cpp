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

#include "run.hpp"

using ConfRow = std::vector<std::string>;

// One section of a conf-file: its header (e.g. "Interactions 192") and its rows.
struct ConfSection {
  ConfRow header;
  std::vector<ConfRow> rows;
};

/**
 *  @brief Reads a conf-file from its 'Particles' line on, cut into the sections
 *         Particles, Interactions and Interfaces. Blank and '#' lines are skipped.
 */
static bool readConfSections(const std::string& fileName, std::vector<ConfSection>& sections) {
  std::ifstream file(fileName);
  if (!file.is_open()) {
    Logger::warn("@compareConf, cannot read file {}", fileName);
    return false;
  }

  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    ConfRow row;
    std::string token;
    while (iss >> token) {
      row.push_back(token);
    }
    if (row.empty() || row[0][0] == '#') {
      continue;
    }
    if (row[0] == "Particles" || row[0] == "Interactions" || row[0] == "Interfaces") {
      sections.push_back({row, {}});
    } else if (!sections.empty()) {
      sections.back().rows.push_back(row);
    }
  }

  if (sections.empty()) {
    Logger::warn("@compareConf, no 'Particles' line in {}", fileName);
    return false;
  }
  return true;
}

static bool toNumber(const std::string& token, double& value) {
  char* end = nullptr;
  value = std::strtod(token.c_str(), &end);
  return end != token.c_str() && *end == '\0';
}

static std::string joinRow(const ConfRow& row) {
  std::string s;
  for (const auto& token : row) {
    s += (s.empty() ? "" : " ") + token;
  }
  return s;
}

/**
 *  @brief Compares a new conf-file with a reference one, from the 'Particles' line on.
 *
 *  Non-numeric tokens (shape names) and the section headers (hence the numbers of
 *  particles, interactions and interfaces) must be identical. Two numbers a (new)
 *  and b (reference) in column j of a section match when
 *
 *      |a - b| <= tolerance * max(|a|, |b|, S_j)
 *
 *  where S_j is the largest magnitude found in column j of the reference section.
 *  S_j keeps the comparison meaningful for values that are the small difference of
 *  large ones, such as the resultant force on a grain at equilibrium, whose
 *  round-off error scales with the contact forces and not with the resultant.
 *
 *  Interactions are saved in the order of their addresses in memory, which is not
 *  reproducible, so they are matched by their key (i, j, type, isub, jsub).
 */
bool compareConf(const std::string& newFileName, const std::string& refFileName, double tolerance) {
  std::vector<ConfSection> newSections, refSections;
  if (!readConfSections(newFileName, newSections) || !readConfSections(refFileName, refSections)) {
    return false;
  }

  if (newSections.size() != refSections.size()) {
    Logger::error("@compareConf, {} sections in {} but {} in {}", newSections.size(), newFileName,
                  refSections.size(), refFileName);
    return false;
  }

  auto byInteractionKey = [](const ConfRow& a, const ConfRow& b) {
    for (size_t k = 0; k < 5 && k < a.size() && k < b.size(); k++) {
      long ka = std::stol(a[k]), kb = std::stol(b[k]);
      if (ka != kb) {
        return ka < kb;
      }
    }
    return false;
  };

  const size_t maxReported = 10;
  size_t nbMismatches = 0;
  double worstDeviation = 0.0;

  for (size_t s = 0; s < refSections.size(); s++) {
    ConfSection& newSec = newSections[s];
    ConfSection& refSec = refSections[s];
    const std::string& name = refSec.header[0];

    if (newSec.header != refSec.header || newSec.rows.size() != refSec.rows.size()) {
      Logger::error("@compareConf, section '{}' ({} rows) != '{}' ({} rows)", joinRow(newSec.header),
                    newSec.rows.size(), joinRow(refSec.header), refSec.rows.size());
      return false;
    }

    if (name == "Interactions") {
      std::sort(newSec.rows.begin(), newSec.rows.end(), byInteractionKey);
      std::sort(refSec.rows.begin(), refSec.rows.end(), byInteractionKey);
    }

    std::vector<double> columnScale;
    for (const auto& row : refSec.rows) {
      columnScale.resize(std::max(columnScale.size(), row.size()), 0.0);
      for (size_t j = 0; j < row.size(); j++) {
        double value;
        if (toNumber(row[j], value)) {
          columnScale[j] = std::max(columnScale[j], std::fabs(value));
        }
      }
    }

    for (size_t r = 0; r < refSec.rows.size(); r++) {
      const ConfRow& newRow = newSec.rows[r];
      const ConfRow& refRow = refSec.rows[r];
      bool rowMatches = (newRow.size() == refRow.size());

      for (size_t j = 0; rowMatches && j < refRow.size(); j++) {
        double a, b;
        if (toNumber(newRow[j], a) && toNumber(refRow[j], b)) {
          double scale = std::max({std::fabs(a), std::fabs(b), columnScale[j]});
          double deviation = (scale > 0.0) ? std::fabs(a - b) / scale : 0.0;
          worstDeviation = std::max(worstDeviation, deviation);
          rowMatches = (deviation <= tolerance);
        } else {
          rowMatches = (newRow[j] == refRow[j]);
        }
        if (!rowMatches && nbMismatches < maxReported) {
          Logger::error("@compareConf, {} row {} column {}: {} (new) != {} (reference)", name, r, j, newRow[j],
                        refRow[j]);
        }
      }
      if (!rowMatches) {
        if (newRow.size() != refRow.size() && nbMismatches < maxReported) {
          Logger::error("@compareConf, {} row {}: '{}' (new) != '{}' (reference)", name, r, joinRow(newRow),
                        joinRow(refRow));
        }
        nbMismatches++;
      }
    }
  }

  Logger::info("@compareConf, largest relative deviation: {:.3e} (tolerance {:.1e})", worstDeviation, tolerance);
  if (nbMismatches > 0) {
    Logger::error("@compareConf, {} rows differ beyond the tolerance", nbMismatches);
    return false;
  }
  return true;
}

/**
 *  @brief Deletes files matching the pattern 'conf*', 'kineticEnergy.txt', 'perf.txt', 'staticBalance.txt',
 *         and 'checkplots.txt' in the current directory.
 */
void cleanSimulationFolder() {
  std::vector<std::string> filesToDelete = {"kineticEnergy.txt", "perf.txt", "staticBalance.txt", "checkplots.txt"};
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

/**
 * @brief This is the command line interface (CLI) for using Rockable
 *
 */
int main(int argc, char const* argv[]) {

  std::string confFileName;
  int nbThreads = 1;
  int verboseLevel = 0;
  bool cleanAndLeave = false;
  bool printBannerAndLeave = false;
  std::string newconf = "";
  std::string regconf = "";
  double tolerance = 1e-5;

  try {

    TCLAP::CmdLine cmd("This is the command line interface for Rockable", ' ', ROCKABLE_GIT_TAG);
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the conf-file", false, "conf0", "conf-file");
    TCLAP::ValueArg<int> nbThreadsArg("j", "nbThreads", "Number of threads to be used", false, 1, "int");
    TCLAP::ValueArg<int> verboseArg(
        "v", "verbose", "Verbose level (0='off', 1='critical', 2='err', 3='warn', 4='info', 5='debug', 6='trace')",
        false, 4, "int");
    TCLAP::SwitchArg cleanArg("c", "clean", "Clean files", false);
    TCLAP::SwitchArg bannerArg("b", "banner", "show banner", false);
    TCLAP::ValueArg<std::string> regConfArg("r", "regressionFile", "archive conf", false, "",
                                            "regression-conf-file");
    TCLAP::ValueArg<std::string> newConfArg("n", "newFile", "New conf file to check", false, "", "new-conf-file");
    TCLAP::ValueArg<double> toleranceArg("t", "tolerance",
                                         "Relative tolerance of the comparison between -n and -r files", false,
                                         tolerance, "double");

    cmd.add(nameArg);
    cmd.add(nbThreadsArg);
    cmd.add(verboseArg);
    cmd.add(cleanArg);
    cmd.add(bannerArg);
    cmd.add(newConfArg);
    cmd.add(regConfArg);
    cmd.add(toleranceArg);

    cmd.parse(argc, argv);

    confFileName = nameArg.getValue();
    nbThreads = nbThreadsArg.getValue();
    verboseLevel = verboseArg.getValue();
    cleanAndLeave = cleanArg.getValue();
    printBannerAndLeave = bannerArg.getValue();
    newconf = newConfArg.getValue();
    regconf = regConfArg.getValue();
    tolerance = toleranceArg.getValue();

  } catch (TCLAP::ArgException& e) {
    std::cerr << "TCLAP error: " << e.error() << " for argument " << e.argId() << std::endl;
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
  box.setOpenMPThreads(nbThreads);
  box.console_run(confFileName);

  // In case -r and -n arguments have been used
  if (!(newconf == "") && !(regconf == "")) {
    bool succeed = compareConf(newconf, regconf, tolerance);
    if (succeed) {
      Logger::info("{} and {} are the same", newconf, regconf);
    } else {
      Logger::critical("Test not passed");
      exit(-1);
    }
  }
  
  return 0;
}
