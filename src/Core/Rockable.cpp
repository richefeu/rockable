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

#define CONF_VERSION_DATE "21-08-2022"
#include "Rockable.hpp"

#include "BodyForces/BodyForce_AttractingPoint.hpp"
#include "BodyForces/BodyForce_PreferredDirection.hpp"
#include "BodyForces/BodyForce_ViscousFluid.hpp"

#include "DataExtractors/DataExtractor_ClusterAABB.hpp"
#include "DataExtractors/DataExtractor_DuoBalance.hpp"
#include "DataExtractors/DataExtractor_MeanVelocity.hpp"
#include "DataExtractors/DataExtractor_TrackBody.hpp"
#include "DataExtractors/DataExtractor_TrackRockfall.hpp"
#include "DataExtractors/DataExtractor_dnStat.hpp"

#include "ForceLaws/ForceLaw_Avalanche.hpp"
#include "ForceLaws/ForceLaw_Default.hpp"
#include "ForceLaws/ForceLaw_StickedLinks.hpp"

#include "PreproCommands/PreproCommand_copyParamsToInterfaces.hpp"
#include "PreproCommands/PreproCommand_homothetyRange.hpp"
#include "PreproCommands/PreproCommand_particlesClonage.hpp"
#include "PreproCommands/PreproCommand_randomlyOrientedVelocities.hpp"
#include "PreproCommands/PreproCommand_randomlyOrientedVelocitiesClusters.hpp"
#include "PreproCommands/PreproCommand_setAllVelocities.hpp"
#include "PreproCommands/PreproCommand_setStiffnessRatioInterfaces.hpp"
#include "PreproCommands/PreproCommand_setVariableStickParams.hpp"
#include "PreproCommands/PreproCommand_stickClusters.hpp"
#include "PreproCommands/PreproCommand_stickVerticesInClusters.hpp"
#include "PreproCommands/PreproCommand_stickVerticesInClustersMoments.hpp"

// #include UNSHARED_FOLDER "/includes.hpp"

// ==============================================================================================================
//  INITIALISATIONS
// ==============================================================================================================

/**
 *   Construct a new Rockable:: Rockable object
 *
 */
Rockable::Rockable() {
  // Some default values (actually, most of them will be reset after)
  t = 0.0;
  tmax = 1.0;
  computationStopAsked = 0;
  dt = 1e-6;
  interVerletC = 0.0;
  interVerlet = 0.01;
  interConfC = 0.0;
  interConf = 0.25;
  iconf = 0;
  DVerlet = 0.0;
  dVerlet = 0.0;
  nDriven = 0;
  shapeFile = "noShapeFile";

  gravity.set(0.0, -9.81, 0.0);
  bodyForce = nullptr;

#ifdef ROCKABLE_ENABLE_PERIODIC
  usePeriodicCell = 0;
#endif

  useSoftParticles = 0;
  Cinv.set(1e10, 0.0);

  dynamicUpdateNL = 0;
  dispUpdateNL = 1.0;
  angleUpdateNL = 1.0;

  numericalDampingCoeff = 0.0;
  velocityBarrier = 0.0;
  angularVelocityBarrier = 0.0;
  velocityBarrierExponent = 1.0;
  angularVelocityBarrierExponent = 1.0;

  paramsInInterfaces = 0;
  idDensity = properties.add("density");

  optionNames["forceLaw"] = "Default";

  AddOrRemoveInteractions = [this](size_t i, size_t j, double dmax) -> int {
    return this->AddOrRemoveInteractions_bruteForce(i, j, dmax);
  };
  optionNames["AddOrRemoveInteractions"] = "bruteForce";
  UpdateNL = [this]() { this->UpdateNL_bruteForce(); };

  optionNames["UpdateNL"] = "bruteForce";
  cellMinSizes.set(1.0, 1.0, 1.0);
  boxForLinkCellsOpt = 0;

  integrationStep = [this]() { this->velocityVerletStep(); };
  optionNames["Integrator"] = "velocityVerlet";

  idKnContact = dataTable.add("knContact");
  idEn2Contact = dataTable.add("en2Contact");
  idKtContact = dataTable.add("ktContact");
  idMuContact = dataTable.add("muContact");
  idKrContact = dataTable.add("krContact");
  idMurContact = dataTable.add("murContact");

  idKnInnerBond = dataTable.add("knInnerBond");
  idKtInnerBond = dataTable.add("ktInnerBond");
  idKrInnerBond = dataTable.add("krInnerBond");
  idEn2InnerBond = dataTable.add("en2InnerBond");
  idFn0InnerBond = dataTable.add("fn0InnerBond");
  idFt0InnerBond = dataTable.add("ft0InnerBond");
  idMom0InnerBond = dataTable.add("mom0InnerBond");
  idPowInnerBond = dataTable.add("powInnerBond");

  idKnOuterBond = dataTable.add("knOuterBond");
  idKtOuterBond = dataTable.add("ktOuterBond");
  idKrOuterBond = dataTable.add("krOuterBond");
  idEn2OuterBond = dataTable.add("en2OuterBond");
  idFn0OuterBond = dataTable.add("fn0OuterBond");
  idFt0OuterBond = dataTable.add("ft0OuterBond");
  idMom0OuterBond = dataTable.add("mom0OuterBond");
  idPowOuterBond = dataTable.add("powOuterBond");

  needUpdate = false;
  interactiveMode = false;
  glue_with_walls = false;

  preventCrossingLength = 0.0;

  ExplicitRegistrations();
  initParser();
}

/**
 * Registers various features explicitly to build a static library.
 *
 * With automatic registration, merging the *.o into a single *.a breaks everything.
 *
 * Note: Compilation is slightly slower due to explicit registrations.
 */
void Rockable::ExplicitRegistrations() {
  // BodyForces
  REGISTRER_BASE_DERIVED(BodyForce, AttractingPoint);
  REGISTRER_BASE_DERIVED(BodyForce, PreferredDirection);
  REGISTRER_BASE_DERIVED(BodyForce, ViscousFluid);

  // DataExtractors
  REGISTRER_BASE_DERIVED(DataExtractor, ClusterAABB);
  REGISTRER_BASE_DERIVED(DataExtractor, dnStat);
  REGISTRER_BASE_DERIVED(DataExtractor, DuoBalance);
  REGISTRER_BASE_DERIVED(DataExtractor, MeanVelocity);
  REGISTRER_BASE_DERIVED(DataExtractor, TrackBody);
  REGISTRER_BASE_DERIVED(DataExtractor, TrackRockfall);

  // ForceLaws
  REGISTRER_BASE_DERIVED(ForceLaw, Default);
  REGISTRER_BASE_DERIVED(ForceLaw, Avalanche);
  REGISTRER_BASE_DERIVED(ForceLaw, StickedLinks);

  // PreproCommands
  REGISTRER_BASE_DERIVED(PreproCommand, copyParamsToInterfaces);
  REGISTRER_BASE_DERIVED(PreproCommand, homothetyRange);
  REGISTRER_BASE_DERIVED(PreproCommand, particlesClonage);
  REGISTRER_BASE_DERIVED(PreproCommand, randomlyOrientedVelocities);
  REGISTRER_BASE_DERIVED(PreproCommand, randomlyOrientedVelocitiesClusters);
  REGISTRER_BASE_DERIVED(PreproCommand, setAllVelocities);
  REGISTRER_BASE_DERIVED(PreproCommand, setStiffnessRatioInterfaces);
  REGISTRER_BASE_DERIVED(PreproCommand, setVariableStickParams);
  REGISTRER_BASE_DERIVED(PreproCommand, stickClusters);
  REGISTRER_BASE_DERIVED(PreproCommand, stickVerticesInClusters);
  REGISTRER_BASE_DERIVED(PreproCommand, stickVerticesInClustersMoments);

  // registerUnsharedModules();
}

/**
 *   The version of setVerboseLevel that takes a integer number as argument
 *
 *   @param v The verbose level (trace = 6, debug = 5, info = 4, warn = 3, err = 2, critical = 1, off = 0)
 */
void Rockable::setVerboseLevel(int v) {
  std::string levelNames[] = {"off", "critical", "err", "warn", "info", "debug", "trace"};
      
  std::unordered_map<int, LogLevel> levelMap = {
      {0, LogLevel::off},  
      {1, LogLevel::critical}, 
      {2, LogLevel::error},  
      {3, LogLevel::warn},
      {4, LogLevel::info}, 
      {5, LogLevel::debug},
      {6, LogLevel::trace}
  };

  if (levelMap.count(v) > 0) {
    Logger::setLevel(levelMap[v]);
    fmt::print("Verbosity level has been set to '{}'\n", levelNames[v]);
  } else {
    Logger::setLevel(LogLevel::info);
    fmt::print("The asked-level of verbosity should be in the range 0 to 6. It has been set to '{}'\n", levelNames[4]);
  }
}

/**
 *   A version of setVerboseLevel that takes a string as argument
 *
 *   @param levelName Name of the verbose level
 */
void Rockable::setVerboseLevel(const std::string& levelName) {
  std::unordered_map<std::string, int> levelMap = {{"off", 0},  {"critical", 1}, {"error", 2},  {"warn", 3},
                                                   {"info", 4}, {"debug", 5},    {"trace", 6}};

  if (levelMap.count(levelName) > 0) {
    setVerboseLevel(levelMap[levelName]);
  } else {
    fmt::print("Unknown verbosity level: '{}'\n", levelName);
  }
}

/**
 *   It opens some files (only if interactiveMode is false)
 */
void Rockable::initOutputFiles() {
  if (interactiveMode == true) return;
  perfFile.open("perf.txt");
  staticBalanceFile.open("staticBalance.txt");
  kineticEnergyFile.open("kineticEnergy.txt");
}

/**
 *   If Rockable is not used to make a simulation (in case of its usage for
 *   postprocessing) we need to set its mode as being interactive. In this case,
 *   the output files (the usual ones and the one for dataExtractors) will not be
 *   openned Also, the method 'integrate' is not usable
 */
void Rockable::setInteractive(bool imode) { interactiveMode = imode; }

/**
 * Checks if the Rockable object is in interactive mode.
 *
 * @return true if the Rockable object is in interactive mode, false otherwise
 */
bool Rockable::isInteractive() const { return interactiveMode; }

/**
 *   Print in the terminal a Banner with some information
 */
void Rockable::showBanner() {
  std::cout << "\n\n";
  std::cout << msg::bold() << msg::fg_lightBlue() << "   Rockable" << msg::fg_default() << msg::normal()
            << "  Copyright (C) 2016-2023  <vincent.richefeu@univ-grenoble-alpes.fr>\n";
  std::cout << "   This program comes with ABSOLUTELY NO WARRANTY.\n";
  std::cout << "   " << msg::bold() << msg::fg_lightBlue() << "This is academic software" << msg::fg_default()
            << msg::normal() << "\n\n";
  std::cout << "   Documentation:       install sphinx-doc\n";
  std::cout << "   e.g., for mac OS X   brew install sphinx-doc\n";
  std::cout << "                        brew link sphinx-doc --force\n";
  std::cout << "                        make html\n";
  std::cout << "                        open build/html/index.html\n\n";
  std::cout << std::endl;

#ifdef FT_CORR
  Logger::info("Compilation option: FT_CORR = \033[1;32mYES\033[0m");
#else
  Logger::info("Compilation option: FT_CORR = \033[1;31mNO\033[0m");
#endif

#ifdef COMPONENTWISE_NUM_DAMPING
  Logger::info("Compilation option: COMPONENTWISE_NUM_DAMPING = \033[1;32mYES\033[0m");
#else
  Logger::info("Compilation option: COMPONENTWISE_NUM_DAMPING = \033[1;31mNO\033[0m");
#endif

#ifdef ROCKABLE_ENABLE_PROFILING
  Logger::info("Compilation option: ROCKABLE_ENABLE_PROFILING = \033[1;32mYES\033[0m");
#else
  Logger::info("Compilation option: ROCKABLE_ENABLE_PROFILING = \033[1;31mNO\033[0m");
#endif

#ifdef ROCKABLE_ENABLE_BOUNDARY
  Logger::info("Compilation option: ROCKABLE_ENABLE_BOUNDARY = \033[1;32mYES\033[0m");
#else
  Logger::info("Compilation option: ROCKABLE_ENABLE_BOUNDARY = \033[1;31mNO\033[0m");
#endif

#ifdef ROCKABLE_ENABLE_PERIODIC
  Logger::info("Compilation option: ROCKABLE_ENABLE_PERIODIC = \033[1;32mYES\033[0m");
#else
  Logger::info("Compilation option: ROCKABLE_ENABLE_PERIODIC = \033[1;31mNO\033[0m");
#endif

#ifdef ROCKABLE_ENABLE_SOFT_PARTICLES
  Logger::info("Compilation option: ROCKABLE_ENABLE_SOFT_PARTICLES = \033[1;32mYES\033[0m");
#else
  Logger::info("Compilation option: ROCKABLE_ENABLE_SOFT_PARTICLES = \033[1;31mNO\033[0m");
#endif
}

/**
 *  Sets the number of OpenMP threads for acceleration.
 *
 *  @param nbThreads the number of threads to set
 */
void Rockable::setOpenMPThreads(int nbThreads) {
#ifdef _OPENMP
  omp_set_num_threads(nbThreads);
  Logger::info("OpenMP acceleration (Number of threads = {})", nbThreads);
#else
  Logger::info("No multithreading");
#endif
}

// ==================================================================================================================
//  CHECK METHODS
// ==================================================================================================================

/**
 * Performs initial checks for the Rockable class.
 */
void Rockable::initialChecks() {

  Logger::debug("Option forceLaw is {}", optionNames["forceLaw"]);
  Logger::debug("Option AddOrRemoveInteractions is {}", optionNames["AddOrRemoveInteractions"]);
  Logger::debug("Option UpdateNL is {}", optionNames["UpdateNL"]);
  Logger::debug("Option Integrator is {}", optionNames["Integrator"]);

  Logger::trace("DataTable size {}x{}", dataTable.ngroup, dataTable.ngroup);

  double dtc;
  estimateCriticalTimeStep(dtc);
  Logger::info("Considering a single contact between two particles, dt_critical / dt = {} (estimated)", dtc / dt);

  getCriticalTimeStep(dtc);
  if (dtc > 0.0) {
    Logger::info("dt_critical / dt = {} (over ALL Interactions)", dtc / dt);
  }

  getCurrentCriticalTimeStep(dtc);
  if (dtc > 0.0) {
    Logger::info("dt_critical / dt = {} (over ACTIVE Interactions)", dtc / dt);
  }

  // TODO: ajouter des verifs par rapport aux groupes définies et le nombre de groupes dans properties et dataTable
}

// ==================================================================================================================
//  SAVE/LOAD METHODS
// ==================================================================================================================

/**
 *   Clear the memory (exepted the shape library)
 */
void Rockable::clearMemory() {
  Particles.clear();
  Interactions.clear();
  Interfaces.clear();
  dataExtractors.clear();
  activeInteractions.clear();
  Tempos.clear();
  // Shape is not erased because it will not be read another time
  // if the filename has not been changed
}

/**
 *   Save a configuration-file named 'conf<i>'
 */
void Rockable::saveConf(int i) {
  char fname[256];
  snprintf(fname, 256, "conf%d", i);
  saveConf(fname);
}

/**
 *   Save a configuration-file
 *
 *   @param[in]  name  The name of the conf-file
 */
void Rockable::saveConf(const char* fname) {
  START_TIMER("saveConf");

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) reducedToRealKinematics();
#endif

  std::ofstream conf(fname);

  conf << "Rockable " << CONF_VERSION_DATE << '\n';  // format: progName version-date(dd-mm-yyyy)
  for (auto it = optionNames.begin(); it != optionNames.end(); ++it) {
    conf << it->first << " " << it->second << '\n';
  }
  conf << "t " << t << '\n';
  conf << "tmax " << tmax << '\n';
  conf << "dt " << dt << '\n';
  conf << "interVerlet " << interVerlet << '\n';
  conf << "interConf " << interConf << '\n';
  conf << "DVerlet " << DVerlet << '\n';
  conf << "dVerlet " << dVerlet << '\n';
  for (size_t grp = 0; grp < properties.ngroup; grp++) {
    double density = properties.get(idDensity, grp);
    if (density > 0.0) conf << "density " << grp << " " << density << '\n';
  }
  conf << "gravity " << gravity << '\n';
  if (bodyForce != nullptr) {
    conf << "BodyForce ";
    bodyForce->write(conf);
  }
  conf << "ParamsInInterfaces " << paramsInInterfaces << '\n';
  conf << "dynamicUpdateNL " << dynamicUpdateNL << '\n';
  if (dynamicUpdateNL != 0) {
    conf << "dispUpdateNL " << dispUpdateNL << '\n';
    conf << "angleUpdateNL " << Mth::rad2deg * angleUpdateNL << '\n';
  }

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {
    conf << "usePeriodicCell 1\n";
    conf << "h " << Cell.h << '\n';
    conf << "vh " << Cell.vh << '\n';
    conf << "ah " << Cell.ah << '\n';
    conf << "mh " << Cell.mass << '\n';
  } else {
    conf << "usePeriodicCell 0\n";
  }
#endif

  conf << "numericalDampingCoeff " << numericalDampingCoeff << '\n';

  conf << "VelocityBarrier " << velocityBarrier << '\n';
  conf << "AngularVelocityBarrier " << angularVelocityBarrier << '\n';

  conf << "VelocityBarrierExponent " << velocityBarrierExponent << '\n';
  conf << "AngularVelocityBarrierExponent " << angularVelocityBarrierExponent << '\n';

  if (preventCrossingLength > 0.0) conf << "preventCrossingLength " << preventCrossingLength << '\n';

  writeLawData(conf, "knContact");
  writeLawData(conf, "en2Contact");
  writeLawData(conf, "ktContact");
  writeLawData(conf, "muContact");
  writeLawData(conf, "krContact");
  writeLawData(conf, "murContact");

  writeLawData(conf, "knInnerBond");
  writeLawData(conf, "ktInnerBond");
  writeLawData(conf, "krInnerBond");
  writeLawData(conf, "fn0InnerBond");
  writeLawData(conf, "ft0InnerBond");
  writeLawData(conf, "mom0InnerBond");
  writeLawData(conf, "powInnerBond");
  writeLawData(conf, "en2InnerBond");

  writeLawData(conf, "knOuterBond");
  writeLawData(conf, "ktOuterBond");
  writeLawData(conf, "krOuterBond");
  writeLawData(conf, "fn0OuterBond");
  writeLawData(conf, "ft0OuterBond");
  writeLawData(conf, "mom0OuterBond");
  writeLawData(conf, "powOuterBond");
  writeLawData(conf, "en2OuterBond");

  conf << "ContactPartnership " << ctcPartnership.name << '\n';
  conf << "cellMinSizes " << cellMinSizes << '\n';
  conf << "boxForLinkCellsOpt " << boxForLinkCellsOpt << '\n';
  conf << "iconf " << iconf << '\n';
  conf << "nDriven " << nDriven << '\n';
  conf << "shapeFile " << shapeFile << '\n';
  conf << "precision " << CommBox().precision << '\n';
  conf << std::scientific << std::setprecision(CommBox().precision);
  conf << "Particles " << Particles.size() << '\n';

  // This commented line can help
  conf << "#name" << ' ' << "group" << ' ' << "cluster" << ' ' << "homothety" << ' ' << "pos.x" << ' ' << "pos.y" << ' '
       << "pos.z" << ' ' << "vel.x" << ' ' << "vel.y" << ' ' << "vel.z" << ' ' << "acc.x" << ' ' << "acc.y" << ' '
       << "acc.z" << ' ' << "Q.w" << ' ' << "Q.x" << ' ' << "Q.y" << ' ' << "Q.z" << ' ' << "vrot.x" << ' ' << "vrot.y"
       << ' ' << "vrot.z" << ' ' << "arot.x" << ' ' << "arot.y" << ' ' << "arot.z" << '\n';

  for (size_t i = 0; i < Particles.size(); i++) {
    conf << Particles[i].shape->name << ' ' << Particles[i].group << ' ' << Particles[i].cluster << ' '
         << Particles[i].homothety << ' ' << Particles[i].pos << ' ' << Particles[i].vel << ' ' << Particles[i].acc
         << ' ' << Particles[i].Q << ' ' << Particles[i].vrot << ' ' << Particles[i].arot;
    if (useSoftParticles == 1) {
      conf << ' ' << Particles[i].uniformTransformation << '\n';
    } else {
      conf << '\n';
    }
  }

  conf << "Interactions " << activeInteractions.size() << '\n';
  // Lexico-sort before saving
  std::sort(activeInteractions.begin(), activeInteractions.end(), std::less<Interaction*>());
  const auto prev_round = std::fegetround();
  for (size_t i = 0; i < activeInteractions.size(); i++) {
    conf << activeInteractions[i]->i << ' ' << activeInteractions[i]->j << ' ' << activeInteractions[i]->type << ' '
         << activeInteractions[i]->isub << ' ' << activeInteractions[i]->jsub << ' ';
    std::fesetround(FE_TOWARDZERO);  // value will be trunc to the given precision
    conf << activeInteractions[i]->n << ' ';
    std::fesetround(prev_round);
    conf << activeInteractions[i]->dn << ' ' << activeInteractions[i]->pos << ' ' << activeInteractions[i]->vel << ' '
         << activeInteractions[i]->fn << ' ' << activeInteractions[i]->ft << ' ' << activeInteractions[i]->mom << ' '
         << activeInteractions[i]->damp << '\n';
  }

  // Get the number of interfaces
  size_t nbInterfaces = 0;
  for (size_t i = 0; i < Interfaces.size(); i++) {
    nbInterfaces += Interfaces[i].size();
  }
  conf << "Interfaces " << nbInterfaces << std::fixed << std::setprecision(2) << '\n';
  for (size_t i = 0; i < Interfaces.size(); i++) {
    std::set<BreakableInterface>::iterator it = Interfaces[i].begin();
    for (; it != Interfaces[i].end(); ++it) {
      conf << it->i << ' ' << it->j << "  " << it->concernedBonds.size() << ' ' << it->dn0;
      if (paramsInInterfaces == 1) {
        conf << ' ' << it->kn << ' ' << it->kt << ' ' << it->kr << ' ' << it->fn0 << ' ' << it->ft0 << ' ' << it->mom0
             << ' ' << it->power << ' ';
      }
      for (size_t b = 0; b < it->concernedBonds.size(); ++b) {
        conf << "  " << it->concernedBonds[b]->type << ' ' << it->concernedBonds[b]->isub << ' '
             << it->concernedBonds[b]->jsub;
      }
      conf << '\n';
    }
  }

  conf << std::flush;

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) realToReducedKinematics();
#endif
}

/**
 *   Read law data from an input stream and set it in the data table.
 *
 *   @param is The input stream to read from.
 *   @param id The identifier for the data to be set in the data table.
 */
void Rockable::readLawData(std::istream& is, size_t id) {
  size_t g1, g2;
  double value;
  is >> g1 >> g2 >> value;
  dataTable.set(id, g1, g2, value);
}

/**
 *   Write law data for a given parameter to an output stream.
 *
 *   @param os       The output stream to write to.
 *   @param parName  The name of the parameter to write.
 */
void Rockable::writeLawData(std::ostream& os, const char* parName) {
  std::string parNameStr(parName);
  size_t id = dataTable.data_id[parNameStr];
  for (size_t g1 = 0; g1 < dataTable.ngroup; ++g1) {
    for (size_t g2 = g1; g2 < dataTable.ngroup; ++g2) {
      if (dataTable.isDefined(id, g1, g2)) {
        os << parNameStr << ' ' << g1 << ' ' << g2 << ' ' << dataTable.tables[id][g1][g2] << '\n';
      }
    }
  }
}

void Rockable::initParser() {

  parser.kwMap["t"] = __GET__(conf, t);
  parser.kwMap["tmax"] = __GET__(conf, tmax);
  parser.kwMap["dt"] = __DO__(conf) {
    conf >> dt;
    dt_2 = 0.5 * dt;
    dt2 = dt * dt;
    dt2_2 = 0.5 * dt2;
    dt2_8 = 0.125 * dt2;
    dt2_6 = dt2 / 6.0;
    dt_6 = dt / 6.0;
  };
  parser.kwMap["interVerlet"] = __GET__(conf, interVerlet);
  parser.kwMap["interConf"] = __GET__(conf, interConf);
  parser.kwMap["DVerlet"] = __GET__(conf, DVerlet);
  parser.kwMap["dVerlet"] = __GET__(conf, dVerlet);
  parser.kwMap["density"] = __DO__(conf) {
    size_t grp;
    double density;
    conf >> grp >> density;
    properties.set(idDensity, grp, density);
  };
  parser.kwMap["gravity"] = __GET__(conf, gravity);
  parser.kwMap["ParamsInInterfaces"] = __GET__(conf, paramsInInterfaces);
  parser.kwMap["dynamicUpdateNL"] = __GET__(conf, dynamicUpdateNL);
  parser.kwMap["dispUpdateNL"] = __GET__(conf, dispUpdateNL);
  parser.kwMap["angleUpdateNL"] = __DO__(conf) {
    conf >> angleUpdateNL;
    angleUpdateNL *= Mth::deg2rad;
  };

#ifdef ROCKABLE_ENABLE_PERIODIC
  parser.kwMap["usePeriodicCell"] = __GET__(conf, usePeriodicCell);

  parser.kwMap["h"] = __GET__(conf, Cell.h);
  parser.kwMap["vh"] = __GET__(conf, Cell.vh);
  parser.kwMap["ah"] = __GET__(conf, Cell.ah);
  parser.kwMap["mh"] = __GET__(conf, Cell.mass);
#endif

  parser.kwMap["useSoftParticles"] = __DO__(conf) {
    useSoftParticles = 1;
    double Y, P;
    conf >> Y >> P;
    Cinv.set(Y, P);
  };

  parser.kwMap["numericalDampingCoeff"] = __GET__(conf, numericalDampingCoeff);
  parser.kwMap["VelocityBarrier"] = __GET__(conf, velocityBarrier);
  parser.kwMap["AngularVelocityBarrier"] = __GET__(conf, angularVelocityBarrier);
  parser.kwMap["VelocityBarrierExponent"] = __GET__(conf, velocityBarrierExponent);
  parser.kwMap["AngularVelocityBarrierExponent"] = __GET__(conf, angularVelocityBarrierExponent);

  parser.kwMap["Tempo"] = [this](std::istream& conf) {
    std::string kw;
    conf >> kw;
    if (kw == "NDCoeff") {
      double t1, t2, val1, val2;
      std::string command;
      conf >> command >> t1 >> t2 >> val1 >> val2;
      Tempo<double> CundallTempo;
      Tempos.push_back(CundallTempo);
      Tempos.back().set(command, t1, t2, val1, val2);
      Tempos.back().plug(&numericalDampingCoeff);
    } else if (kw == "Inter") {
      std::string parNameStr;
      size_t g1, g2;
      conf >> parNameStr >> g1 >> g2;
      double t1, t2, val1, val2;
      std::string command;
      conf >> command >> t1 >> t2 >> val1 >> val2;
      Tempo<double> ParamTempo;
      Tempos.push_back(ParamTempo);
      Tempos.back().set(command, t1, t2, val1, val2);
      size_t id = dataTable.data_id[parNameStr];
      Tempos.back().plug(&(dataTable.tables[id][g1][g2]));
      if (g1 != g2) Tempos.back().plug(&(dataTable.tables[id][g2][g1]));
    } else if (kw == "Body") {
      std::string parNameStr;
      size_t grp;
      conf >> parNameStr >> grp;
      double t1, t2, val1, val2;
      std::string command;
      conf >> command >> t1 >> t2 >> val1 >> val2;
      Tempo<double> ParamTempo;
      Tempos.push_back(ParamTempo);
      Tempos.back().set(command, t1, t2, val1, val2);
      size_t id = properties.prop_id[parNameStr];
      Tempos.back().plug(&(properties.prop[id][grp]));
    }
  };

  parser.kwMap["forceLaw"] = __DO__(conf) {
    std::string lawName;
    conf >> lawName;

    ForceLaw* FL = Factory<ForceLaw>::Instance()->Create(lawName);
    if (FL != nullptr) {
      FL->plug(this);
      FL->init();
      forceLaw = FL;
      optionNames["forceLaw"] = lawName;
      Logger::trace("The ForceLaw named {} has been activated", lawName);
    } else {
      Logger::warn("The ForceLaw named {} is unknown! -> set to Default", lawName);
      forceLaw = Factory<ForceLaw>::Instance()->Create("Default");
      forceLaw->plug(this);
      forceLaw->init();
      optionNames["forceLaw"] = "Default";
    }
  };

  parser.kwMap["initSpringJoint"] = __DO__(conf) {
    size_t ibody, jbody;
    vec3r ipos0, jpos0;
    double stiffness;
    conf >> ibody >> ipos0 >> jbody >> jpos0 >> stiffness;
    SpringJoint SJ(ibody, jbody, ipos0, jpos0);
    SJ.init(Particles, stiffness);
    joints.push_back(SJ);
  };

  parser.kwMap["AddOrRemoveInteractions"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    setAddOrRemoveInteractions(Name);
  };
  parser.kwMap["UpdateNL"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    setUpdateNL(Name);
  };
  parser.kwMap["Integrator"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    setIntegrator(Name);
  };

  parser.kwMap["cellMinSizes"] = __GET__(conf, cellMinSizes);
  parser.kwMap["boxForLinkCellsOpt"] = __GET__(conf, boxForLinkCellsOpt);
  parser.kwMap["ContactPartnership"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    ctcPartnership.setModel(Name);
  };

  parser.kwMap["preventCrossingLength"] = __GET__(conf, preventCrossingLength);

  parser.kwMap["knContact"] = __DO__(conf) { readLawData(conf, idKnContact); };
  parser.kwMap["en2Contact"] = __DO__(conf) { readLawData(conf, idEn2Contact); };
  parser.kwMap["en2ContactFromViscRate"] = __DO__(conf) {
    size_t g1, g2;
    double viscRate;
    conf >> g1 >> g2 >> viscRate;
    double en2 = exp(-viscRate * Mth::pi / sqrt(1.0 - viscRate * viscRate));
    dataTable.set(idEn2Contact, g1, g2, en2);
  };
  parser.kwMap["ktContact"] = __DO__(conf) { readLawData(conf, idKtContact); };
  parser.kwMap["muContact"] = __DO__(conf) { readLawData(conf, idMuContact); };
  parser.kwMap["krContact"] = __DO__(conf) { readLawData(conf, idKrContact); };
  parser.kwMap["murContact"] = __DO__(conf) { readLawData(conf, idMurContact); };

  parser.kwMap["knInnerBond"] = __DO__(conf) { readLawData(conf, idKnInnerBond); };
  parser.kwMap["ktInnerBond"] = __DO__(conf) { readLawData(conf, idKtInnerBond); };
  parser.kwMap["krInnerBond"] = __DO__(conf) { readLawData(conf, idKrInnerBond); };
  parser.kwMap["fn0InnerBond"] = __DO__(conf) { readLawData(conf, idFn0InnerBond); };
  parser.kwMap["ft0InnerBond"] = __DO__(conf) { readLawData(conf, idFt0InnerBond); };
  parser.kwMap["mom0InnerBond"] = __DO__(conf) { readLawData(conf, idMom0InnerBond); };
  parser.kwMap["powInnerBond"] = __DO__(conf) { readLawData(conf, idPowInnerBond); };
  parser.kwMap["en2InnerBond"] = __DO__(conf) { readLawData(conf, idEn2InnerBond); };

  parser.kwMap["knOuterBond"] = __DO__(conf) { readLawData(conf, idKnOuterBond); };
  parser.kwMap["ktOuterBond"] = __DO__(conf) { readLawData(conf, idKtOuterBond); };
  parser.kwMap["krOuterBond"] = __DO__(conf) { readLawData(conf, idKrOuterBond); };
  parser.kwMap["fn0OuterBond"] = __DO__(conf) { readLawData(conf, idFn0OuterBond); };
  parser.kwMap["ft0OuterBond"] = __DO__(conf) { readLawData(conf, idFt0OuterBond); };
  parser.kwMap["mom0OuterBond"] = __DO__(conf) { readLawData(conf, idMom0OuterBond); };
  parser.kwMap["powOuterBond"] = __DO__(conf) { readLawData(conf, idPowOuterBond); };
  parser.kwMap["en2OuterBond"] = __DO__(conf) { readLawData(conf, idEn2OuterBond); };

  parser.kwMap["iconf"] = __GET__(conf, iconf);
  parser.kwMap["nDriven"] = __GET__(conf, nDriven);
  parser.kwMap["shapeFile"] = __DO__(conf) {
    std::string name;
    conf >> name;
    std::string wantedLib = m_path + std::string(name);
    Logger::info("wantedLib is {}", wantedLib);
    if (wantedLib != shapeFile) {  // it means that the library is not already loaded
      shapeFile = wantedLib;
      loadShapes(shapeFile.c_str());
      Logger::info("wantedLib is {}", wantedLib);
    } else {
      Logger::info("wantedLib is {} (already loaded)", wantedLib);
    }
  };
  parser.kwMap["separator"] = __DO__(conf) {  // TODO: remove it !?
    std::string keyword;
    conf >> keyword;
    CommBox().sepFromKeyword(keyword);
  };
  parser.kwMap["glue_with_walls"] = __DO__(conf) {
    std::string YesNo;
    conf >> YesNo;
    if (YesNo == "yes" || YesNo == "YES" || YesNo == "y" || YesNo == "Y" || YesNo == "1") glue_with_walls = true;
  };
  parser.kwMap["precision"] = __DO__(conf) {
    int pr;
    conf >> pr;
    CommBox().set_precision(pr);
  };
  parser.kwMap["Particles"] = __DO__(conf) {
    size_t nb;
    conf >> nb;
    Logger::info("Number of bodies: {}", nb);
    Particles.clear();
    Particle P;
    std::string shpName;
    for (size_t i = 0; i < nb; i++) {
      conf >> shpName;
      if (shpName[0] == '#') {  // a comment line within the particle list
        std::string trash;
        getline(conf, trash);
        --i;
        continue;
      }
      if (shpName[0] == '!') {  // to 'disable' a line of particle (the number of
                                // particles is decreased but you don't need to change
                                // the number of particles in the input file)
        std::string trash;
        getline(conf, trash);
        --i;
        --nb;
        continue;
      }
      conf >> P.group >> P.cluster >> P.homothety >> P.pos >> P.vel >> P.acc >> P.Q >> P.vrot >> P.arot;
      if (useSoftParticles == 1) {
        conf >> P.uniformTransformation;
      } else {
        P.uniformTransformation = mat9r::unit();
      }

      P.shape = &(Shapes[shapeId[shpName]]);  // Plug to the shape
      double h = P.homothety;
      P.mass = (h * h * h * P.shape->volume) * properties.get(idDensity, P.group);
      P.inertia = (h * h * P.shape->inertia_mass) * P.mass;

      Particles.push_back(P);
    }
    if (Interactions.size() != Particles.size()) Interactions.resize(Particles.size());
    if (Interfaces.size() != Particles.size()) Interfaces.resize(Particles.size());
  };
  parser.kwMap["Interactions"] = __DO__(conf) {
    size_t nb;
    conf >> nb;
    Logger::info("Number of interactions: {}", nb);

    if (Interactions.size() != Particles.size()) Interactions.resize(Particles.size());
    for (size_t i = 0; i < Particles.size(); ++i) Interactions[i].clear();
    activeInteractions.clear();

    Interaction I;
    for (size_t k = 0; k < nb; ++k) {
      conf >> I.i >> I.j >> I.type >> I.isub >> I.jsub >> I.n >> I.dn >> I.pos >> I.vel >> I.fn >> I.ft >> I.mom >>
          I.damp;

      // Remark: the vector 'activeInteractions' will be filled in the call of
      // 'acceleration' at the end of this method
      Interactions[I.i].insert(I);
    }
  };
  parser.kwMap["Interfaces"] = __DO__(conf) {
    size_t nbInterf;
    conf >> nbInterf;
    Logger::info("Number of interfaces: {}", nbInterf);

    if (Interfaces.size() != Particles.size()) Interfaces.resize(Particles.size());
    for (size_t i = 0; i < Particles.size(); ++i) Interfaces[i].clear();

    for (size_t k = 0; k < nbInterf; ++k) {
      Interaction ItoFind;
      BreakableInterface BI;
      size_t nbBonds;
      conf >> BI.i >> BI.j >> nbBonds >> BI.dn0;
      ItoFind.i = BI.i;
      ItoFind.j = BI.j;
      BI.concernedBonds.resize(nbBonds);

      if (paramsInInterfaces == 1) {
        conf >> BI.kn >> BI.kt >> BI.kr >> BI.fn0 >> BI.ft0 >> BI.mom0 >> BI.power;
      }

      bool missed = false;
      for (size_t u = 0; u < nbBonds; ++u) {
        conf >> ItoFind.type >> ItoFind.isub >> ItoFind.jsub;
        std::set<Interaction>::iterator it = Interactions[ItoFind.i].find(ItoFind);
        if (it != Interactions[ItoFind.i].end()) {
          BI.concernedBonds[u] = const_cast<Interaction*>(std::addressof(*it));
        } else {
          missed = true;
          Logger::warn("Cannot find interaction: i={} j={} type={} isub={} jsub={}", ItoFind.i, ItoFind.j,
                        ItoFind.type, ItoFind.isub, ItoFind.jsub);
        }
      }
      // Connection of the BreakableInterface with the Interactions
      if (missed == false) {
        BreakableInterface* ptr =
            const_cast<BreakableInterface*>(std::addressof(*((Interfaces[BI.i].insert(BI)).first)));
        for (size_t i = 0; i < BI.concernedBonds.size(); i++) {
          ptr->concernedBonds[i]->stick = ptr;
        }
      }
    }
  };
  parser.kwMap["DataExtractor"] = __DO__(conf) {  // Kept for compatibility (use file dataExtractors.txt instead)

    if (interactiveMode == true) return;  // The dataExtractors are not read in interactive mode

    std::string ExtractorName;
    conf >> ExtractorName;

    Logger::warn(
        "The DataExtractor named {} is defined in the input file! It will not be saved in conf-files\nA better "
        "solution is to put them in a file named 'dataExtractors.txt'",
        ExtractorName);

    DataExtractor* DE = Factory<DataExtractor>::Instance()->Create(ExtractorName);
    if (DE != nullptr) {
      DE->plug(this);
      DE->read(conf);
      dataExtractors.push_back(DE);
    } else {
      Logger::warn("The DataExtractor named {} is unknown", ExtractorName);
    }
  };
  parser.kwMap["BodyForce"] = __DO__(conf) {

    std::string BodyForceName;
    conf >> BodyForceName;

    BodyForce* BF = Factory<BodyForce>::Instance()->Create(BodyForceName);
    if (BF != nullptr) {
      BF->plug(this);
      BF->read(conf);
      bodyForce = BF;
      Logger::trace("The BodyForce named {} has been activated", BodyForceName);
    } else {
      Logger::warn("The BodyForce named {} is unknown!", BodyForceName);
    }
  };

  // ======== FROM HERE, THESE ARE PREPRO COMMANDS (not saved in further conf-files)
  //          They are generally put at the end of the input file
  //          so that they apply on a system already set

  const std::string commands[] = {"stickVerticesInClusters",
                                  "stickVerticesInClustersMoments",
                                  "stickClusters",
                                  "randomlyOrientedVelocities",
                                  "randomlyOrientedVelocitiesClusters",
                                  "copyParamsToInterfaces",
                                  "setStiffnessRatioInterfaces",
                                  "setVariableStickParams",
                                  "setAllVelocities",
                                  "homothetyRange",
                                  "particlesClonage"};

  for (const std::string& command : commands) {
    PreproCommand* PC = Factory<PreproCommand>::Instance()->Create(command);
    if (PC != nullptr) {
      PC->plug(this);
      PC->addCommand();
    } else {
      Logger::warn("The PreproCommand named {} was not added!", command);
    }
  }
}

/**
 *   Load a configuration-file named 'conf<i>'
 */
void Rockable::loadConf(int i) {
  char fname[256];
  snprintf(fname, 256, "conf%d", i);
  loadConf(fname);
}

/**
 *   Load a configuration file named name
 *
 *   @param[in]  name  The name of the conf-file
 */
void Rockable::loadConf(const char* a_name) {
  std::string name = m_path + std::string(a_name);
  std::ifstream conf(name);
  if (!conf.is_open()) {
    Logger::warn("@Rockable::loadConf, cannot read {}", name);
    exit(-1);
  }

  // Check header
  std::string prog;
  conf >> prog;
  if (prog != "Rockable") {
    // used for testing
    if (prog == "redirection") {
      std::string newName = std::string();  // conf name
      conf >> m_path;
      conf >> newName;
      std::cout << m_path << std::endl;
      loadConf(newName.c_str());
      return;
    }
    Logger::warn("@Rockable::loadConf, this doesn't seem to be a file for the code Rockable!");
  }
  std::string date;
  conf >> date;
  if (date != CONF_VERSION_DATE) {
    Logger::warn("@Rockable::loadConf, the version-date should be '{}'\n(in most cases, this should not be a problem)",
                  CONF_VERSION_DATE);
  }

  // This single line actually parses the file
  parser.parse(conf);

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {
    Cell.precomputeInverse();
  }
#endif

  Interactions_from_set_to_vec();

  // This is a kind of fake time-step (the time is not incremented)
  // mainly used to compute resultant forces and moments on the particles.
  // It will also populate the vector 'activeInteractions'.
  if (interactiveMode == false) {
    accelerations();
#ifdef ROCKABLE_ENABLE_PERIODIC
    if (usePeriodicCell == 1) {
      realToReducedKinematics();
    }
#endif
  }
}

/**
 *   Read the DataExtractor commands in a file named 'dataExractors.txt'.
 *   If the file is not found, or if the interactiveMode is true, it simply
 *   exit without making anything.
 */
void Rockable::readDataExtractors() {
  if (interactiveMode == true) return;  // The dataExtractors are not read in interactive mode

  if (!fileTool::fileExists("dataExtractors.txt")) {
    return;
  }

  std::ifstream is("dataExtractors.txt");

  std::string ExtractorName;
  is >> ExtractorName;
  while (is) {
    if (ExtractorName[0] == '/' || ExtractorName[0] == '#' || ExtractorName[0] == '!')
      getline(is, ExtractorName);
    else {
      DataExtractor* DE = Factory<DataExtractor>::Instance()->Create(ExtractorName);
      if (DE != nullptr) {
        DE->plug(this);
        DE->read(is);
        dataExtractors.push_back(DE);
      } else {
        Logger::warn("The DataExtractor named {} is unknown!", ExtractorName);
      }
    }
    is >> ExtractorName;
  }
}

/**
 *   Load all shapes defined in the file 'fileName'.
 *   If the library has already been loaded (ie. the file name is the same as
 *   the one in the previously read conf-file), then it is not re-read.
 */
void Rockable::loadShapes(const char* fileName) {
  // If a library file is in the running folder, so it is preferably used
  std::string ModFileName(fileName);
  std::string LocalFileName = fileTool::GetFileName(ModFileName) + "." + fileTool::GetFileExt(ModFileName);
  if (fileTool::fileExists(LocalFileName.c_str())) {
    ModFileName = LocalFileName;
  }

  if (!fileTool::fileExists(ModFileName.c_str())) {
    Logger::warn("@Rockable::loadShapes, shape library named '{}' has not been found", ModFileName);
    exit(-1);
  }

  Shapes.clear();
  shapeId.clear();

  std::ifstream is(ModFileName.c_str());

  std::string token;
  is >> token;
  while (is) {
    if (token == "<") {
      Shape S;
      S.read(is);
      if (S.preCompDone == 'n') {
        S.massProperties();
      }
      Shapes.push_back(S);
    }
    is >> token;
  }

  Logger::info("Number of shapes found in the library file {}: {}", ModFileName, Shapes.size());

  for (size_t s = 0; s < Shapes.size(); ++s) {
    Shapes[s].buildOBBtree();
    shapeId[Shapes[s].name] = s;
  }
}

/**
 * Runs the console for the Rockable class.
 *
 * @param confFileName the name of the configuration file to load
 */
void Rockable::console_run(const std::string& confFileName) {
  loadConf(confFileName.c_str());
  initOutputFiles();

  initialChecks();

  System.read(true);
  readDataExtractors();

  if (!dataExtractors.empty()) {
    std::ofstream docfile("extractedDataDoc.txt");
    for (size_t d = 0; d < dataExtractors.size(); d++) {
      dataExtractors[d]->generateHelp(docfile);
      // The initialisation of dataExtractors can be necessary after conf is loaded
      // and the System is read
      dataExtractors[d]->init();
    }
    docfile.close();
  }

  std::cout << std::endl << std::endl;

  Logger::info("INITIAL UPDATE OF NEIGHBOR LIST");

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {
    reducedToRealKinematics();
    UpdateNL();
    realToReducedKinematics();
  } else {
    UpdateNL();
  }
#else
  UpdateNL();
#endif

  Logger::info("COMPUTATION STARTS");
  integrate();
  Logger::info("COMPUTATION NORMALLY STOPPED");
}

// ==================================================================================================================
//  ADD OR REMOVE A SINGLE INTERACTION
// ==================================================================================================================

/**
 * Sets the `AddOrRemoveInteractions` function based on the given `Name`.
 *
 * @param Name the name of the function to set
 */
void Rockable::setAddOrRemoveInteractions(std::string& Name) {
  if (Name == "bruteForce") {

    AddOrRemoveInteractions = [this](size_t i, size_t j, double dmax) -> int {
      return this->AddOrRemoveInteractions_bruteForce(i, j, dmax);
    };
    optionNames["AddOrRemoveInteractions"] = "bruteForce";

  } else if (Name == "OBBtree") {

    AddOrRemoveInteractions = [this](size_t i, size_t j, double dmax) -> int {
      return this->AddOrRemoveInteractions_OBBtree(i, j, dmax);
    };
    optionNames["AddOrRemoveInteractions"] = "OBBtree";

  } else {

    Logger::warn("AddOrRemoveInteractions {} is unknown, Option remains: AddOrRemoveInteractions = {}", Name,
                  optionNames["AddOrRemoveInteractions"]);
  }
}

/**
 *   This is the brute-force O(N^2) version of the algorithm.
 *   It means that the proximity of all sub-elements of i is tested with all
 *   sub-elements of j
 */
int Rockable::AddOrRemoveInteractions_bruteForce(size_t i, size_t j, double dmax) {
  START_TIMER("AddOrRemoveInteractions_bruteForce");

  double Damp = 0.0;
  int nbAdd = 0;

  // A helper lambda function
  auto addOrRemoveSingleInteraction =
      [&](size_t i_, size_t j_, size_t isub, size_t type, size_t nbj, const vec3r& jPeriodicShift,
          std::function<bool(Particle&, Particle&, size_t, size_t, double, const vec3r&)> func) {
        Interaction to_find;
        to_find.i = i_;
        to_find.j = j_;
        to_find.type = type;
        to_find.isub = isub;
        for (size_t jsub = 0; jsub < nbj; ++jsub) {
          to_find.jsub = jsub;
          auto exist_it = (Interactions[i_]).find(to_find);
          bool NEW = (exist_it == Interactions[i_].end());
          bool NEAR = func(Particles[i_], Particles[j_], isub, jsub, dmax, jPeriodicShift);
          if (NEAR && NEW) {
            Interaction I = Interaction(i_, j_, type, isub, jsub, Damp);
            I.jPeriodicShift = jPeriodicShift;
            Interactions[i_].insert(I);
            ++nbAdd;
          } else if (!NEAR && !NEW) {
            Interactions[i_].erase(exist_it);
            --nbAdd;
          }
        }
      };

  vec3r jPeriodicShift;
#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) jPeriodicShift = Cell.getBranchCorrection(Particles[i].pos, Particles[j].pos);
#endif

  size_t nvi = Particles[i].shape->vertex.size();
  size_t nvj = Particles[j].shape->vertex.size();
  size_t nei = Particles[i].shape->edge.size();
  size_t nej = Particles[j].shape->edge.size();
  size_t nfi = Particles[i].shape->face.size();
  size_t nfj = Particles[j].shape->face.size();

  double en2 = dataTable.get(idEn2Contact, Particles[i].group, Particles[j].group);
  double kn = dataTable.get(idKnContact, Particles[i].group, Particles[j].group);

  double meff;
  if (i < nDriven)
    meff = Particles[j].mass;
  else if (j < nDriven)
    meff = Particles[i].mass;
  else
    meff = (Particles[i].mass * Particles[j].mass) / (Particles[i].mass + Particles[j].mass);
  if (en2 > 0.0 && en2 < 1.0) {
    double logen = 0.5 * log(en2);
    double dampRate = -logen / sqrt(logen * logen + Mth::piSqr);
    Damp = dampRate * 2.0 * sqrt(kn * meff);
  } else if (en2 <= 0.0)
    Damp = 2.0 * sqrt(kn * meff);
  else
    Damp = 0.0;

  // Remark: we do not put the following 3 loops into omp parallel sections
  //         because the function is called itself from a parallel for loop
  for (size_t isub = 0; isub < nvi; ++isub) {

    // First, we check whether the vertex in polyh i intersects the OBB of polyh j.
    // If not, the next calls are not required (and computation time could
    // be saved)
    OBB subObbi;
    subObbi.center = Particles[i].GlobVertex(isub);
    subObbi.enlarge(Particles[i].MinskowskiRadius() + dmax);
    OBB obbj = Particles[j].obb;
    obbj.translate(jPeriodicShift);
    obbj.enlarge(dmax);
    if (!(obbj.intersect(subObbi))) continue;

    // vertex-vertex (i, vertex isub)->(j, all vertices)
    addOrRemoveSingleInteraction(i, j, isub, vvType, nvj, jPeriodicShift, Particle::VertexIsNearVertex);

    // vertex-edge (i, vertex isub)->(j, all edges)
    addOrRemoveSingleInteraction(i, j, isub, veType, nej, jPeriodicShift, Particle::VertexIsNearEdge);

    // vertex-face (i, vertex isub)->(j, all faces)
    addOrRemoveSingleInteraction(i, j, isub, vfType, nfj, jPeriodicShift, Particle::VertexIsNearFace);

  }  // end for isub (i->j)

  for (size_t jsub = 0; jsub < nvj; ++jsub) {

    OBB subObbj;
    subObbj.center = Particles[j].GlobVertex(jsub) + jPeriodicShift;
    subObbj.enlarge(Particles[j].MinskowskiRadius() + dmax);
    OBB obbi = Particles[i].obb;
    obbi.enlarge(dmax);
    if (!(obbi.intersect(subObbj))) continue;

    // vertex-edge (j, vertex jsub)->(i, all edges)
    addOrRemoveSingleInteraction(j, i, jsub, veType, nei, -jPeriodicShift, Particle::VertexIsNearEdge);

    // vertex-face (j, vertex jsub)->(i, all faces)
    addOrRemoveSingleInteraction(j, i, jsub, vfType, nfi, -jPeriodicShift, Particle::VertexIsNearFace);

  }  // end for jsub (j->i)

  // edge-edge i->j
  for (size_t isub = 0; isub < nei; ++isub) {
    // edge-edge (i, edge isub)->(j, all edges)
    addOrRemoveSingleInteraction(i, j, isub, eeType, nej, jPeriodicShift, Particle::EdgeIsNearEdge);
  }

  return nbAdd;
}

/**
 *   This version should be best suited when the particles have a lot of sub-elements
 *   (a terrain for example). When the number of sub-elements is relatively small,
 *   it seems not to dramatically slow down the computation.
 *   REMARK: obbi and obbj need to be already placed BEFORE calling this method
 */
int Rockable::AddOrRemoveInteractions_OBBtree(size_t i, size_t j, double dmax) {
  START_TIMER("AddOrRemoveInteractions_OBBtree");

  static Interaction to_find;
  static std::set<Interaction>::iterator exist_it;

  double Damp = 0.0;
  int nbAdd = 0;

  // A helper lambda function
  auto addOrRemoveSingleInteraction =
      [&](size_t i_, size_t j_, size_t isub, size_t type, size_t jsub, const vec3r& jPeriodicShift,
          std::function<bool(Particle&, Particle&, size_t, size_t, double, const vec3r&)> func) {
        to_find.i = i_;
        to_find.j = j_;
        to_find.type = type;
        to_find.isub = isub;
        to_find.jsub = jsub;
        exist_it = (Interactions[i_]).find(to_find);
        bool NEW = (exist_it == Interactions[i_].end());
        bool NEAR = func(Particles[i_], Particles[j_], isub, jsub, dmax, jPeriodicShift);
        if (NEAR && NEW) {
          Interaction I = Interaction(i_, j_, type, isub, jsub, Damp);
          I.jPeriodicShift = jPeriodicShift;
          Interactions[i_].insert(I);
          ++nbAdd;
        } else if (!NEAR && !NEW) {
          Interactions[i_].erase(exist_it);
          --nbAdd;
        }
      };

  vec3r jPeriodicShift;
#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell) jPeriodicShift = Cell.getBranchCorrection(Particles[i].pos, Particles[j].pos);
#endif

  // Precompute the viscous damping parameter
  double en2 = dataTable.get(idEn2Contact, Particles[i].group, Particles[j].group);
  double kn = dataTable.get(idKnContact, Particles[i].group, Particles[j].group);

  double meff;
  if (i < nDriven)
    meff = Particles[j].mass;
  else if (j < nDriven)
    meff = Particles[i].mass;
  else
    meff = (Particles[i].mass * Particles[j].mass) / (Particles[i].mass + Particles[j].mass);
  if (en2 > 0.0 && en2 < 1.0) {
    double logen = 0.5 * log(en2);
    double dampRate = -logen / sqrt(logen * logen + Mth::piSqr);
    Damp = dampRate * 2.0 * sqrt(kn * meff);
  } else if (en2 <= 0.0)
    Damp = 2.0 * sqrt(kn * meff);
  else
    Damp = 0.0;

  std::vector<std::pair<subBox, subBox>> intersections;

  quat QAconj = Particles[i].Q.get_conjugated();
  vec3r posB_relativeTo_posA = QAconj * (Particles[j].pos - Particles[i].pos);
  quat QB_relativeTo_QA = QAconj * Particles[j].Q;

  OBBtree<subBox>::TreeIntersectionIds(Particles[i].shape->tree.root, Particles[j].shape->tree.root, intersections,
                                       Particles[i].homothety, Particles[j].homothety, 0.5 * dVerlet,
                                       posB_relativeTo_posA, QB_relativeTo_QA);

  for (size_t c = 0; c < intersections.size(); c++) {
    size_t isub = intersections[c].first.isub;
    size_t jsub = intersections[c].second.isub;
    int i_nbPoints = intersections[c].first.nbPoints;
    int j_nbPoints = intersections[c].second.nbPoints;

    if (i_nbPoints == 1) {
      if (j_nbPoints == 1) {  // vertex (i, isub) -> vertex (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, vvType, jsub, jPeriodicShift, Particle::VertexIsNearVertex);
      } else if (j_nbPoints == 2) {  // vertex (i, isub) -> edge (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, veType, jsub, jPeriodicShift, Particle::VertexIsNearEdge);
      } else if (j_nbPoints >= 3) {  // vertex (i, isub) -> face (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, vfType, jsub, jPeriodicShift, Particle::VertexIsNearFace);
      }
    } else if (i_nbPoints == 2) {
      if (j_nbPoints == 1) {  // vertex (j, jsub) -> edge (i, isub)
        addOrRemoveSingleInteraction(j, i, jsub, veType, isub, -jPeriodicShift, Particle::VertexIsNearEdge);
      } else if (j_nbPoints == 2) {  // vertex (i, isub) -> edge (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, eeType, jsub, jPeriodicShift, Particle::EdgeIsNearEdge);
      }
    } else if (i_nbPoints >= 3) {
      if (j_nbPoints == 1) {  // vertex (j, jsub) -> face (i, isub)
        addOrRemoveSingleInteraction(j, i, jsub, vfType, isub, -jPeriodicShift, Particle::VertexIsNearFace);
      }
    }

  }  // end loop over intersections

  return nbAdd;
}

// ==================================================================================================================
//  UPDATING THE NEIGHBOR LIST
// ==================================================================================================================

/**
 * Sets the update function for Neighbor lists, with a specified name.
 *
 * @param Name the name of the update function
 */
void Rockable::setUpdateNL(std::string& Name) {
  if (Name == "bruteForce") {

    UpdateNL = [this]() { this->UpdateNL_bruteForce(); };
    optionNames["UpdateNL"] = "bruteForce";

  } else if (Name == "linkCells") {

    UpdateNL = [this]() { this->UpdateNL_linkCells(); };
    optionNames["UpdateNL"] = "linkCells";

  } else {
    Logger::warn("UpdateNL {} is unknown, Option remains: UpdateNL = {}", Name, optionNames["UpdateNL"]);
  }
}

/**
 * Updates the dynamic check for the Rockable object's position and rotation.
 */
void Rockable::dynamicCheckUpdateNL() {
  START_TIMER("dynamicCheckUpdateNL");

  for (size_t i = 0; i < deltaPos.size(); ++i) {
    deltaPos[i] += Particles[i].vel * dt;
    double tmp1 = norm2(deltaPos[i]);
    if (tmp1 > maxDeltaPos * maxDeltaPos) maxDeltaPos = sqrt(tmp1);
    deltaQ[i] += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
    double tmp2 = deltaQ[i].get_angle();
    if (tmp2 > maxDeltaRot) maxDeltaRot = tmp2;
  }

  if (maxDeltaPos >= dispUpdateNL || maxDeltaRot >= angleUpdateNL) {
    needUpdate = true;
    for (size_t i = 0; i < deltaPos.size(); ++i) {
      deltaPos[i].reset();
      deltaQ[i].reset();
    }
    maxDeltaPos = 0.0;
    maxDeltaRot = 0.0;
  }
}

void Rockable::Interactions_from_set_to_vec() {
  START_TIMER("Interactions_from_set_to_vec");

  m_vecInteractions.resize(Interactions.size());

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < Interactions.size(); ++i) {
    m_vecInteractions[i].clear();
    for (auto& it : Interactions[i]) m_vecInteractions[i].push_back((Interaction*)&it);
  }
}

/**
 *   The most basic algorithm for building a neighbor list.
 *   That is the O(N^2) complexity
 */
void Rockable::UpdateNL_bruteForce() {
  START_TIMER("UpdateNL_bruteForce");

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].updateObb();
  }

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    BreakableInterface BI_to_find;
    BI_to_find.i = i;

    OBB obbi = Particles[i].obb;
    obbi.enlarge(0.5 * DVerlet);

    size_t jnext = i + 1;
    // Prevent interactions between driven bodies
    if (jnext < nDriven) jnext = nDriven;

    for (size_t j = jnext; j < Particles.size(); j++) {
      // if interface (i,j) is still 'active' (ie. no bond has been broken),
      // then NO contact will be possible between bodies i and j.
      BI_to_find.j = j;
      std::set<BreakableInterface>::iterator BI_it = (Interfaces[i]).find(BI_to_find);
      if (BI_it != Interfaces[i].end()) continue;  // Continue the loop (next j) if an interface is found

      OBB obbj = Particles[j].obb;
      obbj.enlarge(0.5 * DVerlet);
#ifdef ROCKABLE_ENABLE_PERIODIC
      if (usePeriodicCell) {
        obbj.center += Cell.getBranchCorrection(obbi.center, obbj.center);
      }
#endif

      // Check intersection
      if (obbi.intersect(obbj)) {
        // this is thread-safe for an access of vector Interaction in particle i
        AddOrRemoveInteractions(i, j, dVerlet);
      }
    }
  }

  Interactions_from_set_to_vec();
}

/**
 *   A broad-phase collision detection that should have a complexity of nearly O(N).
 */
void Rockable::UpdateNL_linkCells() {
  START_TIMER("UpdateNL_linkCells");

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].updateObb();
  }

  // The method computeAABB defines the AABB that surrounds the free bodies.
  // It also computes the AABB of all particles in the vector paabb.
  if (boxForLinkCellsOpt == 1) {
    computeAABB(nDriven);
  } else {  // default choice
    computeAABB();
  }

  linkCells LinkCells(aabb, cellMinSizes);
  // Populate the cells
  for (size_t i = 0; i < paabb.size(); ++i) {
    // remark: vectors paabb and Particles have the same size
    LinkCells.add_body(i, Particles[i].pos, paabb[i]);
  }

  AABB_Cell *Cc, *Cv;
  for (size_t ix = 0; ix < LinkCells.N.x; ++ix) {
    for (size_t iy = 0; iy < LinkCells.N.y; ++iy) {
      for (size_t iz = 0; iz < LinkCells.N.z; ++iz) {
        Cc = &(LinkCells.cells[ix][iy][iz]);

        for (size_t c = 0; c < Cc->pcells.size(); ++c) {

          Cv = Cc->pcells[c];  // The neighbor cells (including itself)
          for (size_t icc = 0; icc < Cc->bodies.size(); ++icc) {

            size_t i = Cc->bodies[icc];

            BreakableInterface BI_to_find;
            BI_to_find.i = i;

            OBB obbi = Particles[i].obb;
            obbi.enlarge(0.5 * DVerlet);

            for (size_t jcv = 0; jcv < Cv->bodies.size(); ++jcv) {
              size_t j = Cv->bodies[jcv];
              if (j < nDriven && i < nDriven) continue;
              if (j <= i) continue;

              // if interface (i,j) is still 'active' (ie. no bond has been broken),
              // then NO contact will be possible between bodies i and j.
              BI_to_find.j = j;
              std::set<BreakableInterface>::iterator BI_it = (Interfaces[i]).find(BI_to_find);
              if (BI_it != Interfaces[i].end()) continue;  // Continue the loop if an interface is found

              OBB obbj = Particles[j].obb;
              obbj.enlarge(0.5 * DVerlet);
#ifdef ROCKABLE_ENABLE_PERIODIC
              if (usePeriodicCell) {
                obbj.center += Cell.getBranchCorrection(obbi.center, obbj.center);
              }
#endif

              // Check intersection
              if (obbi.intersect(obbj)) {
                AddOrRemoveInteractions(i, j, dVerlet);
              }

            }  // jcv
          }    // icc

        }  // c (neighbor cells)

      }  // iz
    }    // iy
  }      // ix

  // Now we test if a too large bodies can collide another body in the
  Cc = &(LinkCells.oversized_bodies);
  for (size_t ix = 0; ix < LinkCells.N.x; ++ix) {
    for (size_t iy = 0; iy < LinkCells.N.y; ++iy) {
      for (size_t iz = 0; iz < LinkCells.N.z; ++iz) {
        Cv = &(LinkCells.cells[ix][iy][iz]);

        for (size_t icc = 0; icc < Cc->bodies.size(); ++icc) {

          size_t i = Cc->bodies[icc];

          BreakableInterface BI_to_find;
          BI_to_find.i = i;

          OBB obbi = Particles[i].obb;
          obbi.enlarge(0.5 * DVerlet);

          for (size_t jcv = 0; jcv < Cv->bodies.size(); ++jcv) {
            size_t j = Cv->bodies[jcv];

            if (j < nDriven && i < nDriven) continue;
            if (j <= i) continue;

            // if interface (i,j) is still 'active' (ie. no bond has been broken),
            // then NO contact will be possible between bodies i and j.
            BI_to_find.j = j;
            std::set<BreakableInterface>::iterator BI_it = (Interfaces[i]).find(BI_to_find);
            if (BI_it != Interfaces[i].end()) continue;  // Continue the loop if an interface is found

            OBB obbj = Particles[j].obb;
            obbj.enlarge(0.5 * DVerlet);

            // Check intersection
            if (obbi.intersect(obbj)) {
              AddOrRemoveInteractions(i, j, dVerlet);
            }

          }  // jcv
        }    // icc

      }  // iz
    }    // iy
  }      // ix

  Interactions_from_set_to_vec();
}

// ==============================================================================================================
//  INTEGRATORS
// ==============================================================================================================

/**
 * Sets the integrator for the Rockable object.
 *
 * @param Name the name of the integrator to be set
 */
void Rockable::setIntegrator(std::string& Name) {
  if (Name == "velocityVerlet") {
    integrationStep = [this]() { this->velocityVerletStep(); };
    optionNames["Integrator"] = "velocityVerlet";
  } else if (Name == "Euler") {
    integrationStep = [this]() { this->EulerStep(); };
    optionNames["Integrator"] = "Euler";
  } else if (Name == "Beeman") {
    integrationStep = [this]() { this->BeemanStep(); };
    optionNames["Integrator"] = "Beeman";
  } else if (Name == "RungeKutta4") {
    integrationStep = [this]() { this->RungeKutta4Step(); };
    optionNames["Integrator"] = "RungeKutta4";
  } else {
    Logger::warn("Integrator {} is unknown, Option remains: Integrator = {}", Name, optionNames["Integrator"]);
  }
}

/**
 * Initializes the integrator based on the selected option.
 */
void Rockable::initIntegrator() {
  if (optionNames["Integrator"] == "Beeman") {
    for (size_t i = 0; i < Particles.size(); ++i) {
      Particles[i].beemanData = std::make_shared<BeemanMoreData>();
    }
  } else if (optionNames["Integrator"] == "RungeKutta4") {
    for (size_t i = 0; i < Particles.size(); ++i) {
      Particles[i].RK4Data = std::make_shared<RK4MoreData>();
    }
  }
}

/**
 *  Updates the velocity-controlled drive of the Rockable object.
 */
void Rockable::velocityControlledDrive() {
  START_TIMER("velocityControlledDrive");
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_Vel_: {
        Particles[System.controls[c].i].vel.x = System.controls[c].value;
        Particles[System.controls[c].i].pos.x += dt * System.controls[c].value;
      } break;
      case _y_Vel_: {
        Particles[System.controls[c].i].vel.y = System.controls[c].value;
        Particles[System.controls[c].i].pos.y += dt * System.controls[c].value;
      } break;
      case _z_Vel_: {
        Particles[System.controls[c].i].vel.z = System.controls[c].value;
        Particles[System.controls[c].i].pos.z += dt * System.controls[c].value;
      } break;
      case _xrot_Vel_: {
        size_t i = System.controls[c].i;
        vec3r vrot(System.controls[c].value, 0.0, 0.0);
        Particles[i].Q += ((Particles[i].Q.dot(vrot)) *= dt);
        Particles[i].Q.normalize();
      } break;
      case _yrot_Vel_: {
        size_t i = System.controls[c].i;
        vec3r vrot(0.0, System.controls[c].value, 0.0);
        Particles[i].Q += ((Particles[i].Q.dot(vrot)) *= dt);
        Particles[i].Q.normalize();
      } break;
      case _zrot_Vel_: {
        size_t i = System.controls[c].i;
        vec3r vrot(0.0, 0.0, System.controls[c].value);
        Particles[i].Q += ((Particles[i].Q.dot(vrot)) *= dt);
        Particles[i].Q.normalize();
      } break;
      case _xyzrot_Vel_: {
        size_t i = System.controls[c].i;
        vec3r vrot = System.controls[c].vec_value;
        Particles[i].Q += ((Particles[i].Q.dot(vrot)) *= dt);
        Particles[i].Q.normalize();
      } break;
      default:
        break;
    }
  }
}

/**
 *  Makes a single step with the Euler scheme.
 */
void Rockable::EulerStep() {
  START_TIMER("Step (Euler)");

  // If a ServoFunction exists, then it is called
  if (System.ServoFunction != nullptr) System.ServoFunction(*this);

  // Controlled bodies
  velocityControlledDrive();
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.x += dt * Particles[i].vel.x;
        Particles[i].vel.x += dt * Particles[i].acc.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.y += dt * Particles[i].vel.y;
        Particles[i].vel.y += dt * Particles[i].acc.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.z += dt * Particles[i].vel.z;
        Particles[i].vel.z += dt * Particles[i].acc.z;
      } break;
      case _xrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot.x * vec3r::unit_x())) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot.x += dt * Particles[i].arot.x;
      } break;
      case _yrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot.y * vec3r::unit_y())) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot.y += dt * Particles[i].arot.y;
      } break;
      case _zrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot.z * vec3r::unit_z())) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot.z += dt * Particles[i].arot.z;
      } break;
      case _xyzrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot += dt * Particles[i].arot;
      } break;
      default:
        break;
    }
  }

  // Bodies that are free to move
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].pos += dt * Particles[i].vel;
    Particles[i].vel += dt * Particles[i].acc;

    // Rotation: Q <- Q + (dQ / dt) * dt
    // It reads like this with quaternions
    Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
    Particles[i].Q.normalize();

    Particles[i].vrot += dt * Particles[i].arot;
  }

  accelerations();
}

/**
 *  Makes a single step with the velocity-Verlet scheme.
 */
void Rockable::velocityVerletStep() {
  START_TIMER("Step (velocity-Verlet)");

  // If a ServoFunction exists, then it is called
  if (System.ServoFunction != nullptr) System.ServoFunction(*this);

  // Controlled bodies
  velocityControlledDrive();
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.x += dt * Particles[i].vel.x + dt2_2 * Particles[i].acc.x;
        Particles[i].vel.x += dt_2 * Particles[i].acc.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.y += dt * Particles[i].vel.y + dt2_2 * Particles[i].acc.y;
        Particles[i].vel.y += dt_2 * Particles[i].acc.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.z += dt * Particles[i].vel.z + dt2_2 * Particles[i].acc.z;
        Particles[i].vel.z += dt_2 * Particles[i].acc.z;
      } break;
      case _xrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot.x * vec3r::unit_x())) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot.x += dt_2 * Particles[i].arot.x;
      } break;
      case _yrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot.y * vec3r::unit_y())) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot.y += dt_2 * Particles[i].arot.y;
      } break;
      case _zrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot.z * vec3r::unit_z())) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot.z += dt_2 * Particles[i].arot.z;
      } break;
      case _xyzrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
        Particles[i].Q.normalize();
        Particles[i].vrot += dt_2 * Particles[i].arot;
      } break;
      default:
        break;
    }
  }

  // Bodies that are free to move
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].pos += dt * Particles[i].vel + dt2_2 * Particles[i].acc;
    Particles[i].vel += dt_2 * Particles[i].acc;

#ifdef ROCKABLE_ENABLE_PERIODIC
    // remember here that we use reduced coordinates here, in case of periodic cell
    if (usePeriodicCell == 1) Cell.forceToStayInside(Particles[i].pos);
#endif

    // Rotation: Q(k+1) = Q(k) + dQ(k) * dt + ddQ(k) * dt2/2
    // It reads like this with quaternions
    Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
    Particles[i].Q += ((Particles[i].Q.ddot(Particles[i].vrot, Particles[i].arot)) *= dt2_2);
    Particles[i].Q.normalize();

    Particles[i].vrot += dt_2 * Particles[i].arot;
  }

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {

    for (size_t c = 0; c < 9; c++) {  // loop over components
      if (System.cellControl.Drive[c] == ForceDriven) {
        Cell.h[c] += dt * Cell.vh[c] + dt2_2 * Cell.ah[c];
        Cell.vh[c] += dt_2 * Cell.ah[c];
      } else {
        Cell.h[c] += dt * System.cellControl.v[c];
        Cell.vh[c] = System.cellControl.v[c];
        Cell.ah[c] = 0.0;
      }
    }
    Cell.precomputeInverse();  // because Cell.h has just been updated

    reducedToRealKinematics();
    accelerations();
    realToReducedKinematics();

  } else {
    accelerations();
  }
#else
  accelerations();
#endif

  // Bodies where a component is controlled by force (like free in the corresponding direction)
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.x += dt_2 * Particles[i].acc.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.y += dt_2 * Particles[i].acc.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.z += dt_2 * Particles[i].acc.z;
      } break;
      case _xrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].vrot.x += dt_2 * Particles[i].arot.x;
      } break;
      case _yrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].vrot.y += dt_2 * Particles[i].arot.y;
      } break;
      case _zrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].vrot.z += dt_2 * Particles[i].arot.z;
      } break;
      case _xyzrot_Mom_: {
        size_t i = System.controls[c].i;
        Particles[i].vrot += dt_2 * Particles[i].arot;
      } break;
      default:
        break;
    }
  }

  // Free bodies
#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {
    vec3r vmean;
    for (size_t i = nDriven; i < Particles.size(); i++) {
      Particles[i].vel += dt_2 * Particles[i].acc;
      vmean += Particles[i].vel;
      Particles[i].vrot += dt_2 * Particles[i].arot;
    }
    vmean /= (double)(Particles.size());
    for (size_t i = 0; i < Particles.size(); i++) {
      Particles[i].vel -= vmean;
    }

    for (size_t c = 0; c < 9; c++) {
      if (System.cellControl.Drive[c] == ForceDriven) Cell.vh[c] += dt_2 * Cell.ah[c];
    }

  } else {

#pragma omp parallel for default(shared)
    for (size_t i = nDriven; i < Particles.size(); ++i) {
      Particles[i].vel += dt_2 * Particles[i].acc;
      Particles[i].vrot += dt_2 * Particles[i].arot;
    }
  }
#else

#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].vel += dt_2 * Particles[i].acc;
    Particles[i].vrot += dt_2 * Particles[i].arot;
  }

#endif
}

/**
 * Makes a single step with the Beeman scheme.
 */
void Rockable::BeemanStep() {
  START_TIMER("Step (Beeman)");

  static const double CST_2_DIV_3 = 2.0 / 3.0;
  static const double CST_1_DIV_6 = 1.0 / 6.0;
  static const double CST_1_DIV_12 = 1.0 / 12.0;
  static const double CST_5_DIV_12 = 5.0 / 12.0;

  // If a ServoFunction exists, then it is called
  if (System.ServoFunction != nullptr) System.ServoFunction(*this);

  // Controlled bodies
  velocityControlledDrive();

  accelerations();

  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].beemanData->accCurrent.x = Particles[i].acc.x;
        Particles[i].pos.x += dt * Particles[i].vel.x + dt2 * (CST_2_DIV_3 * Particles[i].acc.x -
                                                               CST_1_DIV_6 * Particles[i].beemanData->accPrevious.x);
        Particles[i].beemanData->velPrevious.x = Particles[i].vel.x;
        Particles[i].vel.x += dt * (1.5 * Particles[i].acc.x - 0.5 * Particles[i].beemanData->accPrevious.x);
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].beemanData->accCurrent.y = Particles[i].acc.y;
        Particles[i].pos.y += dt * Particles[i].vel.y + dt2 * (CST_2_DIV_3 * Particles[i].acc.y -
                                                               CST_1_DIV_6 * Particles[i].beemanData->accPrevious.y);
        Particles[i].beemanData->velPrevious.y = Particles[i].vel.y;
        Particles[i].vel.y += dt * (1.5 * Particles[i].acc.y - 0.5 * Particles[i].beemanData->accPrevious.y);
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].beemanData->accCurrent.z = Particles[i].acc.z;
        Particles[i].pos.z += dt * Particles[i].vel.z + dt2 * (CST_2_DIV_3 * Particles[i].acc.z -
                                                               CST_1_DIV_6 * Particles[i].beemanData->accPrevious.z);
        Particles[i].beemanData->velPrevious.z = Particles[i].vel.z;
        Particles[i].vel.z += dt * (1.5 * Particles[i].acc.z - 0.5 * Particles[i].beemanData->accPrevious.z);
      } break;
      // TODO imposed moments with Beeman scheme
      default:
        break;
    }
  }

  // Bodies that are free to move
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].beemanData->accCurrent = Particles[i].acc;

    Particles[i].pos += dt * Particles[i].vel +
                        dt2 * (CST_2_DIV_3 * Particles[i].acc - CST_1_DIV_6 * Particles[i].beemanData->accPrevious);
    Particles[i].beemanData->velPrevious = Particles[i].vel;
    Particles[i].vel += dt * (1.5 * Particles[i].acc - 0.5 * Particles[i].beemanData->accPrevious);

    Particles[i].beemanData->arotCurrent = Particles[i].arot;
    // Rotation: Q(k+1) = Q(k) + dQ(k) * dt + (2/3 ddQ(k) - 1/6 ddQ(k-1)) * dt2
    // It reads like this with quaternions
    Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
    quat Qacc = ((Particles[i].Q.ddot(Particles[i].vrot, Particles[i].arot)) *= CST_2_DIV_3);
    Qacc -= ((Particles[i].Q.ddot(Particles[i].vrot, Particles[i].beemanData->arotPrevious)) *= CST_1_DIV_6);
    Particles[i].Q += (Qacc *= dt2);
    Particles[i].Q.normalize();

    Particles[i].beemanData->vrotPrevious = Particles[i].vrot;
    Particles[i].vrot += dt * (1.5 * Particles[i].arot - 0.5 * Particles[i].beemanData->arotPrevious);
  }

  accelerations();

  // Bodies where a component is controlled by force (like free in the corresponding direction)
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.x =
            Particles[i].beemanData->velPrevious.x +
            dt * (CST_5_DIV_12 * Particles[i].acc.x + CST_2_DIV_3 * Particles[i].beemanData->accCurrent.x -
                  CST_1_DIV_12 * Particles[i].beemanData->accPrevious.x);
        Particles[i].beemanData->accPrevious.x = Particles[i].beemanData->accCurrent.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.y =
            Particles[i].beemanData->velPrevious.y +
            dt * (CST_5_DIV_12 * Particles[i].acc.y + CST_2_DIV_3 * Particles[i].beemanData->accCurrent.y -
                  CST_1_DIV_12 * Particles[i].beemanData->accPrevious.y);
        Particles[i].beemanData->accPrevious.y = Particles[i].beemanData->accCurrent.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.z =
            Particles[i].beemanData->velPrevious.z +
            dt * (CST_5_DIV_12 * Particles[i].acc.z + CST_2_DIV_3 * Particles[i].beemanData->accCurrent.z -
                  CST_1_DIV_12 * Particles[i].beemanData->accPrevious.z);
        Particles[i].beemanData->accPrevious.z = Particles[i].beemanData->accCurrent.z;
      } break;
      // TODO imposed moments with Beeman scheme
      default:
        break;
    }
  }

  // Free bodies
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].vel = Particles[i].beemanData->velPrevious +
                       dt * (CST_5_DIV_12 * Particles[i].acc + CST_2_DIV_3 * Particles[i].beemanData->accCurrent -
                             CST_1_DIV_12 * Particles[i].beemanData->accPrevious);
    Particles[i].vrot = Particles[i].beemanData->vrotPrevious +
                        dt * (CST_5_DIV_12 * Particles[i].arot + CST_2_DIV_3 * Particles[i].beemanData->arotCurrent -
                              CST_1_DIV_12 * Particles[i].beemanData->arotPrevious);

    Particles[i].beemanData->accPrevious = Particles[i].beemanData->accCurrent;
    Particles[i].beemanData->arotPrevious = Particles[i].beemanData->arotCurrent;
  }
}

/**
 * Makes a single step with the Runge-Kutta-Nystrom (4th order) scheme.
 */
void Rockable::RungeKutta4Step() {
  START_TIMER("Step (RK4)");

  // If a ServoFunction exists, then it is called
  if (System.ServoFunction != nullptr) System.ServoFunction(*this);

  // velocity-controlled bodies
  velocityControlledDrive();

  // Bodies that are free to move

  // k1 ........
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].RK4Data->pos0.x = Particles[i].pos.x;
        Particles[i].RK4Data->vel0.x = Particles[i].vel.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].RK4Data->pos0.y = Particles[i].pos.y;
        Particles[i].RK4Data->vel0.y = Particles[i].vel.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].RK4Data->pos0.z = Particles[i].pos.z;
        Particles[i].RK4Data->vel0.z = Particles[i].vel.z;
      } break;
      // TODO imposed moments with RK4 scheme
      default:
        break;
    }
  }
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].RK4Data->pos0 = Particles[i].pos;
    Particles[i].RK4Data->Q0 = Particles[i].Q;
    Particles[i].RK4Data->vel0 = Particles[i].vel;
    Particles[i].RK4Data->vrot0 = Particles[i].vrot;
  }
  accelerations();
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].RK4Data->k1acc = Particles[i].acc;
    Particles[i].RK4Data->k1arot = Particles[i].arot;
  }

  // k2 ........
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.x =
            Particles[i].RK4Data->pos0.x + dt_2 * Particles[i].RK4Data->vel0.x + dt2_8 * Particles[i].RK4Data->k1acc.x;
        Particles[i].vel.x = Particles[i].RK4Data->vel0.x + dt_2 * Particles[i].RK4Data->k1acc.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.y =
            Particles[i].RK4Data->pos0.y + dt_2 * Particles[i].RK4Data->vel0.y + dt2_8 * Particles[i].RK4Data->k1acc.y;
        Particles[i].vel.y = Particles[i].RK4Data->vel0.y + dt_2 * Particles[i].RK4Data->k1acc.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.z =
            Particles[i].RK4Data->pos0.z + dt_2 * Particles[i].RK4Data->vel0.z + dt2_8 * Particles[i].RK4Data->k1acc.z;
        Particles[i].vel.z = Particles[i].RK4Data->vel0.z + dt_2 * Particles[i].RK4Data->k1acc.z;
      } break;
      // TODO imposed moments with RK4 scheme
      default:
        break;
    }
  }
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    // pos0 + h_2 * vel0 + h2_8 * k1acc
    Particles[i].pos =
        Particles[i].RK4Data->pos0 + dt_2 * Particles[i].RK4Data->vel0 + dt2_8 * Particles[i].RK4Data->k1acc;
    // Q0 + h_2 * dot(Q0,vrot0) + h2_8 * ddot(Q0,vrot0,k1arot)
    Particles[i].Q = Particles[i].RK4Data->Q0;
    Particles[i].Q += (Particles[i].RK4Data->Q0.dot(Particles[i].RK4Data->vrot0) *= dt_2);
    Particles[i].Q +=
        (Particles[i].RK4Data->Q0.ddot(Particles[i].RK4Data->vrot0, Particles[i].RK4Data->k1arot) *= dt2_8);
    // vel0 + h_2 * k1acc
    Particles[i].vel = Particles[i].RK4Data->vel0 + dt_2 * Particles[i].RK4Data->k1acc;
    // vrot0 + h_2 * k1arot
    Particles[i].vrot = Particles[i].RK4Data->vrot0 + dt_2 * Particles[i].RK4Data->k1arot;
  }
  accelerations();
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].RK4Data->k2acc = Particles[i].acc;
    Particles[i].RK4Data->k2arot = Particles[i].arot;
  }

  // k3 ........
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.x = Particles[i].RK4Data->vel0.x + dt_2 * Particles[i].RK4Data->k2acc.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.y = Particles[i].RK4Data->vel0.y + dt_2 * Particles[i].RK4Data->k2acc.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].vel.z = Particles[i].RK4Data->vel0.z + dt_2 * Particles[i].RK4Data->k2acc.z;
      } break;
      // TODO imposed moments with RK4 scheme
      default:
        break;
    }
  }
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    // pos and Q not changed
    // vel0 + h_2 * k2acc
    Particles[i].vel = Particles[i].RK4Data->vel0 + dt_2 * Particles[i].RK4Data->k2acc;
    // vrot0 + h_2 * k2arot
    Particles[i].vrot = Particles[i].RK4Data->vrot0 + dt_2 * Particles[i].RK4Data->k2arot;
  }
  accelerations();
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].RK4Data->k3acc = Particles[i].acc;
    Particles[i].RK4Data->k3arot = Particles[i].arot;
  }

  // k4 ........
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.x =
            Particles[i].RK4Data->pos0.x + dt * Particles[i].RK4Data->vel0.x + dt2_2 * Particles[i].RK4Data->k3acc.x;
        Particles[i].vel.x = Particles[i].RK4Data->vel0.x + dt * Particles[i].RK4Data->k3acc.x;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.y =
            Particles[i].RK4Data->pos0.y + dt * Particles[i].RK4Data->vel0.y + dt2_2 * Particles[i].RK4Data->k3acc.y;
        Particles[i].vel.y = Particles[i].RK4Data->vel0.y + dt * Particles[i].RK4Data->k3acc.y;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.z =
            Particles[i].RK4Data->pos0.z + dt * Particles[i].RK4Data->vel0.z + dt2_2 * Particles[i].RK4Data->k3acc.z;
        Particles[i].vel.z = Particles[i].RK4Data->vel0.z + dt * Particles[i].RK4Data->k3acc.z;
      } break;
      // TODO imposed moments with RK4 scheme
      default:
        break;
    }
  }
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    // pos0 + h * vel0 + h2_2 * k3acc
    Particles[i].pos =
        Particles[i].RK4Data->pos0 + dt * Particles[i].RK4Data->vel0 + dt2_2 * Particles[i].RK4Data->k3acc;
    // Q0 + h * dot(Q0,vrot0) + h2_2 * ddot(Q0,vrot0,k3arot)
    Particles[i].Q = Particles[i].RK4Data->Q0;
    Particles[i].Q += (Particles[i].RK4Data->Q0.dot(Particles[i].RK4Data->vrot0) *= dt);
    Particles[i].Q +=
        (Particles[i].RK4Data->Q0.ddot(Particles[i].RK4Data->vrot0, Particles[i].RK4Data->k3arot) *= dt2_2);
    // vel0 + h * k3acc
    Particles[i].vel = Particles[i].RK4Data->vel0 + dt * Particles[i].RK4Data->k3acc;
    // vrot0 + h * k3arot
    Particles[i].vrot = Particles[i].RK4Data->vrot0 + dt * Particles[i].RK4Data->k3arot;
  }
  accelerations();
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].RK4Data->k4acc = Particles[i].acc;
    Particles[i].RK4Data->k4arot = Particles[i].arot;
  }

  // Weighting
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.x =
            Particles[i].RK4Data->pos0.x + dt * Particles[i].RK4Data->vel0.x +
            dt2_6 * (Particles[i].RK4Data->k1acc.x + Particles[i].RK4Data->k2acc.x + Particles[i].RK4Data->k3acc.x);
        Particles[i].vel.x =
            Particles[i].RK4Data->vel0.x + dt_6 * (Particles[i].RK4Data->k1acc.x + 2.0 * Particles[i].RK4Data->k2acc.x +
                                                   2.0 * Particles[i].RK4Data->k3acc.x + Particles[i].RK4Data->k4acc.x);
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.y =
            Particles[i].RK4Data->pos0.y + dt * Particles[i].RK4Data->vel0.y +
            dt2_6 * (Particles[i].RK4Data->k1acc.y + Particles[i].RK4Data->k2acc.y + Particles[i].RK4Data->k3acc.y);
        Particles[i].vel.y =
            Particles[i].RK4Data->vel0.y + dt_6 * (Particles[i].RK4Data->k1acc.y + 2.0 * Particles[i].RK4Data->k2acc.y +
                                                   2.0 * Particles[i].RK4Data->k3acc.y + Particles[i].RK4Data->k4acc.y);
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        Particles[i].pos.z =
            Particles[i].RK4Data->pos0.z + dt * Particles[i].RK4Data->vel0.z +
            dt2_6 * (Particles[i].RK4Data->k1acc.z + Particles[i].RK4Data->k2acc.z + Particles[i].RK4Data->k3acc.z);
        Particles[i].vel.z =
            Particles[i].RK4Data->vel0.z + dt_6 * (Particles[i].RK4Data->k1acc.z + 2.0 * Particles[i].RK4Data->k2acc.z +
                                                   2.0 * Particles[i].RK4Data->k3acc.z + Particles[i].RK4Data->k4acc.z);
      } break;
      // TODO imposed moments with RK4 scheme
      default:
        break;
    }
  }
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    // pos = pos0 + dt * vel0 + dt2_6 * (k1acc + k2acc + k3acc);
    Particles[i].pos =
        Particles[i].RK4Data->pos0 + dt * Particles[i].RK4Data->vel0 +
        dt2_6 * (Particles[i].RK4Data->k1acc + Particles[i].RK4Data->k2acc + Particles[i].RK4Data->k3acc);
    // Q = Q0 + dt * dot(Q0,vel0) + dt2_6 * (ddot(Q0,vel0,k1arot) + ddot(Q0,vel0,k2arot) + ddot(Q0,vel0,k3arot));
    Particles[i].Q = Particles[i].RK4Data->Q0;
    Particles[i].Q += (Particles[i].RK4Data->Q0.dot(Particles[i].RK4Data->vrot0) *= dt);
    quat Qadd = Particles[i].RK4Data->Q0.ddot(Particles[i].RK4Data->vrot0, Particles[i].RK4Data->k1arot);
    Qadd += Particles[i].RK4Data->Q0.ddot(Particles[i].RK4Data->vrot0, Particles[i].RK4Data->k2arot);
    Qadd += Particles[i].RK4Data->Q0.ddot(Particles[i].RK4Data->vrot0, Particles[i].RK4Data->k3arot);
    Particles[i].Q += (Qadd *= dt2_6);
    // vel = vel0 + dt_6 * (k1acc + 2.0 * k2acc + 2.0 * k3acc + k4acc);
    Particles[i].vel =
        Particles[i].RK4Data->vel0 + dt_6 * (Particles[i].RK4Data->k1acc + 2.0 * Particles[i].RK4Data->k2acc +
                                             2.0 * Particles[i].RK4Data->k3acc + Particles[i].RK4Data->k4acc);
    // vrot = vrot0 + dt_6 * (k1arot + 2.0 * k2arot + 2.0 * k3arot + k4arot);
    Particles[i].vrot =
        Particles[i].RK4Data->vrot0 + dt_6 * (Particles[i].RK4Data->k1arot + 2.0 * Particles[i].RK4Data->k2arot +
                                              2.0 * Particles[i].RK4Data->k3arot + Particles[i].RK4Data->k4arot);
  }
}

// ==============================================================================================================
//  CORE METHODS OF THE DEM ALGORITHM
// ==============================================================================================================

/**
 * Integrate the system over time.
 *
 */
void Rockable::integrate() {
  START_TIMER("Integrate");

  if (interactiveMode == true) {
    Logger::warn("It is not possible to invoke Rockable::integrate if interactiveMode is true");
    return;
  }

  if (dynamicUpdateNL != 0) {
    deltaPos.resize(Particles.size());
    deltaQ.resize(Particles.size());
    for (size_t i = 0; i < deltaPos.size(); ++i) {
      deltaPos[i].reset();
      deltaQ[i].reset();
    }
    maxDeltaPos = 0.0;
    maxDeltaRot = 0.0;
  }

  // (re)-compute some time-step constants in case they were not yet set
  // or if they have been modified for any reason
  dt_2 = 0.5 * dt;
  dt2_2 = dt_2 * dt;
  dt2 = dt * dt;
  dt2_8 = 0.125 * dt2;
  dt2_6 = dt2 / 6.0;
  dt_6 = dt / 6.0;

  Logger::trace("initIntegrator");
  initIntegrator();

  // Save the current configuration
  int frameWidth = 80;
  fmt::print("┌{0:─^{1}}┐\n", "", frameWidth);
  fmt::print("│{0: <{1}}│\n", "  Initial configuration", frameWidth);
  fmt::print("│{0: <{1}}│\n", fmt::format("  iconf: {:<6d}   Time: {:<13.8e}", iconf, t), frameWidth);
  fmt::print("└{0:─^{1}}┘\n", "", frameWidth);
  Logger::trace("start saving first conf");
  saveConf(iconf);
  Logger::trace("end saving first conf");

  PerfTimer ptimer;
  size_t step = (size_t)(t / dt);  // suppose that the time-step has not been changed
  timeInUpdateNL = 0.0;
  timeInForceComputation = 0.0;

  while (t <= tmax) {

    t += dt;
    interConfC += dt;
    interVerletC += dt;

    // It will use the selected integration scheme
    integrationStep();

    if (computationStopAsked > 0) {
      break;
    }

    if (interConfC >= interConf - dt_2) {
      iconf++;

      double elapsedTime = ptimer.getIntermediateElapsedTimeSeconds();

      // following 2 lines is a way to obtain a 1/10 precision (in percents)
      double NLPercent = std::round(1000.0 * timeInUpdateNL / elapsedTime) / 10.0;
      double ForcePercent = std::round(1000.0 * timeInForceComputation / elapsedTime) / 10.0;
      timeInUpdateNL = 0.0;
      timeInForceComputation = 0.0;

      double efficiency = interConf / elapsedTime;
      perfFile << t << ' ' << efficiency << ' ' << NLPercent << ' ' << ForcePercent << '\n' << std::flush;

      double Fmax, F_fnmax, Fmean, Fstddev;
      getResultantQuickStats(Fmax, F_fnmax, Fmean, Fstddev, nDriven);

      double fnMin, fnMax, fnMean, fnStddev;
      getInteractionQuickStats(fnMin, fnMax, fnMean, fnStddev);

      // Display a frame with some values
      fmt::print("┌{0:─^{1}}┐\n", "", frameWidth);
      fmt::print("│{0: <{1}}│\n", fmt::format("  iconf: {:<6d}   Time: {:<13.8e}", iconf, t), frameWidth);
      fmt::print(
          "│{0: <{1}}│\n",
          fmt::format("  Elapsed time since last configuration: {:<20s}", msg::HumanReadableSeconds(elapsedTime)),
          frameWidth);
      fmt::print("│{0: <{1}}│\n",
                 fmt::format("  Neighbor List: {:>6.2f}%   Forces: {:>6.2f}%   Other: {:>6.2f}% ", NLPercent,
                             ForcePercent, 100.0 - NLPercent - ForcePercent),
                 frameWidth);

#ifdef ROCKABLE_ENABLE_PERIODIC
      if (usePeriodicCell == 1) {
        fmt::print("│{0: <{1}}│\n", "  Periodic Cell:", frameWidth);
        fmt::print("│{0: <{1}}│\n",
                   fmt::format("    {:<12.4e}  {:<12.4e}  {:<12.4e} ", Cell.h.xx, Cell.h.xy, Cell.h.xz), frameWidth);
        fmt::print("│{0: <{1}}│\n",
                   fmt::format("    {:<12.4e}  {:<12.4e}  {:<12.4e} ", Cell.h.yx, Cell.h.yy, Cell.h.yz), frameWidth);
        fmt::print("│{0: <{1}}│\n",
                   fmt::format("    {:<12.4e}  {:<12.4e}  {:<12.4e} ", Cell.h.zx, Cell.h.zy, Cell.h.zz), frameWidth);

        fmt::print("│{0: <{1}}│\n", "  Stress tensor:", frameWidth);
        fmt::print("│{0: <{1}}│\n",
                   fmt::format("    {:<12.4e}  {:<12.4e}  {:<12.4e} ", Cell.Sig.xx, Cell.Sig.xy, Cell.Sig.xz),
                   frameWidth);
        fmt::print("│{0: <{1}}│\n",
                   fmt::format("    {:<12.4e}  {:<12.4e}  {:<12.4e} ", Cell.Sig.yx, Cell.Sig.yy, Cell.Sig.yz),
                   frameWidth);
        fmt::print("│{0: <{1}}│\n",
                   fmt::format("    {:<12.4e}  {:<12.4e}  {:<12.4e} ", Cell.Sig.zx, Cell.Sig.zy, Cell.Sig.zz),
                   frameWidth);
      }
#endif

      fmt::print("│{0: <{1}}│\n", "  Resultant forces on particles:", frameWidth);
      fmt::print("│{0: <{1}}│\n",
                 fmt::format("    Fmax: {:<13.8e}   Fmean: {:<13.8e}   Fstddev: {:<13.8e} ", Fmax, Fmean, Fstddev),
                 frameWidth);

      fmt::print("│{0: <{1}}│\n", "  Interaction forces:", frameWidth);
      fmt::print("│{0: <{1}}│\n", fmt::format("    fnMin: {:<13.8e}   fnMax: {:<13.8e} ", fnMin, fnMax), frameWidth);
      fmt::print("│{0: <{1}}│\n", fmt::format("    fnMean: {:<13.8e}   fnStddev: {:<13.8e} ", fnMean, fnStddev),
                 frameWidth);

      if (fnMax != 0.0 && fnMean != 0.0) {
        double Fmax_fnMean = fabs(Fmax / fnMean);
        fmt::print("│{0: <{1}}│\n", "  Static balance:", frameWidth);
        fmt::print("│{0: <{1}}│\n",
                   fmt::format("    max(F/fnMax): {:<13.8e}   Fmax/fnMean: {:<13.8e}", F_fnmax, Fmax_fnMean),
                   frameWidth);
        staticBalanceFile << t << ' ' << F_fnmax << ' ' << Fmax_fnMean << '\n' << std::flush;
      }

      double Etrans, Erot;
      getKineticEnergy(Etrans, Erot);

      fmt::print("│{0: <{1}}│\n", "  Kinetic energy:", frameWidth);
      fmt::print("│{0: <{1}}│\n", fmt::format("    Etrans: {:<13.8e}   Erot: {:<13.8e}", Etrans, Erot), frameWidth);
      kineticEnergyFile << t << ' ' << Etrans << ' ' << Erot << '\n' << std::flush;

      fmt::print("└{0:─^{1}}┘\n", "", frameWidth);

      saveConf(iconf);

      // re-read the file "drivingSystem.txt" without warning if the file does not exist.
      // Then, recompute the values imposed if a servoFunction is defined
      // REMARK: maybe a system of signal (SIGHUP) could be designed instead
      System.read(false);
      if (System.ServoFunction != nullptr) System.ServoFunction(*this);

      interConfC = 0.0;
    }

    if (needUpdate || interVerletC >= interVerlet - dt_2) {
      PerfTimer tm;

#ifdef ROCKABLE_ENABLE_PERIODIC
      if (usePeriodicCell == 1) {
        reducedToRealKinematics();
        UpdateNL();
        realToReducedKinematics();
      } else {
        UpdateNL();
      }
#else
      UpdateNL();
#endif

      timeInUpdateNL += tm.getElapsedTimeSeconds();
      if (!needUpdate)
        interVerletC = 0.0;
      else
        needUpdate = false;
    }

    // For all dataExtractors, we enventually execute some processing,
    // and then we record some data in the files
    for (size_t d = 0; d < dataExtractors.size(); ++d) {
      if (step % dataExtractors[d]->nstep == 0) dataExtractors[d]->exec();
      if (step % dataExtractors[d]->nrec == 0) dataExtractors[d]->record();
    }

    for (size_t p = 0; p < Tempos.size(); p++) {
      if (Tempos[p].update != nullptr) Tempos[p].update(t);
    }

    if (dynamicUpdateNL != 0) dynamicCheckUpdateNL();

    step++;
  }

  for (size_t d = 0; d < dataExtractors.size(); d++) dataExtractors[d]->end();

  return;
}

/**
 *  Increments the resultant forces and moments of the interacting bodies
 *  with the local forces and moments.
 */
void Rockable::incrementResultants(Interaction& I) {
  START_TIMER("incrementResultants");

  // Forces
  vec3r f = I.fn * I.n + I.ft;
  Particles[I.i].force += f;
  Particles[I.j].force -= f;

  // Moments
  vec3r Ci = (I.pos - Particles[I.i].pos);
  vec3r Cj = (I.pos - (Particles[I.j].pos + I.jPeriodicShift));
  Particles[I.i].moment += cross(Ci, f) + I.mom;
  Particles[I.j].moment += cross(Cj, -f) - I.mom;

  if (useSoftParticles == 1) {
    Particles[I.i].stress.xx += f.x * Ci.x;
    Particles[I.i].stress.xy += f.x * Ci.y;
    Particles[I.i].stress.xz += f.x * Ci.z;
    Particles[I.i].stress.yx += f.y * Ci.x;
    Particles[I.i].stress.yy += f.y * Ci.y;
    Particles[I.i].stress.yz += f.y * Ci.z;
    Particles[I.i].stress.zx += f.z * Ci.x;
    Particles[I.i].stress.zy += f.z * Ci.y;
    Particles[I.i].stress.zz += f.z * Ci.z;

    Particles[I.j].stress.xx -= f.x * Cj.x;
    Particles[I.j].stress.xy -= f.x * Cj.y;
    Particles[I.j].stress.xz -= f.x * Cj.z;
    Particles[I.j].stress.yx -= f.y * Cj.x;
    Particles[I.j].stress.yy -= f.y * Cj.y;
    Particles[I.j].stress.yz -= f.y * Cj.z;
    Particles[I.j].stress.zx -= f.z * Cj.x;
    Particles[I.j].stress.zy -= f.z * Cj.y;
    Particles[I.j].stress.zz -= f.z * Cj.z;
  }
}

#ifdef ROCKABLE_ENABLE_PERIODIC
/**
 * This function increment the moment tensor. To be used in real world
 *
 * @param I concerned Interaction
 */
void Rockable::incrementPeriodicCellTensorialMoment(Interaction& I) {
  START_TIMER("incrementPeriodicCellTensorialMoment");

  vec3r branch = (Particles[I.j].pos + I.jPeriodicShift) - Particles[I.i].pos;

  vec3r f = I.fn * I.n + I.ft;
  Cell.Sig.xx += f.x * branch.x;
  Cell.Sig.xy += f.x * branch.y;
  Cell.Sig.xz += f.x * branch.z;

  Cell.Sig.yx += f.y * branch.x;
  Cell.Sig.yy += f.y * branch.y;
  Cell.Sig.yz += f.y * branch.z;

  Cell.Sig.zx += f.z * branch.x;
  Cell.Sig.zy += f.z * branch.y;
  Cell.Sig.zz += f.z * branch.z;
}

/**
 * Convert the real coordinates of the particles to reduced coordinates.
 *
 */
void Rockable::realToReducedKinematics() {
  START_TIMER("realToReducedKinematics");

  // s = hinv . r
  // dots = hinv . (dotr - doth . s)

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    // this is supposed to be real coordinates
    vec3r r = Particles[i].pos;
    vec3r dr = Particles[i].vel;
    // vec3r ddr = Particles[i].acc;

    // now expressed in reduced coordinates
    Particles[i].pos = Cell.hinv * r;
    Particles[i].vel = Cell.hinv * (dr - Cell.vh * Particles[i].pos);
    // Particles[i].acc = Cell.hinv * ddr;
  }
}

/**
 * Transforms the reduced coordinates of the particles to real coordinates.
 *
 */
void Rockable::reducedToRealKinematics() {
  START_TIMER("reducedToRealKinematics");

  // r = h . s
  // dotr = doth . s + h . dots

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    // this is supposed to be reduced coordinates
    vec3r s = Particles[i].pos;
    vec3r ds = Particles[i].vel;
    // vec3r dds = Particles[i].acc;

    // now expressed in real coordinates
    Particles[i].pos = Cell.h * s;
    Particles[i].vel = Cell.vh * s + Cell.h * ds;
    // Particles[i].acc = Cell.h * dds;
  }
}
#endif

/**
 * Initializes the forces and moments for the Rockable class.
 *
 */
void Rockable::initialise_particle_forces_and_moments() {
  START_TIMER("initialise_particle_forces_and_moments");

  // Set resultant forces and moments to zero
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < nDriven; ++i) {
    Particles[i].force.reset();
    Particles[i].acc.reset();
    Particles[i].moment.reset();
    Particles[i].arot.reset();
    Particles[i].stress.reset();
  }
#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) Cell.Sig.reset();
#endif

#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].force = Particles[i].mass * gravity;
    Particles[i].moment.reset();
    if (bodyForce != nullptr) {
      vec3r force_more;
      vec3r moment_more;
      bodyForce->getForceAndMoment(i, force_more, moment_more);
      Particles[i].force += force_more;
      Particles[i].moment += moment_more;
    }
  }

  // Set eventually the imposed forces or moments (no body force for driven-components)
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_:
        Particles[System.controls[c].i].force.x = System.controls[c].value;
        break;
      case _y_For_:
        Particles[System.controls[c].i].force.y = System.controls[c].value;
        break;
      case _z_For_:
        Particles[System.controls[c].i].force.z = System.controls[c].value;
        break;
      case _xrot_Mom_:
        Particles[System.controls[c].i].moment.x = System.controls[c].value;
        break;
      case _yrot_Mom_:
        Particles[System.controls[c].i].moment.y = System.controls[c].value;
        break;
      case _zrot_Mom_:
        Particles[System.controls[c].i].moment.z = System.controls[c].value;
        break;
      case _xyzrot_Mom_:
        Particles[System.controls[c].i].moment = System.controls[c].vec_value;
        break;
    }
  }
}

/**
 *  Updates the interactions in the Rockable class.
 */
void Rockable::update_interactions() {
  START_TIMER("update_interactions");

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {

    // In the following loop, ALL interactions are updated,
    // including the interactions that are sticked or with positive dn
#pragma omp parallel for default(shared)
    for (size_t k = 0; k < m_vecInteractions.size(); ++k) {
      std::vector<Interaction*>& InterLoc = m_vecInteractions[k];
      for (auto it = InterLoc.begin(); it != InterLoc.end(); ++it) {
        (*it)->jPeriodicShift = Cell.getBranchCorrection(Particles[(*it)->i].pos, Particles[(*it)->j].pos);
        Interaction::UpdateDispatcherPeriodic[(*it)->type](**it, Particles[(*it)->i], Particles[(*it)->j]);
      }
    }

  } else {

    // In the following loop, ALL interactions are updated,
    // including the interactions that are sticked or with positive dn
#pragma omp parallel for default(shared)
    for (size_t k = 0; k < m_vecInteractions.size(); ++k) {
      std::vector<Interaction*>& InterLoc = m_vecInteractions[k];
      for (auto it = InterLoc.begin(); it != InterLoc.end(); ++it) {
        Interaction::UpdateDispatcher[(*it)->type](**it, Particles[(*it)->i], Particles[(*it)->j]);
      }
    }
  }
#else

// In the following loop, ALL interactions are updated,
// including the interactions that are sticked or with positive dn
#pragma omp parallel for default(shared)
  for (size_t k = 0; k < m_vecInteractions.size(); ++k) {
    std::vector<Interaction*>& InterLoc = m_vecInteractions[k];
    for (auto it = InterLoc.begin(); it != InterLoc.end(); ++it) {
      Interaction::UpdateDispatcher[(*it)->type](**it, Particles[(*it)->i], Particles[(*it)->j]);
    }
  }

#endif
}

/**
 * Updates the interactions boundary for the Rockable class.
 *
 */
void Rockable::update_interactions_boundary() {
  START_TIMER("update_interactions_boundary");

  // Remark: there can be only one Boundary
  //         For now, periodic boundary conditions cannot be used together with Boundary

#ifdef ROCKABLE_ENABLE_BOUNDARY
  activeInteractionsBoundary.clear();

#pragma omp parallel for default(shared)
  for (size_t k = 0; k < InteractionsBoundary.size(); ++k) {
    for (auto it = InteractionsBoundary[k].begin(); it != InteractionsBoundary[k].end(); ++it) {
      it->update(boundary, Particles[it->i]);
    }
  }
#endif
}

/**
 * Builds the active interactions for the Rockable class.
 *
 */
void Rockable::build_activeInteractions() {
  START_TIMER("build_activeInteractions");

  size_t part_size = 0;
  std::vector<size_t> offsets(m_vecInteractions.size(), 0);
  for (size_t k = 0; k < m_vecInteractions.size(); ++k) {
    offsets[k] = part_size;
    part_size += m_vecInteractions[k].size();
  }

  activeInteractions.resize(part_size, nullptr);

#pragma omp for schedule(static, 1)
  for (size_t k = 0; k < m_vecInteractions.size(); ++k) {
    size_t index = 0;
    std::vector<Interaction*>& InterLoc = m_vecInteractions[k];
    for (auto it = InterLoc.begin(); it != InterLoc.end(); ++it, ++index) {
      if ((*it)->dn < 0.0 || (*it)->stick != nullptr) {
        activeInteractions[offsets[k] + index] = *it;
      }
    }
  }
  activeInteractions.erase(std::remove(activeInteractions.begin(), activeInteractions.end(), nullptr),
                           activeInteractions.end());
}

/**
 * Compute the forces and moments for the Rockable object.
 *
 */
void Rockable::compute_forces_and_moments() {
  START_TIMER("compute_forces_and_moments");

  //  Now the forces and moments are computed
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < activeInteractions.size(); ++i) {
    Interaction& I = *activeInteractions[i];
    if (I.dn < 0.0 || I.stick != nullptr) {
      forceLaw->computeInteraction(I);
    }
  }
}

/**
 * Compute the resultants (forces and moments) on the bodies.
 *
 * The increment of resultants cannot be parallelised because of possible conflicts.
 * Pointer to interactions are stored in the vector 'activeInteractions'.
 *
 */
void Rockable::compute_resultants() {
  START_TIMER("compute_resultants");

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {

    for (size_t i = 0; i < activeInteractions.size(); ++i) {
      Interaction& I = *activeInteractions[i];
      if (I.dn < 0.0 || I.stick != nullptr) {
        incrementResultants(I);
        incrementPeriodicCellTensorialMoment(I);
      }
    }

  } else {

    for (size_t i = 0; i < activeInteractions.size(); ++i) {
      Interaction& I = *activeInteractions[i];
      if (I.dn < 0.0 || I.stick != nullptr) {
        incrementResultants(I);
      }
    }
  }

#else

  for (size_t i = 0; i < activeInteractions.size(); ++i) {
    Interaction& I = *activeInteractions[i];
    if (I.dn < 0.0 || I.stick != nullptr) {
      incrementResultants(I);
    }
  }

#endif
}

/**
 * Compute the SpringJoints for the Rockable object.
 *
 */
void Rockable::compute_SpringJoints() {
  START_TIMER("compute_SpringJoints");

  // SpringJoints
  for (size_t sj = 0; sj < joints.size(); ++sj) {
    vec3r forceOnj;
    joints[sj].getForceOnj(Particles, forceOnj);  // WARNING!!! this is not yet ok with periodic boundary conditions
    vec3r forceOni = -forceOnj;
    Particle* ip = &(Particles[joints[sj].ibody]);
    Particle* jp = &(Particles[joints[sj].jbody]);

    // Forces
    ip->force += forceOni;
    jp->force += forceOnj;

    // Moments
    vec3r Ci = ip->Q * joints[sj].ipos0;
    vec3r Cj = jp->Q * joints[sj].jpos0;
    ip->moment += cross(Ci, forceOni);
    jp->moment += cross(Cj, forceOnj);
  }
}

/**
 * Breaks the identified bonds in the interfaces.
 *
 */
void Rockable::breakage_of_interfaces() {
  START_TIMER("breakage_of_interfaces");

  // In this loop, all the bonds that are identified to be broken will actually be broken now
  for (std::set<BreakableInterface*>::iterator BI = interfacesToBreak.begin(); BI != interfacesToBreak.end(); ++BI) {
    std::string whichBond;
    if ((*BI)->isInner == 1)
      whichBond = "Inner";
    else
      whichBond = "Outer";
    std::cout << "Sticked (" << whichBond << ") interface between " << (*BI)->i << " and " << (*BI)->j
              << " has broken ";
    std::cout << "(" << (*BI)->concernedBonds.size() << " bonds)\n\n";

    int g1 = Particles[(*BI)->i].group;
    int g2 = Particles[(*BI)->j].group;
    for (size_t b = 0; b < (*BI)->concernedBonds.size(); ++b) {  // concernedBonds includes the current bond
      (*BI)->concernedBonds[b]->stick = nullptr;                 // it unplugs the interface from the interaction
      if ((*BI)->concernedBonds[b]->dn >= 0.0) {
        (*BI)->concernedBonds[b]->fn = 0.0;
        (*BI)->concernedBonds[b]->ft.reset();
      } else {
        // Limit the tangential force with coulomb friction
        vec3r T = (*BI)->concernedBonds[b]->ft;
        double FT = T.normalize();
        double mu = dataTable.get(idMuContact, g1, g2);
        double ft_threshold = fabs(mu * (*BI)->concernedBonds[b]->fn);
        if (FT > ft_threshold) {
          (*BI)->concernedBonds[b]->ft = ft_threshold * T;
        }
      }
    }

    // Remove the interface
    size_t BIi = (*BI)->i;
    std::set<BreakableInterface>::iterator BI_it = (Interfaces[BIi]).find(*(*BI));
    if (BI_it != Interfaces[BIi].end()) {  // if it has been found
      Interfaces[BIi].erase(BI_it);
    } else {
      size_t BIj = (*BI)->j;
      BI_it = (Interfaces[BIj]).find(*(*BI));
      if (BI_it != Interfaces[BIj].end()) {  // if it has been found
        Interfaces[BIj].erase(BI_it);
      } else
        __shouldNeverHappen;
    }

    needUpdate = true;
  }
}

/**
 * Compute the accelerations of the particles based on the resultants.
 *
 */
void Rockable::compute_accelerations_from_resultants() {
  START_TIMER("compute_accelerations_from_resultants");

#ifdef ROCKABLE_ENABLE_PERIODIC
  if (usePeriodicCell == 1) {

    double Vcell = fabs(Cell.h.det());
    Cell.Sig *= (1.0 / Vcell);

    // Acceleration of the periodic cell
    mat9r Vhinv = Vcell * Cell.hinv;
    mat9r deltaSig = Cell.Sig - System.cellControl.Sig;
    double invMass = 1.0 / Cell.mass;

    for (size_t row = 0; row < 3; row++) {
      for (size_t col = 0; col < 3; col++) {
        size_t c = 3 * row + col;
        if (System.cellControl.Drive[c] == ForceDriven) {
          Cell.ah[c] = 0.0;
          for (size_t s = 0; s < 3; s++) {
            Cell.ah[c] += Vhinv[3 * row + s] * deltaSig[3 * s + col];
          }
          Cell.ah[c] *= invMass;
        }
      }
    }
  }
#endif

  // Acceleration of controlled bodies
  for (size_t c = 0; c < System.controls.size(); ++c) {
    size_t i = System.controls[c].i;
    switch (System.controls[c].type) {
      case _x_For_:
        Particles[i].acc.x = Particles[i].force.x / Particles[i].mass;
        break;
      case _y_For_:
        Particles[i].acc.y = Particles[i].force.y / Particles[i].mass;
        break;
      case _z_For_:
        Particles[i].acc.z = Particles[i].force.z / Particles[i].mass;
        break;
      case _xrot_Mom_: {
        quat Qinv = Particles[i].Q.get_conjugated();
        vec3r vrotx(Particles[i].vrot.x, 0.0, 0.0);
        vec3r omega = Qinv * vrotx;            // Express omega in the body framework
        vec3r M = Qinv * Particles[i].moment;  // Express torque in the body framework
        vec3r domega((M[0] - (Particles[i].inertia[2] - Particles[i].inertia[1]) * omega[1] * omega[2]) /
                         Particles[i].inertia[0],
                     (M[1] - (Particles[i].inertia[0] - Particles[i].inertia[2]) * omega[2] * omega[0]) /
                         Particles[i].inertia[1],
                     (M[2] - (Particles[i].inertia[1] - Particles[i].inertia[0]) * omega[0] * omega[1]) /
                         Particles[i].inertia[2]);
        Particles[i].arot.x = (Particles[i].Q * domega).x;  // Express arot in the global framework
      } break;
      case _yrot_Mom_: {
        quat Qinv = Particles[i].Q.get_conjugated();
        vec3r vroty(0.0, Particles[i].vrot.y, 0.0);
        vec3r omega = Qinv * vroty;            // Express omega in the body framework
        vec3r M = Qinv * Particles[i].moment;  // Express torque in the body framework
        vec3r domega((M[0] - (Particles[i].inertia[2] - Particles[i].inertia[1]) * omega[1] * omega[2]) /
                         Particles[i].inertia[0],
                     (M[1] - (Particles[i].inertia[0] - Particles[i].inertia[2]) * omega[2] * omega[0]) /
                         Particles[i].inertia[1],
                     (M[2] - (Particles[i].inertia[1] - Particles[i].inertia[0]) * omega[0] * omega[1]) /
                         Particles[i].inertia[2]);
        Particles[i].arot.y = (Particles[i].Q * domega).y;  // Express arot in the global framework
      } break;
      case _zrot_Mom_: {
        quat Qinv = Particles[i].Q.get_conjugated();
        vec3r vrotz(0.0, 0.0, Particles[i].vrot.z);
        vec3r omega = Qinv * vrotz;            // Express omega in the body framework
        vec3r M = Qinv * Particles[i].moment;  // Express torque in the body framework
        vec3r domega((M[0] - (Particles[i].inertia[2] - Particles[i].inertia[1]) * omega[1] * omega[2]) /
                         Particles[i].inertia[0],
                     (M[1] - (Particles[i].inertia[0] - Particles[i].inertia[2]) * omega[2] * omega[0]) /
                         Particles[i].inertia[1],
                     (M[2] - (Particles[i].inertia[1] - Particles[i].inertia[0]) * omega[0] * omega[1]) /
                         Particles[i].inertia[2]);
        Particles[i].arot.z = (Particles[i].Q * domega).z;  // Express arot in the global framework
      } break;
      case _xyzrot_Mom_: {
        quat Qinv = Particles[i].Q.get_conjugated();
        vec3r omega = Qinv * Particles[i].vrot;  // Express omega in the body framework
        vec3r M = Qinv * Particles[i].moment;    // Express torque in the body framework
        vec3r domega((M[0] - (Particles[i].inertia[2] - Particles[i].inertia[1]) * omega[1] * omega[2]) /
                         Particles[i].inertia[0],
                     (M[1] - (Particles[i].inertia[0] - Particles[i].inertia[2]) * omega[2] * omega[0]) /
                         Particles[i].inertia[1],
                     (M[2] - (Particles[i].inertia[1] - Particles[i].inertia[0]) * omega[0] * omega[1]) /
                         Particles[i].inertia[2]);
        Particles[i].arot = (Particles[i].Q * domega);  // Express arot in the global framework
      } break;
    }
  }

// Finally compute the accelerations (translation and rotation) of the particles
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].acc = Particles[i].force / Particles[i].mass;

#ifdef ROCKABLE_ENABLE_PERIODIC
    // If there's a periodic cell, the accelerations are rescaled toward reduced coordinates
    if (usePeriodicCell == 1) Particles[i].acc = Cell.hinv * Particles[i].acc;
#endif

    quat Qinv = Particles[i].Q.get_conjugated();
    vec3r omega = Qinv * Particles[i].vrot;  // Express omega in the body framework
    vec3r M = Qinv * Particles[i].moment;    // Express torque in the body framework
    vec3r domega(
        (M[0] - (Particles[i].inertia[2] - Particles[i].inertia[1]) * omega[1] * omega[2]) / Particles[i].inertia[0],
        (M[1] - (Particles[i].inertia[0] - Particles[i].inertia[2]) * omega[2] * omega[0]) / Particles[i].inertia[1],
        (M[2] - (Particles[i].inertia[1] - Particles[i].inertia[0]) * omega[0] * omega[1]) / Particles[i].inertia[2]);
    Particles[i].arot = Particles[i].Q * domega;  // Express arot in the global framework
  }
}

/**
 *   Compute the particle accelerations (translations and rotations)
 *
 *   The method computes actually 3 things:\n
 *     1. the interaction forces and moments with the force laws,\n
 *     2. the resultant forces and moments at the body centers,\n
 *     3. the axial and angular accelerations of the bodies
 */
void Rockable::accelerations() {
  START_TIMER("Accelerations");

  initialise_particle_forces_and_moments();

  // Update the interactions (n, dn, pos and vel)
  // and apply the force law
  PerfTimer tm;
  activeInteractions.clear();
  interfacesToBreak.clear();

  update_interactions();

  // Some weighting relations can be established for the stiffnesses
  // of the cantacts that share the same body pair
  if (ctcPartnership.update != nullptr) ctcPartnership.update(*this);

  build_activeInteractions();

  compute_forces_and_moments();

  compute_resultants();

  compute_SpringJoints();

#ifdef ROCKABLE_ENABLE_SOFT_PARTICLES
  // compute_SoftParticles_transformation();
  if (useSoftParticles == 1) {
    // pour le moment on code en dure la loi. ***** DEVEL

    for (size_t i = 0; i < Particles.size(); i++) {
      double volume = /* Particles[i].uniformTransformation * */ (Particles[i].homothety * Particles[i].homothety *
                                                                  Particles[i].homothety * Particles[i].shape->volume);
      Particles[i].stress /= volume;
      Particles[i].stress.symmetrize();
      // mat9r strain = Cinv.getStrain(Particles[i].stress);
      Particles[i].uniformTransformation = mat9r::unit();  // + strain;
    }
  }
#endif

  breakage_of_interfaces();

  timeInForceComputation += tm.getElapsedTimeSeconds();

  if (numericalDampingCoeff > 0.0) numericalDamping();

  // damping solutions based on the weighting of accelerations
  if (velocityBarrier > 0.0) applyVelocityBarrier();
  if (angularVelocityBarrier > 0.0) applyAngularVelocityBarrier();

  compute_accelerations_from_resultants();
}

/**
 *   This is the so-called Cundall damping solution
 *
 *   @attention It acts on forces, so call this method after the summation of forces/moment
 *              and before the computation of accelerations
 */
void Rockable::numericalDamping() {
  START_TIMER("numericalDamping");

  double factor;
  double factorMinus = 1.0 - numericalDampingCoeff;
  double factorPlus = 1.0 + numericalDampingCoeff;

  // Damping applied to the force-driven components of the driven bodies
  for (size_t c = 0; c < System.controls.size(); ++c) {
    switch (System.controls[c].type) {
      case _x_For_: {
        size_t i = System.controls[c].i;
        factor = (Particles[i].force.x * Particles[i].vel.x > 0.0) ? factorMinus : factorPlus;
        Particles[i].force.x *= factor;
      } break;
      case _y_For_: {
        size_t i = System.controls[c].i;
        factor = (Particles[i].force.y * Particles[i].vel.y > 0.0) ? factorMinus : factorPlus;
        Particles[i].force.y *= factor;
      } break;
      case _z_For_: {
        size_t i = System.controls[c].i;
        factor = (Particles[i].force.z * Particles[i].vel.z > 0.0) ? factorMinus : factorPlus;
        Particles[i].force.z *= factor;
      } break;
      default:
        break;
    }
  }

#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {

#ifdef COMPONENTWISE_NUM_DAMPING
    // Translation
    for (int c = 0; c < 3; ++c) {
      factor = (Particles[i].force[c] * Particles[i].vel[c] > 0.0) ? factorMinus : factorPlus;
      Particles[i].force[c] *= factor;
    }

    // Rotation
    for (int c = 0; c < 3; ++c) {
      factor = (Particles[i].moment[c] * Particles[i].vrot[c] > 0.0) ? factorMinus : factorPlus;
      Particles[i].moment[c] *= factor;
    }
#else
    // Translation
    factor = (Particles[i].force * Particles[i].vel > 0.0) ? factorMinus : factorPlus;
    Particles[i].force *= factor;

    // Rotation
    factor = (Particles[i].moment * Particles[i].vrot > 0.0) ? factorMinus : factorPlus;
    Particles[i].moment *= factor;
#endif
  }
}

/**
 *   Component-wise weighting of translation acceleration to limit the velocity components to
 *   the value 'VelocityBarrier'
 *
 *   The weighting factor is equal to 1 when the velocity is of the order of zero,
 *   it tends towards 0 when the velocity approaches 'VelocityBarrier',
 *   and it becomes negative when the velocity exceeds 'VelocityBarrier'
 *   (it tends towards -1 when v tends towards infinity).
 */
void Rockable::applyVelocityBarrier() {
  START_TIMER("applyVelocityBarrier");

  for (size_t i = 0; i < Particles.size(); ++i) {
    double vxratio = pow(fabs(Particles[i].vel.x / velocityBarrier), velocityBarrierExponent);
    Particles[i].acc.x *= (1.0 - vxratio) / (1.0 + vxratio);

    double vyratio = pow(fabs(Particles[i].vel.y / velocityBarrier), velocityBarrierExponent);
    Particles[i].acc.y *= (1.0 - vyratio) / (1.0 + vyratio);

    double vzratio = pow(fabs(Particles[i].vel.z / velocityBarrier), velocityBarrierExponent);
    Particles[i].acc.z *= (1.0 - vzratio) / (1.0 + vzratio);
  }
}

/**
 *   Component-wise weighting of rotation acceleration to limit the velocity components to
 *   the value 'AngularVelocityBarrier'
 */
void Rockable::applyAngularVelocityBarrier() {
  START_TIMER("applyAngularVelocityBarrier");

  for (size_t i = 0; i < Particles.size(); ++i) {
    double vrotxratio = pow(fabs(Particles[i].vrot.x / angularVelocityBarrier), angularVelocityBarrierExponent);
    Particles[i].arot.x *= (1.0 - vrotxratio) / (1.0 + vrotxratio);

    double vrotyratio = pow(fabs(Particles[i].vrot.y / angularVelocityBarrier), angularVelocityBarrierExponent);
    Particles[i].arot.y *= (1.0 - vrotyratio) / (1.0 + vrotyratio);

    double vrotzratio = pow(fabs(Particles[i].vrot.z / angularVelocityBarrier), angularVelocityBarrierExponent);
    Particles[i].arot.z *= (1.0 - vrotzratio) / (1.0 + vrotzratio);
  }
}

// ==============================================================================================================
//  PROCESSING METHODS
// ==============================================================================================================

/**
 *   Computes the Axis Aligned Bounding Box (AABB) of a part of the scene.
 *
 *   @remark  The AABB (paabb) of each particle is also updated in this method
 *
 *   @param[in]   first   Smallest ID of particles (default value is 0)
 *   @param[in]   last    Largest ID of particles (default value corresponds to the last particle)
 */
void Rockable::computeAABB(size_t first, size_t last) {
  if (last == 0) last = Particles.size() - 1;

  paabb.clear();
  paabb.resize(Particles.size());

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); i++) {
    double radius = Particles[i].MinskowskiRadius();
    paabb[i].set_single(Particles[i].GlobVertex(0));
    for (size_t v = 1; v < Particles[i].shape->vertex.size(); v++) {
      paabb[i].add(Particles[i].GlobVertex(v));
    }
    paabb[i].enlarge(radius);
  }

  aabb = paabb[first];

  for (size_t i = first + 1; i <= last; i++) {
    aabb.enlarge(paabb[i]);
  }
}

/**
 *   Get the global kinetic energy for translations and for rotations.
 *
 *   @param[out]  Etrans  The translation kinetic energy (for particles id in range[first last])
 *   @param[out]  Erot    The ratational kinetic energy (for particles id in range[first last])
 *   @param[in]   first   Smallest ID of particles (default value is 0)
 *   @param[in]   last    Largest ID of particles (default value corresponds to the last particle)
 */
void Rockable::getKineticEnergy(double& Etrans, double& Erot, size_t first, size_t last) {
  if (last == 0) last = Particles.size() - 1;
  Etrans = 0.0;
  Erot = 0.0;
  for (size_t i = first; i <= last; i++) {
    Etrans += Particles[i].mass * (Particles[i].vel * Particles[i].vel);
    vec3r omega = Particles[i].Q.get_conjugated() * Particles[i].vrot;  // Express angular velocity in the body frame
    Erot += (Particles[i].inertia[0] * omega.x * omega.x + Particles[i].inertia[1] * omega.y * omega.y +
             Particles[i].inertia[2] * omega.z * omega.z);
  }
  Etrans *= 0.5;
  Erot *= 0.5;
}

/**
 *   Estimate the critical time step according to the stiffness values in dataTable.
 *
 *   @param[out]  dtc  minimum value of square root of meff/kn.
 */
void Rockable::estimateCriticalTimeStep(double& dtc) {
  // find the particle with the smallest mass
  // (this particle can also be the driven one)
  double massMin = 1e20;
  for (size_t i = 0; i < Particles.size(); ++i) {
    if (Particles[i].mass < massMin) massMin = Particles[i].mass;
  }

  // find the largest kn within the parameters that have been set
  // by the user.
  double knMax = -1e20;
  for (size_t igrp = 0; igrp < dataTable.ngroup; ++igrp) {
    for (size_t jgrp = 0; jgrp < dataTable.ngroup; ++jgrp) {
      if (dataTable.isDefined(idKnContact, igrp, jgrp) && dataTable.get(idKnContact, igrp, jgrp) > knMax) {
        knMax = dataTable.get(idKnContact, igrp, jgrp);
      }
      if (dataTable.isDefined(idKnInnerBond, igrp, jgrp) && dataTable.get(idKnInnerBond, igrp, jgrp) > knMax) {
        knMax = dataTable.get(idKnInnerBond, igrp, jgrp);
      }
      if (dataTable.isDefined(idKnOuterBond, igrp, jgrp) && dataTable.get(idKnOuterBond, igrp, jgrp) > knMax) {
        knMax = dataTable.get(idKnOuterBond, igrp, jgrp);
      }
    }
  }

  dtc = M_PI * sqrt(massMin / knMax);
}

/**
 *   Compute the critical time step by looping over all potential
 *   interactions, even those that are not active.
 *
 *   @param[out]  dtc  minimum value of square root of meff/kn.
 */
void Rockable::getCriticalTimeStep(double& dtc) {
  dtc = 0.0;
  double sqrdtcMin = 1.0e20;
  double sqrdtc = 0.0;
  bool okay = false;

  // remark: Interaction.size() equals Particles.size()
  //         so an empty list of interaction can not rapidly be tested.
  //         We use the variable 'okay' for that reason
  for (size_t k = 0; k < Interactions.size(); ++k) {
    for (auto it = Interactions[k].begin(); it != Interactions[k].end(); ++it) {
      int igrp = Particles[it->i].group;
      int jgrp = Particles[it->j].group;
      if (!(dataTable.isDefined(idKnContact, igrp, jgrp) || dataTable.isDefined(idKnInnerBond, igrp, jgrp) ||
            dataTable.isDefined(idKnOuterBond, igrp, jgrp)))
        continue;

      double meff;
      if (it->i < nDriven)
        meff = Particles[it->j].mass;
      else if (it->j < nDriven)
        meff = Particles[it->i].mass;
      else
        meff = (Particles[it->i].mass * Particles[it->j].mass) / (Particles[it->i].mass + Particles[it->j].mass);

      double kn;
      if (it->stick != nullptr) {
        if (Particles[it->i].cluster == Particles[it->j].cluster)
          kn = dataTable.get(idKnInnerBond, Particles[it->i].group, Particles[it->j].group);
        else
          kn = dataTable.get(idKnOuterBond, Particles[it->i].group, Particles[it->j].group);
      } else
        kn = dataTable.get(idKnContact, Particles[it->i].group, Particles[it->j].group);
      sqrdtc = meff / kn;

      if (sqrdtc < sqrdtcMin) {
        sqrdtcMin = sqrdtc;
        okay = true;
      }
    }
  }

  // compute the critical time step
  if (okay)
    dtc = M_PI * sqrt(sqrdtcMin);
  else
    dtc = -1.0;
}

/**
 *   Compute the critical time step by looping over all currently active interactions
 *
 *   @param  dtc  minimum value of square root of meff/kn.
 */
void Rockable::getCurrentCriticalTimeStep(double& dtc) {
  if (activeInteractions.empty()) {
    dtc = -1.0;
    return;
  }

  dtc = 0.0;
  double sqrdtcMin = 1.0e20;

  for (size_t i = 0; i < activeInteractions.size(); i++) {
    Interaction* it = activeInteractions[i];
    double meff;
    if (it->i < nDriven)
      meff = Particles[it->j].mass;
    else if (it->j < nDriven)
      meff = Particles[it->i].mass;
    else
      meff = (Particles[it->i].mass * Particles[it->j].mass) / (Particles[it->i].mass + Particles[it->j].mass);

    double kn;
    if (it->stick != nullptr) {
      if (Particles[it->i].cluster == Particles[it->j].cluster)
        kn = dataTable.get(idKnInnerBond, Particles[it->i].group, Particles[it->j].group);
      else
        kn = dataTable.get(idKnOuterBond, Particles[it->i].group, Particles[it->j].group);
    } else
      kn = dataTable.get(idKnContact, Particles[it->i].group, Particles[it->j].group);
    double sqrdtc = meff / kn;

    if (sqrdtc < sqrdtcMin) sqrdtcMin = sqrdtc;
  }

  // compute the critical time step
  dtc = M_PI * sqrt(sqrdtcMin);
}

/**
 * Calculates the resultant quick statistics for a range of particles in the Rockable class.
 *
 * @param Fmax    output parameter to store the maximum resultant force
 * @param F_fnmax output parameter to store the maximum resultant force normalized by the maximum interaction force
 * @param Fmean   output parameter to store the mean resultant force
 * @param Fstddev output parameter to store the standard deviation of resultant forces
 * @param first   the index of the first particle in the range
 * @param last    the index of the last particle in the range. If 0, the last particle is the last in the Particles
 * vector
 *
 * @throws None
 */
void Rockable::getResultantQuickStats(double& Fmax, double& F_fnmax, double& Fmean, double& Fstddev, size_t first,
                                      size_t last) {
  std::vector<double> fnMax(Particles.size(), 0.0);
  std::vector<size_t> nbCtc(Particles.size(), 0);
  for (size_t k = 0; k < activeInteractions.size(); k++) {
    size_t i = activeInteractions[k]->i;
    size_t j = activeInteractions[k]->j;
    nbCtc[i] += 1;
    nbCtc[j] += 1;
    double fn = activeInteractions[k]->fn;
    if (fnMax[i] < fn) fnMax[i] = fn;
    if (fnMax[j] < fn) fnMax[j] = fn;
  }

  if (last == 0) last = Particles.size() - 1;
  double F = norm(Particles[first].force);
  Fmax = F;
  F_fnmax = 0.0;
  Fmean = 0.0;
  size_t n = 0;
  for (size_t i = first; i <= last; i++) {
    if (nbCtc[i] <= 1) continue;
    n += 1;
    F = norm(Particles[i].force);
    Fmean += F;
    if (F > Fmax) Fmax = F;
    if (fnMax[i] > 0.0) {
      F /= fnMax[i];
      if (F > F_fnmax) F_fnmax = F;
    }
  }

  // Computation of the standard deviation of resultant forces
  Fstddev = 0.0;
  if (n > 1) {
    Fmean /= (double)n;
    for (size_t i = first; i <= last; i++) {
      if (nbCtc[i] <= 1) continue;
      F = norm(Particles[i].force);
      F = F - Fmean;
      F = F * F;
      Fstddev += F;
    }
    Fstddev /= (double)(n - 1);  // remark: done only if (n >= 2)
    Fstddev = sqrt(Fstddev);
  }
}

/**
 *   Get some statistics of the normal interaction forces
 *
 *   @param fnMin     Smallest normal contact force
 *   @param fnMin     Biggest normal contact force
 *   @param fnMean    Averaged value of normal contact forces
 *   @param fnStddev  Standard deviation of normal contact forces
 *
 *   @remarks Notice that this method needs that 'activeInteractions' is not empty.
 */
void Rockable::getInteractionQuickStats(double& fnMin, double& fnMax, double& fnMean, double& fnStddev) {
  if (activeInteractions.empty()) {
    fnMin = fnMax = fnMean = fnStddev = 0.0;
    return;
  }

  double fn = activeInteractions[0]->fn;
  fnMax = fnMin = fn;
  fnMean = 0.0;
  size_t n = 1;
  for (size_t i = 1; i < activeInteractions.size(); i++) {
    n += 1;
    fn = activeInteractions[i]->fn;
    fnMean += fn;
    if (fn < fnMin) fnMin = fn;
    if (fn > fnMax) fnMax = fn;
  }
  if (n > 1) {
    fnMean /= (double)(n - 1);
  }

  fnStddev = 0.0;
  if (n >= 2) {
    for (size_t i = 0; i < activeInteractions.size(); i++) {
      fn = activeInteractions[i]->fn;
      fn = fn - fnMean;
      fnStddev += fn * fn;
    }
    fnStddev /= (double)(n - 1);  // remark: done only if (n >= 2)
    fnStddev = sqrt(fnStddev);
  }
}
