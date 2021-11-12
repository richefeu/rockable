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

#define CONF_VERSION_DATE "29-11-2018"
#include "Rockable.hpp"

// ==============================================================================================================
//  INITIALIZATIONS
// ==============================================================================================================

Rockable::Rockable() {
  // Some default values (actually, most of them will be reset after)
  t = 0.0;
  tmax = 1.0;
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

  dynamicUpdateNL = 0;
  dispUpdateNL = 1.0;
  angleUpdateNL = 1.0;

  numericalDampingCoeff = 0.0;
  VelocityBarrier = 0.0;
  AngularVelocityBarrier = 0.0;
  VelocityBarrierExponent = 1.0;
  AngularVelocityBarrierExponent = 1.0;

  ParamsInInterfaces = 0;
  idDensity = properties.add("density");

  forceLawPtr = std::bind(&Rockable::forceLawAvalanches, this, std::placeholders::_1);
  optionNames["forceLaw"] = "Default";

  AddOrRemoveInteractions = std::bind(&Rockable::AddOrRemoveInteractions_bruteForce, this, std::placeholders::_1,
                                      std::placeholders::_2, std::placeholders::_3);
  optionNames["AddOrRemoveInteractions"] = "bruteForce";

  UpdateNL = std::bind(&Rockable::UpdateNL_bruteForce, this);
  optionNames["UpdateNL"] = "bruteForce";
  cellMinSizes.set(1.0, 1.0, 1.0);
  boxForLinkCellsOpt = 0;

  IntegrationStep = std::bind(&Rockable::velocityVerletStep, this);
  optionNames["Integrator"] = "velocityVerlet";

  idKnContact = dataTable.add("knContact");
  idEn2Contact = dataTable.add("en2Contact");
  idKtContact = dataTable.add("ktContact");
  idMuContact = dataTable.add("muContact");
  idKrContact = dataTable.add("krContact");
  idMurContact = dataTable.add("murContact");

  idKnInnerBond = dataTable.add("knInnerBond");
  idKtInnerBond = dataTable.add("ktInnerBond");
  idEn2InnerBond = dataTable.add("en2InnerBond");
  idFn0InnerBond = dataTable.add("fn0InnerBond");
  idFt0InnerBond = dataTable.add("ft0InnerBond");
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
  verbose = 0;
  glue_with_walls = false;
}

void Rockable::setVerboseLevel(int v) {
  if (v > 0) {
    verbose = v;
  } else {
    verbose = 0;
  }
}

void Rockable::initOutputFiles() {
  if (interactiveMode == true) return;
  perfFile.open("perf.txt");
  staticBalanceFile.open("staticBalance.txt");
  kineticEnergyFile.open("kineticEnergy.txt");
}

/**
    @brief  If Rockable is not used to make a simulation (in case of its usage for
            postprocessing) we need to set its mode as being interactive. In this case,
            the output files (the usual ones and the one for dataExtractors) will not be
            openned Also, the method 'integrate' is not usable
*/
void Rockable::setInteractive(bool imode) { interactiveMode = imode; }

/**
    @return interactiveMode
*/
bool Rockable::isInteractive() const { return interactiveMode; }

/**
    @brief Print in the terminal a Banner with some information
*/
void Rockable::showBanner() {
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "Rockable  Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>\n";
  std::cout << "This program comes with ABSOLUTELY NO WARRANTY.\n";
  std::cout << "This is academic software\n";
  std::cout << "Documentation: "
               "https://richefeu.gitbook.io/cdm/\n\n";

  std::cout << std::endl;
  std::cout << "Compilation options:\n";

#ifdef QUAT_ACC
  std::cout << "    QUAT_ACC.............YES\n";
#else
  std::cout << "    QUAT_ACC.............NO\n";
#endif

#ifdef FT_CORR
  std::cout << "    FT_CORR..............YES\n";
#else
  std::cout << "    FT_CORR..............NO\n";
#endif

#ifdef ROT_MATRIX
  std::cout << "    ROT_MATRIX...........YES\n";
#else
  std::cout << "    ROT_MATRIX...........NO\n";
#endif

#ifdef BREAK_ONCE
  std::cout << "    BREAK_ONCE...........YES (** SHOULD NOT BE USED! **)\n";
#else
  std::cout << "    BREAK_ONCE...........NO\n";
#endif

  std::cout << std::endl;
}

// ==================================================================================================================
//  CHECK METHODS
// ==================================================================================================================

void Rockable::initialChecks() {

  if (verbose >= 1) {
    std::cout << "\n OPTIONS:\n";
    std::cout << "forceLaw                = " << optionNames["forceLaw"] << '\n';
    std::cout << "AddOrRemoveInteractions = " << optionNames["AddOrRemoveInteractions"] << '\n';
    std::cout << "UpdateNL                = " << optionNames["UpdateNL"] << '\n';
    std::cout << "Integrator              = " << optionNames["Integrator"] << '\n';
  }

  if (verbose >= 2) {
    std::cout << "\n DataTable size " << dataTable.ngroup << "x" << dataTable.ngroup << " for parameters:\n";
    for (auto& d : dataTable.data_id) {
      std::cout << "      " << d.first << '\n';
    }
  }

  // Always displayed:
  double dtc;
  estimateCriticalTimeStep(dtc);
  std::cout << msg::info()
            << "Considering a single contact between two particles\n"
               "        (this can be done by using a 'ContactPartnership')\n"
               "  dt_critical / dt = "
            << dtc / dt << " (estimated)" << std::endl;

  getCriticalTimeStep(dtc);
  if (dtc > 0.0) {
    std::cout << "  dt_critical / dt = " << dtc / dt << " (over ALL Interactions)" << std::endl;
  }

  getCurrentCriticalTimeStep(dtc);
  if (dtc > 0.0) {
    std::cout << "  dt_critical / dt = " << dtc / dt << " (over ACTIVE Interactions)" << msg::normal() << std::endl;
  }
}

// ==================================================================================================================
//  SAVE/LOAD METHODS
// ==================================================================================================================

/**
    @bried Clear the memory (exepted the shape library)
*/
void Rockable::clearMemory() {
  Particles.clear();
  Interactions.clear();
  Interfaces.clear();
  dataExtractors.clear();
  activeInteractions.clear();
  Tempos.clear();
  // Shape is not erased because it will not be reread
  // if the filename has not been changed
}

/**
    @brief Save a configuration-file named 'conf<i>'
*/
void Rockable::saveConf(int i) {
  char fname[256];
  sprintf(fname, "conf%d", i);
  saveConf(fname);
}

/**
    @brief Save a configuration-file
    @param[in]  name  The name of the conf-file
*/
void Rockable::saveConf(const char* fname) {
  //char fname[256];
  //sprintf(fname, "conf%d", i);
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
    if (density > 0.0)  // TODO implement properties.isDefined(idDensity, grp)
      conf << "density " << grp << " " << density << '\n';
  }
  conf << "gravity " << gravity << '\n';
  if (bodyForce != nullptr) {
    conf << "BodyForce ";
    bodyForce->write(conf);
  }
  conf << "ParamsInInterfaces " << ParamsInInterfaces << '\n';
  conf << "dynamicUpdateNL " << dynamicUpdateNL << '\n';
  if (dynamicUpdateNL != 0) {
    conf << "dispUpdateNL " << dispUpdateNL << '\n';
    conf << "angleUpdateNL " << Mth::rad2deg * angleUpdateNL << '\n';
  }
  conf << "numericalDampingCoeff " << numericalDampingCoeff << '\n';

  conf << "VelocityBarrier " << VelocityBarrier << '\n';
  conf << "AngularVelocityBarrier " << AngularVelocityBarrier << '\n';

  conf << "VelocityBarrierExponent " << VelocityBarrierExponent << '\n';
  conf << "AngularVelocityBarrierExponent " << AngularVelocityBarrierExponent << '\n';

  writeLawData(conf, "knContact");
  writeLawData(conf, "en2Contact");
  writeLawData(conf, "ktContact");
  writeLawData(conf, "muContact");
  writeLawData(conf, "krContact");
  writeLawData(conf, "murContact");

  writeLawData(conf, "knInnerBond");
  writeLawData(conf, "ktInnerBond");
  writeLawData(conf, "fn0InnerBond");
  writeLawData(conf, "ft0InnerBond");
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
  if (CommBox().sep != ' ') conf << "separator " << CommBox().keywordFromSep() << '\n';
  conf << "precision " << CommBox().precision << '\n';
  conf << std::scientific << std::setprecision(CommBox().precision);
  conf << "Particles " << Particles.size() << '\n';
  conf << "#name" << CommBox().sep << "group" << CommBox().sep << "cluster" << CommBox().sep << "homothety"
       << CommBox().sep << "pos.x" << CommBox().sep << "pos.y" << CommBox().sep << "pos.z" << CommBox().sep << "vel.x"
       << CommBox().sep << "vel.y" << CommBox().sep << "vel.z" << CommBox().sep << "acc.x" << CommBox().sep << "acc.y"
       << CommBox().sep << "acc.z" << CommBox().sep << "Q.w" << CommBox().sep << "Q.x" << CommBox().sep << "Q.y"
       << CommBox().sep << "Q.z" << CommBox().sep << "vrot.x" << CommBox().sep << "vrot.y" << CommBox().sep << "vrot.z"
       << CommBox().sep << "arot.x" << CommBox().sep << "arot.y" << CommBox().sep << "arot.z" << '\n';
  for (size_t i = 0; i < Particles.size(); i++) {
    conf << Particles[i].shape->name << CommBox().sep << Particles[i].group << CommBox().sep << Particles[i].cluster
         << CommBox().sep << Particles[i].homothety << CommBox().sep << Particles[i].pos << CommBox().sep
         << Particles[i].vel << CommBox().sep << Particles[i].acc << CommBox().sep << Particles[i].Q << CommBox().sep
         << Particles[i].vrot << CommBox().sep << Particles[i].arot << '\n';
  }

  conf << "Interactions " << activeInteractions.size() << '\n';
  // Lexico-sort before saving
  std::sort(activeInteractions.begin(), activeInteractions.end(), std::less<Interaction*>());
  const auto prev_round = std::fegetround();
  for (size_t i = 0; i < activeInteractions.size(); i++) {
    conf << activeInteractions[i]->i << CommBox().sep << activeInteractions[i]->j << CommBox().sep
         << activeInteractions[i]->type << CommBox().sep << activeInteractions[i]->isub << CommBox().sep
         << activeInteractions[i]->jsub << CommBox().sep;
    std::fesetround(FE_TOWARDZERO);  // value will be trunc to the given precision
    conf << activeInteractions[i]->n << CommBox().sep;
    std::fesetround(prev_round);
    conf << activeInteractions[i]->dn << CommBox().sep << activeInteractions[i]->pos << CommBox().sep
         << activeInteractions[i]->vel << CommBox().sep << activeInteractions[i]->fn << CommBox().sep
         << activeInteractions[i]->ft << CommBox().sep << activeInteractions[i]->mom << CommBox().sep
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
      if (ParamsInInterfaces == 1) {
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
}

void Rockable::readLawData(std::istream& is, size_t id) {
  size_t g1, g2;
  double value;
  is >> g1 >> g2 >> value;
  dataTable.set(id, g1, g2, value);
}

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


/**
    @brief Load a configuration-file named 'conf<i>'
*/
void Rockable::loadConf(int i) {
  char fname[256];
  sprintf(fname, "conf%d", i);
  loadConf(fname);
}

/**
    @brief Load a configuration file named name
    @param[in]  name  The name of the conf-file
*/
void Rockable::loadConf(const char* name) {
  std::ifstream conf(name);
  if (!conf.is_open()) {
    std::cerr << msg::warn() << "@Rockable::loadConf, Cannot read " << name << msg::normal() << std::endl;
    return;
  }

  // Check header
  std::string prog;
  conf >> prog;
  if (prog != "Rockable") {
    std::cerr << msg::warn() << "@Rockable::loadConf, This is not a file for the code Rockable!" << msg::normal()
              << std::endl;
  }
  std::string date;
  conf >> date;
  if (date != CONF_VERSION_DATE) {
    std::cerr << msg::warn() << "@Rockable::loadConf, The version-date should be '" << CONF_VERSION_DATE << "'!\n"
              << "                     in most cases, this should not be a problem." << msg::normal() << std::endl;
  }

  kwParser parser;

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
  parser.kwMap["ParamsInInterfaces"] = __GET__(conf, ParamsInInterfaces);
  parser.kwMap["dynamicUpdateNL"] = __GET__(conf, dynamicUpdateNL);
  parser.kwMap["dispUpdateNL"] = __GET__(conf, dispUpdateNL);
  parser.kwMap["angleUpdateNL"] = __DO__(conf) {
    conf >> angleUpdateNL;
    angleUpdateNL *= Mth::deg2rad;
  };
  parser.kwMap["numericalDampingCoeff"] = __GET__(conf, numericalDampingCoeff);
  parser.kwMap["VelocityBarrier"] = __GET__(conf, VelocityBarrier);
  parser.kwMap["AngularVelocityBarrier"] = __GET__(conf, AngularVelocityBarrier);
  parser.kwMap["VelocityBarrierExponent"] = __GET__(conf, VelocityBarrierExponent);
  parser.kwMap["AngularVelocityBarrierExponent"] = __GET__(conf, AngularVelocityBarrierExponent);

  parser.kwMap["Tempo"] = [this](std::istream& conf) {
    std::string kw;
    conf >> kw;
    if (kw == "numericalDampingCoeff") {
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
    if (lawName == "Default") {
      forceLawPtr = std::bind(&Rockable::forceLawDefault, this, std::placeholders::_1);
      optionNames["forceLaw"] = "Default";
    } else if (lawName == "Avalanches") {
      forceLawPtr = std::bind(&Rockable::forceLawAvalanches, this, std::placeholders::_1);
      optionNames["forceLaw"] = "Avalanches";
    } else if (lawName == "StickedLinks") {
      forceLawPtr = std::bind(&Rockable::forceLawStickedLinks, this, std::placeholders::_1);
      optionNames["forceLaw"] = "StickedLinks";
    } else {
      std::cout << msg::warn() << "forceLaw " << lawName << " is unknown" << msg::normal() << std::endl;
      std::cout << "Option remains: forceLaw = " << optionNames["forceLaw"] << std::endl;
    }
  };
  parser.kwMap["AddOrRemoveInteractions"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    if (Name == "bruteForce") {
      AddOrRemoveInteractions = std::bind(&Rockable::AddOrRemoveInteractions_bruteForce, this, std::placeholders::_1,
                                          std::placeholders::_2, std::placeholders::_3);
      optionNames["AddOrRemoveInteractions"] = "bruteForce";
    } else if (Name == "OBBtree") {
      AddOrRemoveInteractions = std::bind(&Rockable::AddOrRemoveInteractions_OBBtree, this, std::placeholders::_1,
                                          std::placeholders::_2, std::placeholders::_3);
      optionNames["AddOrRemoveInteractions"] = "OBBtree";
    } else {
      std::cout << msg::warn() << "AddOrRemoveInteractions " << Name << " is unknown" << msg::normal() << std::endl;
      std::cout << "Option remains: AddOrRemoveInteractions = " << optionNames["AddOrRemoveInteractions"] << std::endl;
    }
  };
  parser.kwMap["UpdateNL"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    if (Name == "bruteForce") {
      UpdateNL = std::bind(&Rockable::UpdateNL_bruteForce, this);
      optionNames["UpdateNL"] = "bruteForce";
    } else if (Name == "linkCells") {
      UpdateNL = std::bind(&Rockable::UpdateNL_linkCells, this);
      optionNames["UpdateNL"] = "linkCells";
    } else {
      std::cout << msg::warn() << "UpdateNL " << Name << " is unknown" << msg::normal() << std::endl;
      std::cout << "Option remains: UpdateNL = " << optionNames["UpdateNL"] << std::endl;
    }
  };
  parser.kwMap["Integrator"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    if (Name == "velocityVerlet") {
      IntegrationStep = std::bind(&Rockable::velocityVerletStep, this);
      optionNames["Integrator"] = "velocityVerlet";
    } else if (Name == "Euler") {
      IntegrationStep = std::bind(&Rockable::EulerStep, this);
      optionNames["Integrator"] = "Euler";
    } else if (Name == "Beeman") {
      IntegrationStep = std::bind(&Rockable::BeemanStep, this);
      optionNames["Integrator"] = "Beeman";
    } else if (Name == "RungeKutta4") {
      IntegrationStep = std::bind(&Rockable::RungeKutta4Step, this);
      optionNames["Integrator"] = "RungeKutta4";
    } else {
      std::cout << msg::warn() << "Integrator " << Name << " is unknown" << msg::normal() << std::endl;
      std::cout << "Option remains: Integrator = " << optionNames["Integrator"] << std::endl;
    }
  };
  parser.kwMap["cellMinSizes"] = __GET__(conf, cellMinSizes);
  parser.kwMap["boxForLinkCellsOpt"] = __GET__(conf, boxForLinkCellsOpt);
  parser.kwMap["ContactPartnership"] = __DO__(conf) {
    std::string Name;
    conf >> Name;
    ctcPartnership.setModel(Name);
  };

  parser.kwMap["knContact"] = __DO__(conf) { readLawData(conf, idKnContact); };
  parser.kwMap["en2Contact"] = __DO__(conf) { readLawData(conf, idEn2Contact); };
  parser.kwMap["ktContact"] = __DO__(conf) { readLawData(conf, idKtContact); };
  parser.kwMap["muContact"] = __DO__(conf) { readLawData(conf, idMuContact); };
  parser.kwMap["krContact"] = __DO__(conf) { readLawData(conf, idKrContact); };
  parser.kwMap["murContact"] = __DO__(conf) { readLawData(conf, idMurContact); };

  parser.kwMap["knInnerBond"] = __DO__(conf) { readLawData(conf, idKnInnerBond); };
  parser.kwMap["ktInnerBond"] = __DO__(conf) { readLawData(conf, idKtInnerBond); };
  parser.kwMap["fn0InnerBond"] = __DO__(conf) { readLawData(conf, idFn0InnerBond); };
  parser.kwMap["ft0InnerBond"] = __DO__(conf) { readLawData(conf, idFt0InnerBond); };
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
    std::string wantedLib;
    conf >> wantedLib;
    std::cout << "wantedLib is " << wantedLib << std::endl;
    if (wantedLib != shapeFile) {  // it means that the library is not already loaded
      shapeFile = wantedLib;
      loadShapes(shapeFile.c_str());
    }
  };
  parser.kwMap["separator"] = __DO__(conf) {
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
    std::cout << "Number of bodies: " << nb << std::endl;
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
      if (shpName[0] == '!') {  // to 'remove' a line of particle (the number of
                                // particles is decreased but you don't need to change
                                // the number of particles in the input file)
        std::string trash;
        getline(conf, trash);
        --i;
        --nb;
        continue;
      }
      conf >> P.group >> P.cluster >> P.homothety >> P.pos >> P.vel >> P.acc >> P.Q >> P.vrot >> P.arot;

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
    std::cout << "Number of interactions: " << nb << std::endl;

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
    std::cout << "Number of interfaces: " << nbInterf << std::endl;

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

      if (ParamsInInterfaces == 1) {
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
          std::cerr << msg::warn() << "Cannot find interaction: " << ItoFind.i << " " << ItoFind.j << " "
                    << ItoFind.type << " " << ItoFind.isub << " " << ItoFind.jsub << msg::normal() << std::endl;
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

    std::cerr << msg::warn() << "The DataExtractor named " << ExtractorName
              << " is defined in the input file! It will not be saved in conf-files\n"
              << "A better solution is to put them in a file named 'dataExtractors.txt'" << msg::normal() << std::endl;

    DataExtractor* DE = Factory<DataExtractor>::Instance()->Create(ExtractorName);
    if (DE != nullptr) {
      DE->plug(this);
      DE->read(conf);
      dataExtractors.push_back(DE);
    } else {
      std::cerr << msg::warn() << "The DataExtractor named " << ExtractorName << " is unknown!" << msg::normal()
                << std::endl;
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
      if (verbose == 1) {
        std::cout << "The BodyForce named " << BodyForceName << " has been activated\n";
      }
    } else {
      std::cerr << msg::warn() << "The BodyForce named " << BodyForceName << " is unknown!" << msg::normal()
                << std::endl;
    }
  };

  // ======== FROM HERE, THESE ARE PREPRO COMMANDS (not saved in further conf-files)
  //          They are generally put at the end of the input file

  parser.kwMap["stickVerticesInClusters"] = __DO__(conf) {
    double epsilonDist;
    conf >> epsilonDist;
    stickVerticesInClusters(epsilonDist);
  };
  parser.kwMap["stickClusters"] = __DO__(conf) {
    double epsilonDist;
    conf >> epsilonDist;
    stickClusters(epsilonDist);
  };
  parser.kwMap["copyParamsToInterfaces"] = __DO__(conf) {
    std::string isInnerStr;
    conf >> isInnerStr;
    copyParamsToInterfaces(isInnerStr);
  };
  parser.kwMap["setStiffnessRatioInterfaces"] = __DO__(conf) {
    double ratio;
    conf >> ratio;
    setStiffnessRatioInterfaces(ratio);
  };
  parser.kwMap["setVariableStickParams"] = __DO__(conf) {
    std::string paramName;
    std::string isInnerStr;
    double lambda, m;
    int timeSeeded;
    conf >> paramName >> isInnerStr >> lambda >> m >> timeSeeded;
    setVariableStickParams(paramName, isInnerStr, lambda, m, (bool)timeSeeded);
  };
  parser.kwMap["randomlyOrientedVelocities"] = __DO__(conf) {
    double vel;
    conf >> vel;
    randomlyOrientedVelocities(vel);
  };
  parser.kwMap["randomlyOrientedVelocitiesClusters"] = __DO__(conf) {
    double vel;
    int opt;
    conf >> vel >> opt;
    randomlyOrientedVelocitiesClusters(vel, opt);
  };
  parser.kwMap["homothetyRange"] = __DO__(conf) {
    size_t ifirst, ilast;
    double hmin, hmax;
    int timeSeeded;
    conf >> ifirst >> ilast >> hmin >> hmax >> timeSeeded;
    homothetyRange(ifirst, ilast, hmin, hmax, (bool)timeSeeded);
  };
  parser.kwMap["particlesClonage"] = __DO__(conf) {
    size_t ifirst, ilast;
    vec3r translation;
    conf >> ifirst >> ilast >> translation;
    particlesClonage(ifirst, ilast, translation);
  };

  // This single line actually parses the file
  parser.parse(conf);

  // This is a kind of fake time-step (the time is not incremented)
  // mainly used to compute resultant forces and moments on the particles.
  // It will also populate the vector 'activeInteractions'.
  if (interactiveMode == false) accelerations();
}

/**
    Read the DataExtractor commands in a file named 'dataExractors.txt'.
    If the file is not found, or if the interactiveMode is true, it simply
    exit without making anything.
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
        std::cerr << msg::warn() << "The DataExtractor named " << ExtractorName << " is unknown!" << msg::normal()
                  << std::endl;
      }
    }
    is >> ExtractorName;
  }
}

/**
    @brief Load all shapes defined in the file 'fileName'.
           If the library has already been loaded (ie. the file name is the same as
           the one in the previously read conf-file), then it is not re-read.
*/
void Rockable::loadShapes(const char* fileName) {
  // If a library file is in the running folder, so it is preferably used
  std::string ModFileName(fileName);
  std::string LocalFileName = fileTool::GetFileName(ModFileName) + "." + fileTool::GetFileExt(ModFileName);
  if (fileTool::fileExists(LocalFileName.c_str())) {
    ModFileName = LocalFileName;
  }

  if (!fileTool::fileExists(ModFileName.c_str())) {
    std::cerr << msg::warn() << "@Rockable::loadShapes, shape library named '" << ModFileName << "' has not been found."
              << msg::normal() << std::endl
              << std::endl;
    return;
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

  std::cout << "Number of shapes found in the library file " << ModFileName << ": " << Shapes.size() << std::endl;

  for (size_t s = 0; s < Shapes.size(); ++s) {
    Shapes[s].buildOBBtree();
    shapeId[Shapes[s].name] = s;
  }
}

// ==================================================================================================================
//  ADD OR REMOVE A SINGLE INTERACTION
// ==================================================================================================================

/**
    @brief This is the brute-force O(n^2) version of the algorithm.
           It means that the proximity of all sub-elements of i is tested with all
           sub-elements of j
*/
int Rockable::AddOrRemoveInteractions_bruteForce(size_t i, size_t j, double dmax) {

  double Damp = 0.0;
  int nbAdd = 0;

  // A helper lambda function
  auto addOrRemoveSingleInteraction = [&](size_t i, size_t j, size_t isub, size_t type, size_t nbj,
                                          std::function<bool(Particle&, Particle&, size_t, size_t, double)> func) {
    Interaction to_find;
    to_find.i = i;
    to_find.j = j;
    to_find.type = type;
    to_find.isub = isub;
    for (size_t jsub = 0; jsub < nbj; ++jsub) {
      to_find.jsub = jsub;
      auto exist_it = (Interactions[i]).find(to_find);
      bool NEW = (exist_it == Interactions[i].end());
      bool NEAR = func(Particles[i], Particles[j], isub, jsub, dmax);
      if (NEAR && NEW) {
        Interactions[i].insert(Interaction(i, j, type, isub, jsub, Damp));
        ++nbAdd;
      } else if (!NEAR && !NEW) {
        Interactions[i].erase(exist_it);
        --nbAdd;
      }
    }
  };

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
    obbj.enlarge(dmax);
    if (!(obbj.intersect(subObbi))) continue;

    // vertex-vertex i->j
    addOrRemoveSingleInteraction(i, j, isub, vvType, nvj, Particle::VertexIsNearVertex);

    // vertex-edge i->j
    addOrRemoveSingleInteraction(i, j, isub, veType, nej, Particle::VertexIsNearEdge);

    // vertex-face i->j
    addOrRemoveSingleInteraction(i, j, isub, vfType, nfj, Particle::VertexIsNearFace);

  }  // end for isub (i->j)

  for (size_t jsub = 0; jsub < nvj; ++jsub) {

    OBB subObbj;
    subObbj.center = Particles[j].GlobVertex(jsub);
    subObbj.enlarge(Particles[j].MinskowskiRadius() + dmax);
    OBB obbi = Particles[i].obb;
    obbi.enlarge(dmax);
    if (!(obbi.intersect(subObbj))) continue;

    // vertex-edge j->i
    addOrRemoveSingleInteraction(j, i, jsub, veType, nei, Particle::VertexIsNearEdge);

    // vertex-face j->i
    addOrRemoveSingleInteraction(j, i, jsub, vfType, nfi, Particle::VertexIsNearFace);

  }  // end for jsub (j->i)

  // edge-edge i->j
  for (size_t isub = 0; isub < nei; ++isub) {
    addOrRemoveSingleInteraction(i, j, isub, eeType, nej, Particle::EdgeIsNearEdge);
  }

  return nbAdd;
}

/**
    This version should be best adapted when the particles have a lot of sub-elements
    (a terrain for example). When the number of sub-elements is relatively small,
    it seems not to dramatically slow down the computation.
    Remark: obbi and obbj need to be already placed BEFORE calling this method
*/
int Rockable::AddOrRemoveInteractions_OBBtree(size_t i, size_t j, double dmax) {
  static Interaction to_find;
  static std::set<Interaction>::iterator exist_it;

  double Damp = 0.0;
  int nbAdd = 0;

  // A helper lambda function
  auto addOrRemoveSingleInteraction = [&](size_t i, size_t j, size_t isub, size_t type, size_t jsub,
                                          std::function<bool(Particle&, Particle&, size_t, size_t, double)> func) {
    to_find.i = i;
    to_find.j = j;
    to_find.type = type;
    to_find.isub = isub;
    to_find.jsub = jsub;
    exist_it = (Interactions[i]).find(to_find);
    bool NEW = (exist_it == Interactions[i].end());
    bool NEAR = func(Particles[i], Particles[j], isub, jsub, dmax);
    if (NEAR && NEW) {
      Interactions[i].insert(Interaction(i, j, type, isub, jsub, Damp));
      ++nbAdd;
    } else if (!NEAR && !NEW) {
      Interactions[i].erase(exist_it);
      --nbAdd;
    }
  };

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
      if (j_nbPoints == 1) { // vertex (i, isub) -> vertex (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, vvType, jsub, Particle::VertexIsNearVertex);
      } else if (j_nbPoints == 2) { // vertex (i, isub) -> edge (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, veType, jsub, Particle::VertexIsNearEdge);
      } else if (j_nbPoints == 3) { // vertex (i, isub) -> face (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, vfType, jsub, Particle::VertexIsNearFace);
      }
    } else if (i_nbPoints == 2) { 
      if (j_nbPoints == 1) { // vertex (j, jsub) -> edge (i, isub)
        addOrRemoveSingleInteraction(j, i, jsub, veType, isub, Particle::VertexIsNearEdge);
      } else if (j_nbPoints == 2) { // vertex (i, isub) -> edge (j, jsub)
        addOrRemoveSingleInteraction(i, j, isub, eeType, jsub, Particle::EdgeIsNearEdge);
      }
    } else if (i_nbPoints == 3) {
      addOrRemoveSingleInteraction(j, i, jsub, vfType, isub, Particle::VertexIsNearFace); 
      // vertex (j, jsub) -> face (i, isub)
    }

  }  // end loop over intersections

  return nbAdd;
}

// ==================================================================================================================
//  UPDATING THE NEIGHBOR LIST
// ==================================================================================================================

void Rockable::dynamicCheckUpdateNL() {

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

/**
    @brief  The most basic algorithm for building a neighbor list.
            That is the O(N^2) complexity
*/
void Rockable::UpdateNL_bruteForce() {
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].updateObb();
  }

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < Particles.size(); ++i) {
#ifndef BREAK_ONCE
    BreakableInterface BI_to_find;
    BI_to_find.i = i;
#endif

    OBB obbi = Particles[i].obb;
    obbi.enlarge(0.5 * DVerlet);

    size_t jnext = i + 1;
    // Prevent interactions between driven bodies
    if (jnext < nDriven) jnext = nDriven;

    for (size_t j = jnext; j < Particles.size(); j++) {

#ifndef BREAK_ONCE
      // if interface (i,j) is still 'active' (ie. no bond has been broken),
      // then NO contact will be possible between bodies i and j.
      BI_to_find.j = j;
      std::set<BreakableInterface>::iterator BI_it = (Interfaces[i]).find(BI_to_find);
      if (BI_it != Interfaces[i].end()) continue;  // Continue the loop (next j) if an interface is found
#endif

      OBB obbj = Particles[j].obb;
      obbj.enlarge(0.5 * DVerlet);

      // Check intersection
      if (obbi.intersect(obbj)) {
        AddOrRemoveInteractions(i, j, dVerlet);
      }
    }
  }
}

/**
    @brief  A broad-phase collision detection that should have a complexity of nearly O(N).
*/
void Rockable::UpdateNL_linkCells() {
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
#ifndef BREAK_ONCE
            BreakableInterface BI_to_find;
            BI_to_find.i = i;
#endif
            OBB obbi = Particles[i].obb;
            obbi.enlarge(0.5 * DVerlet);

            for (size_t jcv = 0; jcv < Cv->bodies.size(); ++jcv) {
              size_t j = Cv->bodies[jcv];
              if (j < nDriven && i < nDriven) continue;
              if (j <= i) continue;

#ifndef BREAK_ONCE
              // if interface (i,j) is still 'active' (ie. no bond has been broken),
              // then NO contact will be possible between bodies i and j.
              BI_to_find.j = j;
              std::set<BreakableInterface>::iterator BI_it = (Interfaces[i]).find(BI_to_find);
              if (BI_it != Interfaces[i].end()) continue;  // Continue the loop if an interface is found
#endif

              OBB obbj = Particles[j].obb;
              obbj.enlarge(0.5 * DVerlet);

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
#ifndef BREAK_ONCE
          BreakableInterface BI_to_find;
          BI_to_find.i = i;
#endif
          OBB obbi = Particles[i].obb;
          obbi.enlarge(0.5 * DVerlet);

          for (size_t jcv = 0; jcv < Cv->bodies.size(); ++jcv) {
            size_t j = Cv->bodies[jcv];

            if (j < nDriven && i < nDriven) continue;
            if (j <= i) continue;

#ifndef BREAK_ONCE
            // if interface (i,j) is still 'active' (ie. no bond has been broken),
            // then NO contact will be possible between bodies i and j.
            BI_to_find.j = j;
            std::set<BreakableInterface>::iterator BI_it = (Interfaces[i]).find(BI_to_find);
            if (BI_it != Interfaces[i].end()) continue;  // Continue the loop if an interface is found
#endif

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
}

// ======================
// FORCE LAWS
// ======================

/**
    @brief Since the group number can be zoned in a particle,
           this method is meant to get these right numbers
*/
void Rockable::getInteractingGroups(Interaction& I, int& g1, int& g2) {
  g1 = Particles[I.i].group;
  g2 = Particles[I.j].group;
}

/**
    @brief   Maybe the simplest law to use.
             It includes Linear normal repulsion, normal viscosity, Coulomb Friction, and moment resistance
    @return  Return true if the interaction is active (ie. with a non-zero force)
*/
bool Rockable::forceLawDefault(Interaction& I) {
  if (I.dn > 0.0) {
    I.fn = 0.0;
    I.ft.reset();
    I.mom.reset();
    return false;
  }

  int g1 = Particles[I.i].group;
  int g2 = Particles[I.j].group;
  double kn = dataTable.get(idKnContact, g1, g2);
  double kt = dataTable.get(idKtContact, g1, g2);
  double mu = dataTable.get(idMuContact, g1, g2);
  double kr = dataTable.get(idKrContact, g1, g2);
  double mur = dataTable.get(idMurContact, g1, g2);
  double damp = I.damp;

  if (ctcPartnership.getWeight != nullptr) {
    double w = ctcPartnership.getWeight(I);
    kn *= w;
    kt *= w;
    kr *= w;
    damp *= sqrt(w);
  }

  // === Normal force (elatic contact + viscous damping)
  double vn = I.vel * I.n;
  double fne = -kn * I.dn;
  double fnv = damp * vn;
  I.fn = fne + fnv;

  // === Tangential force (friction)
  vec3r vt = I.vel - vn * I.n;
#ifdef FT_CORR
  vec3r ft_corr = I.ft;
  ft_corr -= cross(ft_corr, cross(I.prev_n, I.n));
  ft_corr -= cross(ft_corr, (dt_2 * (Particles[I.i].vrot + Particles[I.j].vrot) * I.n) * I.n);
  I.ft = ft_corr + kt * (vt * dt);
#else
  I.ft += kt * (vt * dt);
#endif

  double threshold = fabs(mu * fne);
  double ft_square = I.ft * I.ft;
  if (ft_square > 0.0 && ft_square >= threshold * threshold) I.ft *= threshold / sqrt(ft_square);
  // Remark: in fact, the test (ft_square > 0.0) means that ft_square is not
  // zero, because ft_square >= 0 by definition.

  // === Resistant moment
  if (kr > 0.0) {
    I.mom += kr * (Particles[I.j].vrot - Particles[I.i].vrot) * dt;
    double threshold_mom = fabs(mur * I.fn);  // in this case mur is a *length*
    double mom_square = I.mom * I.mom;
    if (mom_square > 0.0 && mom_square >= threshold_mom * threshold_mom) I.mom *= threshold_mom / sqrt(mom_square);
  }

  return true;
}

/**
    @brief   Force-law used for rock avalanches at Laboratoire 3SR
    @return  Return true if the interaction is active (ie. with a non-zero force)
*/
bool Rockable::forceLawAvalanches(Interaction& I) {
  if (I.dn > 0.0) {
    I.fn = 0.0;
    I.ft.reset();
    I.mom.reset();
    return false;
  }

  // getInteractingGroups(I, g1, g2);
  int g1 = Particles[I.i].group;
  int g2 = Particles[I.j].group;
  double kn = dataTable.get(idKnContact, g1, g2);
  double en2 = dataTable.get(idEn2Contact, g1, g2);
  double kt = dataTable.get(idKtContact, g1, g2);
  double mu = dataTable.get(idMuContact, g1, g2);
  double kr = dataTable.get(idKrContact, g1, g2);
  double mur = dataTable.get(idMurContact, g1, g2);

  if (ctcPartnership.getWeight != nullptr) {
    double w = ctcPartnership.getWeight(I);
    kn *= w;
    kt *= w;
    kr *= w;
  }

  // === Normal force
  if (I.prev_dn > 0.0) I.prev_dn = 0.0;
  double delta_Dn = I.dn - I.prev_dn;
  if (delta_Dn > 0.0)
    I.fn = -kn * en2 * I.dn;  // Unloading
  else if (delta_Dn < 0.0)
    I.fn += -kn * delta_Dn;  // Loading
  if (I.fn < 0.0) I.fn = 0.0;

  // === Tangential force (friction)
  vec3r vt = I.vel - (I.vel * I.n) * I.n;
#ifdef FT_CORR
  vec3r ft_corr = I.ft;
  ft_corr -= cross(ft_corr, cross(I.prev_n, I.n));
  ft_corr -= cross(ft_corr, (dt_2 * (Particles[I.i].vrot + Particles[I.j].vrot) * I.n) * I.n);
  I.ft = ft_corr + kt * (vt * dt);
#else
  I.ft += kt * (dt * vt);
#endif
  double threshold_ft = fabs(mu * I.fn);  // even without fabs the value should be positive
  double ft_square = I.ft * I.ft;
  if (ft_square > 0.0 && ft_square > threshold_ft * threshold_ft) I.ft *= threshold_ft / sqrt(ft_square);
  // Remark: in fact, the test (ft * ft > 0.0) means that ft_square is NOT null,
  // because ft * ft >= 0 by definition.

  // === Resistant moment
  I.mom += kr * (Particles[I.j].vrot - Particles[I.i].vrot) * dt;
  vec3r branch;
  if (I.i < nDriven) {  // j is the free body (a rock block)
    branch = I.pos - Particles[I.j].pos;
    double r = (branch * Particles[I.j].vrot) / (Particles[I.j].vrot * Particles[I.j].vrot);
    branch -= r * Particles[I.j].vrot;
  } else {  // i is the free body (a rock block)
    branch = I.pos - Particles[I.i].pos;
    double r = (branch * Particles[I.i].vrot) / (Particles[I.i].vrot * Particles[I.i].vrot);
    branch -= r * Particles[I.i].vrot;
  }
  double threshold_mom = fabs(mur * norm(branch) * I.fn);  // even without fabs, the value should
                                                           // be positive
  double mom_square = I.mom * I.mom;
  if (mom_square > 0.0 && mom_square > threshold_mom * threshold_mom) I.mom *= threshold_mom / sqrt(mom_square);

  return true;
}

/**
   This is the force-law initiated by ANDRA's study (PhD of Marta Stasiak)
   Bodies can be glued and when the glue is 'broken', it is irreversibly
   switched to frictional contact
*/
bool Rockable::forceLawStickedLinks(Interaction& I) {
  if (I.stick != nullptr) {  // =========== Cohesive bond
    double kn = 0.0, kt = 0.0, kr = 0.0;
    double fn0 = 1.0, ft0 = 1.0, mom0 = 1.0, power = 1.0;
    double dn0 = I.stick->dn0;
    bool isInner = false;

    // First we need to get the parameters
    if (ParamsInInterfaces == 1) {
      isInner = I.stick->isInner;
      kn = I.stick->kn;
      kt = I.stick->kt;
      kr = I.stick->kr;
      fn0 = I.stick->fn0;
      ft0 = I.stick->ft0;
      mom0 = I.stick->mom0;
      power = I.stick->power;
    } else {  // set parameters according to the group-numbers
      int g1 = Particles[I.i].group;
      int g2 = Particles[I.j].group;

      if (Particles[I.i].cluster == Particles[I.j].cluster) {  // Inner
        isInner = true;
        kn = dataTable.get(idKnInnerBond, g1, g2);
        kt = dataTable.get(idKtInnerBond, g1, g2);
        fn0 = dataTable.get(idFn0InnerBond, g1, g2);
        ft0 = dataTable.get(idFt0InnerBond, g1, g2);
        power = dataTable.get(idPowInnerBond, g1, g2);
        dn0 = 0.0;
      } else {  // Outer
        isInner = false;
        kn = dataTable.get(idKnOuterBond, g1, g2);
        kt = dataTable.get(idKtOuterBond, g1, g2);
        kr = dataTable.get(idKrOuterBond, g1, g2);
        fn0 = dataTable.get(idFn0OuterBond, g1, g2);
        ft0 = dataTable.get(idFt0OuterBond, g1, g2);
        mom0 = dataTable.get(idMom0OuterBond, g1, g2);
        power = dataTable.get(idPowOuterBond, g1, g2);
      }
    }

    // === Normal force (elastic contact + viscous damping)
    double vn = I.vel * I.n;
    double fne = -kn * (I.dn - dn0);
    double fnv = I.damp * vn;
    I.fn = fne + fnv;

    // === Tangential force (friction)
    vec3r vt = (I.vel - (vn * I.n));
#ifdef FT_CORR
    vec3r ft_corr = I.ft;
    ft_corr -= cross(ft_corr, cross(I.prev_n, I.n));
    ft_corr -= cross(ft_corr, (dt_2 * (Particles[I.i].vrot + Particles[I.j].vrot) * I.n) * I.n);
    I.ft = ft_corr + kt * (vt * dt);
#else
    I.ft += kt * (vt * dt);
#endif

    // Tangential viscosity ========
    // This term should be added only on the elastic part of ft
    // So it is somehow wrong because the viscosity is cumulated... Be carreful!
    vec3r ftv = I.damp * vt;
    I.ft += ftv;
    // =============================

    // === Rupture criterion (and resistant moment for outer bonds)
    double f;  // it defines the yield surface
    if (isInner == true) {
      f = pow(norm(I.ft) / ft0, power) - I.fn / fn0 - 1.0;
    } else {
      I.mom += kr * (Particles[I.j].vrot - Particles[I.i].vrot) * dt;
      f = pow(norm(I.ft) / ft0, power) + pow(norm(I.mom) / mom0, power) - I.fn / fn0 - 1.0;
    }

    if (f > 0.0) {
#ifdef BREAK_ONCE
#pragma message "It is not recommanded to use BREAK_ONCE"
      I.fn = 0.0;
      I.ft.reset();
      I.mom.reset();
      I.stick = nullptr;
      needUpdate = true;
#else
      // All the bonds (Interactions) of the interface are broken.
      // The Interactions pointer inserted in this std::set
      // will be 'broken' just after all the forces are computed
      interfacesToBreak.insert(I.stick);
      // The breakage is postponed to avoid assymetry in the tangential forces
      // that are computed incrementally
#endif
    }
  } else {  // =============== Contact
    if (I.dn > 0.0) {
      I.fn = 0.0;
      I.ft.reset();
      I.mom.reset();
      return false;
    }

    int g1 = Particles[I.i].group;
    int g2 = Particles[I.j].group;
    double kn = dataTable.get(idKnContact, g1, g2);
    double kt = dataTable.get(idKtContact, g1, g2);
    double mu = dataTable.get(idMuContact, g1, g2);
    double kr = dataTable.get(idKrContact, g1, g2);
    double mur = dataTable.get(idMurContact, g1, g2);
    double damp = I.damp;

    if (ctcPartnership.getWeight != nullptr) {
      double w = ctcPartnership.getWeight(I);
      kn *= w;
      kt *= w;
      kr *= w;
      damp *= sqrt(w);
    }

    // === Normal force (elatic contact + viscous damping)
    double vn = I.vel * I.n;
    double fne = -kn * I.dn;
    double fnv = damp * vn;
    I.fn = fne + fnv;
    if (I.fn < 0.0) I.fn = 0.0;  // Because viscous damping can make the normal force negative

    // === Tangential force (friction)
    vec3r vt = I.vel - vn * I.n;
#ifdef FT_CORR
    vec3r ft_corr = I.ft;
    ft_corr -= cross(ft_corr, cross(I.prev_n, I.n));
    ft_corr -= cross(ft_corr, (dt_2 * (Particles[I.i].vrot + Particles[I.j].vrot) * I.n) * I.n);
    I.ft = ft_corr + kt * (vt * dt);
#else
    I.ft += kt * (vt * dt);
#endif

    double threshold = fabs(mu * I.fn);
    double ft_square = I.ft * I.ft;
    if (ft_square > 0.0 && ft_square >= threshold * threshold) I.ft *= threshold / sqrt(ft_square);
    // Remark: in fact, the test (ft_square > 0.0) means that ft_square is not
    // zero, because ft_square >= 0 by definition.

    // === Resistant moment
    I.mom += kr * (Particles[I.j].vrot - Particles[I.i].vrot) * dt;
    double threshold_mom = fabs(mur * I.fn);  // in this case mur is a *length*
    double mom_square = I.mom * I.mom;
    if (mom_square > 0.0 && mom_square >= threshold_mom * threshold_mom) I.mom *= threshold_mom / sqrt(mom_square);
  }
  return true;
}

// ==============================================================================================================
//  CORE METHODS OF THE DEM ALGORITHM
// ==============================================================================================================

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

void Rockable::velocityControlledDrive() {
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
      default:
        break;
    }
  }
}

/**
    @brief Makes ONE step with the Euler scheme.
*/
void Rockable::EulerStep() {
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
      default:
        break;
    }
  }

  // Bodies that are free to move
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].pos += dt * Particles[i].vel;
    Particles[i].vel += dt * Particles[i].acc;

    // Rotation: Q <- Q + (dQ/dt)*dt
    // It reads like this with quaternions
    Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
    Particles[i].Q.normalize();

    Particles[i].vrot += dt * Particles[i].arot;
  }

  accelerations();
}

/**
    @brief Makes ONE step with the velocity-Verlet scheme.
*/
void Rockable::velocityVerletStep() {
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
      default:
        break;
    }
  }

  // Bodies that are free to move
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].pos += dt * Particles[i].vel + dt2_2 * Particles[i].acc;
    Particles[i].vel += dt_2 * Particles[i].acc;

    // Rotation: Q(k+1) = Q(k) + dQ(k) * dt + ddQ(k) * dt2/2
    // It reads like this with quaternions
    Particles[i].Q += ((Particles[i].Q.dot(Particles[i].vrot)) *= dt);
#ifdef QUAT_ACC
    Particles[i].Q += ((Particles[i].Q.ddot(Particles[i].vrot, Particles[i].arot)) *= dt2_2);
#endif
    Particles[i].Q.normalize();

    Particles[i].vrot += dt_2 * Particles[i].arot;
  }

  accelerations();

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
      default:
        break;
    }
  }

  // Free bodies
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].vel += dt_2 * Particles[i].acc;
    Particles[i].vrot += dt_2 * Particles[i].arot;
  }
}

/**
    @brief Makes ONE step with the Beeman scheme.
*/
void Rockable::BeemanStep() {
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
    @brief Makes ONE step with the Runge-Kutta-Nystrom (4th order) scheme.
*/
void Rockable::RungeKutta4Step() {
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

/**
    @brief The integration LOOP
*/
void Rockable::integrate() {
  if (interactiveMode == true) {
    std::cout << "It is not possible to invoke Rockable::integrate if "
                 "interactiveMode is true"
              << std::endl;
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

  initIntegrator();

  // Save the current configuration
  saveConf(iconf);
  std::cout << "#  iconf = " << iconf << ", Time = " << t << std::endl;

  PerfTimer ptimer;
  size_t step = (size_t)(t / dt);
  timeInUpdateNL = 0.0;
  timeInForceComputation = 0.0;
  while (t < tmax) {

    if (interConfC >= interConf - dt_2) {
      iconf++;

      std::cout << "+-----------------------------------------------------------------------------------\n";
      std::cout << "|  #iconf = " << iconf << ", Time = " << t << '\n';
      double elapsedTime = ptimer.getIntermediateElapsedTimeSeconds();
      double NLPercent = std::round(1000.0 * timeInUpdateNL / elapsedTime) /
                         10.0;  // It's a way to obtain a 1/10 precision (in percents)
      double ForcePercent = std::round(1000.0 * timeInForceComputation / elapsedTime) / 10.0;
      std::cout << "|   Elapsed time since last configuration: " << msg::HumanReadableSeconds(elapsedTime) << '\n';
      std::cout << "|   " << __printNamedCell(10, "%NL:", 13, NLPercent)
                << __printNamedCell(13, "%Forces:", 13, ForcePercent)
                << __printNamedCell(13, "%Rest:", 13, 100.0 - NLPercent - ForcePercent) << '\n';

      double efficiency = interConf / elapsedTime;
      perfFile << t << " " << efficiency << " " << NLPercent << " " << ForcePercent << '\n' << std::flush;

      timeInUpdateNL = 0.0;
      timeInForceComputation = 0.0;

      double Fmax, F_fnmax, Fmean, Fstddev;
      getResultantQuickStats(Fmax, F_fnmax, Fmean, Fstddev, nDriven);
      std::cout << "|   Resultant forces on particles:\n|   " << __printNamedCell(13, "Fmax:", 13, Fmax)
                << __printNamedCell(13, "Fmean:", 13, Fmean) << __printNamedCell(13, "Fstddev:", 13, Fstddev) << '\n';
      double fnMin, fnMax, fnMean, fnStddev;
      getInteractionQuickStats(fnMin, fnMax, fnMean, fnStddev);
      std::cout << "|   Interaction forces:\n|   " << __printNamedCell(13, "fnMin:", 13, fnMin)
                << __printNamedCell(13, "fnMax:", 13, fnMax) << "\n|   " << __printNamedCell(13, "fnMean:", 13, fnMean)
                << __printNamedCell(13, "fnStddev:", 13, fnStddev) << '\n';

      if (fnMax != 0.0 && fnMean != 0.0) {
        double Fmax_fnMean = fabs(Fmax / fnMean);
        std::cout << "|   Static balance:\n|   " << __printNamedCell(13, "max(F/fnMax):", 13, F_fnmax)
                  << __printNamedCell(13, "Fmax/fnMean:", 13, Fmax_fnMean) << '\n';
        staticBalanceFile << t << " " << F_fnmax << " " << Fmax_fnMean << '\n' << std::flush;
      }

      double Etrans, Erot;
      getKineticEnergy(Etrans, Erot);
      std::cout << "|   Kinetic energy:\n|   " << __printNamedCell(13, "Etrans:", 13, Etrans)
                << __printNamedCell(13, "Erot:", 13, Erot) << '\n';
      kineticEnergyFile << t << " " << Etrans << " " << Erot << '\n' << std::flush;
      std::cout << "+-----------------------------------------------------------------------------------\n\n";

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
      UpdateNL();
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

    // velocityVerletStep();
    IntegrationStep();

    t += dt;
    interConfC += dt;
    interVerletC += dt;

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
    @brief Increments the resultant forces and moments of the interacting bodies
           with the local forces and moments.
*/
void Rockable::incrementResultants(Interaction& I) {
  // Forces
  vec3r f = I.fn * I.n + I.ft;
  Particles[I.i].force += f;
  Particles[I.j].force -= f;

  // Moments
  vec3r Ci = (I.pos - Particles[I.i].pos);
  vec3r Cj = (I.pos - Particles[I.j].pos);
  Particles[I.i].moment += cross(Ci, f) + I.mom;
  Particles[I.j].moment += cross(Cj, -f) - I.mom;
}

/**
    @brief  Compute the particle accelerations (translations and rotations)

    The method computes actually 3 things:\n
      1. the interaction forces and moments with the force laws,\n
      2. the resultant forces and moments at the body centers,\n
      3. the axial and angular accelerations of the bodies
*/
void Rockable::accelerations() {
#ifdef ROT_MATRIX
  static mat9r P;
  static mat9r Pt;
#endif

  // Set resultant forces and moments to zero
#pragma omp parallel for default(shared)
  for (size_t i = 0; i < nDriven; ++i) {
    Particles[i].force.reset();
    Particles[i].acc.reset();
    Particles[i].moment.reset();
    Particles[i].arot.reset();
  }

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
    }
  }

  // Update the interactions (n, dn, pos and vel)
  // and apply the force law
  PerfTimer tm;
  activeInteractions.clear();
  interfacesToBreak.clear();

  // In the following loop, ALL interactions are updated,
  // including the interactions that are sticked or with positive dn
#pragma omp parallel for default(shared)
  for (size_t k = 0; k < Interactions.size(); ++k) {
    for (auto it = Interactions[k].begin(); it != Interactions[k].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      Interaction::UpdateDispatcher[it->type](*I, Particles[it->i], Particles[it->j]);
    }
  }

  // Some weighting relations can be established for the stiffnesses
  // of the cantacts that share the same body pair
  if (ctcPartnership.update != nullptr) ctcPartnership.update(*this);

    // Now the forces and moments are computed
#pragma omp parallel for default(shared)
  for (size_t k = 0; k < Interactions.size(); ++k) {
    for (auto it = Interactions[k].begin(); it != Interactions[k].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      if (it->dn < 0.0 || it->stick != nullptr) {
        forceLawPtr(*I);
      }
    }
  }

  // The increment of resultants (forces and moments) on the bodies
  // cannot be parallelised (because of possible conflicts).
  // Pointer to interactions are stored in the vector 'activeInteractions'
  for (size_t k = 0; k < Interactions.size(); ++k) {
    for (auto it = Interactions[k].begin(); it != Interactions[k].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      if (it->dn < 0.0 || it->stick != nullptr) {
        incrementResultants(*I);
        activeInteractions.push_back(I);
      }
    }
  }

#ifndef BREAK_ONCE
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
#endif

  timeInForceComputation += tm.getElapsedTimeSeconds();

  if (numericalDampingCoeff > 0.0) numericalDamping();

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
        // case Moments... _xrot_Mom_ _yrot_Mom_ _zrot_Mom_
    }
  }

// Finally compute the accelerations (translation and rotation) of the particles
#pragma omp parallel for default(shared)
  for (size_t i = nDriven; i < Particles.size(); ++i) {
    Particles[i].acc = Particles[i].force / Particles[i].mass;

#ifdef ROT_MATRIX
    Particles[i].Q.get_rot_matrix(P);
    Pt = P;
    Pt.transpose();
    vec3r omega = Pt * Particles[i].vrot;  // Express omega in the body framework
    vec3r M = Pt * Particles[i].moment;    // Express torque in the body framework
    vec3r domega(
        (M[0] - (Particles[i].inertia[2] - Particles[i].inertia[1]) * omega[1] * omega[2]) / Particles[i].inertia[0],
        (M[1] - (Particles[i].inertia[0] - Particles[i].inertia[2]) * omega[2] * omega[0]) / Particles[i].inertia[1],
        (M[2] - (Particles[i].inertia[1] - Particles[i].inertia[0]) * omega[0] * omega[1]) / Particles[i].inertia[2]);
    Particles[i].arot = P * domega;  // Express arot in the global framework
#else
    quat Qinv = Particles[i].Q.get_conjugated();
    vec3r omega = Qinv * Particles[i].vrot;  // Express omega in the body framework
    vec3r M = Qinv * Particles[i].moment;    // Express torque in the body framework
    vec3r domega(
        (M[0] - (Particles[i].inertia[2] - Particles[i].inertia[1]) * omega[1] * omega[2]) / Particles[i].inertia[0],
        (M[1] - (Particles[i].inertia[0] - Particles[i].inertia[2]) * omega[2] * omega[0]) / Particles[i].inertia[1],
        (M[2] - (Particles[i].inertia[1] - Particles[i].inertia[0]) * omega[0] * omega[1]) / Particles[i].inertia[2]);
    Particles[i].arot = Particles[i].Q * domega;  // Express arot in the global framework
#endif
  }

  // damping solutions based on the weighting of accelerations
  if (VelocityBarrier > 0.0) velocityBarrier();
  if (AngularVelocityBarrier > 0.0) angularVelocityBarrier();
}

/**
    @brief This is the Cundall damping solution
    @todo Change the implementation so that it acts on acceleration rather than forces/moments
*/
void Rockable::numericalDamping() {
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
    @brief Component-wise weighting of translation acceleration to limit the velocity components to
           the value 'VelocityBarrier'

    The weighting factor is equal to 1 when the velocity is of the order of zero,
    it tends towards 0 when the velocity approaches 'VelocityBarrier',
    and it becomes negative when the velocity exceeds 'VelocityBarrier'
    (it tends towards -1 when v tends towards infinity).
*/
void Rockable::velocityBarrier() {
  for (size_t i = 0; i < Particles.size(); ++i) {
    double vxratio = pow(fabs(Particles[i].vel.x / VelocityBarrier), VelocityBarrierExponent);
    Particles[i].acc.x *= (1.0 - vxratio) / (1.0 + vxratio);

    double vyratio = pow(fabs(Particles[i].vel.y / VelocityBarrier), VelocityBarrierExponent);
    Particles[i].acc.y *= (1.0 - vyratio) / (1.0 + vyratio);

    double vzratio = pow(fabs(Particles[i].vel.z / VelocityBarrier), VelocityBarrierExponent);
    Particles[i].acc.z *= (1.0 - vzratio) / (1.0 + vzratio);
  }
}

/**
    @brief Component-wise weighting of rotation acceleration to limit the velocity components to
           the value 'AngularVelocityBarrier'
*/
void Rockable::angularVelocityBarrier() {
  for (size_t i = 0; i < Particles.size(); ++i) {
    double vrotxratio = pow(fabs(Particles[i].vrot.x / AngularVelocityBarrier), AngularVelocityBarrierExponent);
    Particles[i].arot.x *= (1.0 - vrotxratio) / (1.0 + vrotxratio);

    double vrotyratio = pow(fabs(Particles[i].vrot.y / AngularVelocityBarrier), AngularVelocityBarrierExponent);
    Particles[i].arot.y *= (1.0 - vrotyratio) / (1.0 + vrotyratio);

    double vrotzratio = pow(fabs(Particles[i].vrot.z / AngularVelocityBarrier), AngularVelocityBarrierExponent);
    Particles[i].arot.z *= (1.0 - vrotzratio) / (1.0 + vrotzratio);
  }
}

// ==============================================================================================================
//  PROCESSING METHODS
// ==============================================================================================================

/**
    @brief   Computes the Axis Aligned Bounding Box (AABB) of a part of the scene.
    @remark  The AABB (paabb) of each particle is also updated in this method

    @param[in]   first   Smallest ID of particles (default value is 0)
    @param[in]   last    Largest ID of particles (default value corresponds to the last particle)
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

  // aabb.set_single(Particles[first].GlobVertex(0));
  // aabb.enlarge(Particles[first].MinskowskiRadius());
  aabb = paabb[first];

  for (size_t i = first + 1; i <= last; i++) {
    aabb.enlarge(paabb[i]);
  }
}

double Rockable::probeSolidFraction(AABB& aabb, size_t MCnstep) {
  if (MCnstep == 0) return -1.0;

  // select the concerned particles
  OBB zone;
  zone.center = 0.5 * (aabb.min + aabb.max);
  zone.extent.set(0.5 * (aabb.max.x - aabb.min.x), 0.5 * (aabb.max.y - aabb.min.y), 0.5 * (aabb.max.z - aabb.min.z));
  std::vector<size_t> pid;
  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].updateObb();
    if (zone.intersect(Particles[i].obb)) {
      pid.push_back(i);
    }
  }

  vec3r pt3;
  std::vector<double> vv(3);
  Mth::sobolSequence(-3, vv);  // Initialize the Sobol sequence
  size_t count = 0;
  for (size_t imc = 0; imc < MCnstep; ++imc) {
    Mth::sobolSequence(3, vv);
    pt3.set(aabb.min.x + vv[0] * (aabb.max.x - aabb.min.x), aabb.min.y + vv[1] * (aabb.max.y - aabb.min.y),
            aabb.min.z + vv[2] * (aabb.max.z - aabb.min.z));

    bool inSolid = false;
    for (size_t ii = 0; ii < pid.size(); ii++) {
      size_t i = pid[ii];
      vec3r ptTest = pt3 - Particles[i].pos;
      quat Qinv = Particles[i].Q.get_conjugated();
      ptTest = Qinv * ptTest;
      ptTest /= Particles[i].homothety;

      if (Particles[i].shape->inside(ptTest)) {
        inSolid = true;
        break;
      }
    }
    if (inSolid) count++;
  }

  return ((double)count / (double)MCnstep);
}

/**
    @brief  Get the global kinetic energy for translations and for rotations.

    @param[out]  Etrans  The translation kinetic energy (for particles id in range[first last])
    @param[out]  Erot    The ratational kinetic energy (for particles id in range[first last])
    @param[in]   first   Smallest ID of particles (default value is 0)
    @param[in]   last    Largest ID of particles (default value corresponds to the last particle)
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

void Rockable::getInteractionGroups(std::vector<size_t>& nbInt) {
  nbInt.clear();
  std::sort(activeInteractions.begin(), activeInteractions.end(), std::less<Interaction*>());
  size_t igrp = activeInteractions[0]->i;
  size_t jgrp = activeInteractions[0]->j;
  size_t nb = 1;
  for (size_t k = 1; k < activeInteractions.size(); k++) {
    if (activeInteractions[k]->i == igrp && activeInteractions[k]->j == jgrp) {
      nb++;
    } else {
      nbInt.push_back(nb);
      igrp = activeInteractions[k]->i;
      jgrp = activeInteractions[k]->j;
      nb = 1;
    }
  }
  nbInt.push_back(nb);
}

///
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
    @brief  Compute the critical time step by looping over all potential
            interactions, even those that are not active.
    @param[out]  dtc  minimum value of square root of meff/kn.
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
    @brief Compute the critical time step by looping over all currently active interactions
    @param[out]  dtc  minimum value of square root of meff/kn.
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
   @brief  Get the range of masses

   @param[out]  massMin  The lightest mass (for particles id in range[first last])
   @param[out]  massMax  The heaviest mass (for particles id in range[first last])
   @param[in]   first    Smallest ID of particles (default value is 0)
   @param[in]   last     Largest ID of particles (default value corresponds to the last particle)
*/
void Rockable::getMassRange(double& massMin, double& massMax, size_t first, size_t last) {
  if (last == 0) last = Particles.size() - 1;
  massMin = massMax = Particles[first].mass;
  for (size_t i = first + 1; i <= last; i++) {
    if (Particles[i].mass < massMin) massMin = Particles[i].mass;
    if (Particles[i].mass > massMax) massMax = Particles[i].mass;
  }
}

/**
   Fmax is the maximum value of resultant-force norms
   F_fnmax is the maximum value over all resultant-force norms divided by the maximum normal force acting on them
   Fmean is the mean value of the resultant-force norms
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
  if (n > 1) {
    Fmean /= (double)n;

    Fstddev = 0.0;
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
    @brief Get some statistics of the normal interaction forces

    @param[out] fnMin     Smallest normal contact force
    @param[out] fnMin     Biggest normal contact force
    @param[out] fnMean    Averaged value of normal contact forces
    @param[out] fnStddev  Standard deviation of normal contact forces

    @remarks Notice that this method needs that 'activeInteractions' is not empty.
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

/**
   @brief Get the set of clusters according to the cluster-ID that each particle has.
          Each cluster holds a vector of particle numbers (index in the vector @c Rockable::Particles)

   @param[out]  clusters  A vector of 'clusterParticles' that hold the clusterId and the list of involved particleId
*/
void Rockable::getClusters(std::vector<clusterParticles>& clusters) {
  clusters.clear();  // clear clusters if not empty
  std::set<clusterParticles> clusterSet;
  clusterParticles C;
  for (size_t i = nDriven; i < Particles.size(); i++) {
    C.clusterId = Particles[i].cluster;
    auto itClust = clusterSet.find(C);
    if (itClust == clusterSet.end()) {  // if not found
      auto p = clusterSet.insert(C);
      itClust = p.first;
    }

    // An element in the set clusterSet cannot be modified (because a set is
    // ordered) so we get a pointer to it, and change only the parameters that
    // will not affect the order
    clusterParticles* cp = const_cast<clusterParticles*>(std::addressof(*itClust));
    cp->particleId.push_back(i);
  }

  clusters.assign(clusterSet.begin(), clusterSet.end());

#if 0
  for (size_t c = 0 ; c < clusters.size() ; c++) {
    std::cout << std::endl;
    __SHOW(c);
    __SHOW( clusters[c].particleId.size() );
    for (size_t i = 0 ; i < clusters[c].particleId.size() ; i++) {
      __SHOW( clusters[c].particleId[i] );
    }
  }
#endif
}

/**
    @brief      Get the set of sub-parts (broken clusters).
    @attention  The clusterId in subParts is NOT the original cluster number (when it was not broken).
                Let say clusterId should be called partId in 'clusterParticles'

    @param[out]   subParts   A vector of 'clusterParticles'

    @remark This method has been carfully checked with Marta. It seems to work correctly
*/
void Rockable::getBrokenSubClusters(std::vector<clusterParticles>& subParts) {
  subParts.clear();  // clear clusters in case it's not empty

  std::vector<std::set<int>> subs;
  std::set<int> addedParticles;
  for (size_t i = 0; i < Interfaces.size(); i++) {
    for (auto it = Interfaces[i].begin(); it != Interfaces[i].end(); ++it) {
      size_t I = it->i;
      size_t J = it->j;
      if (Particles[I].cluster == Particles[J].cluster) {

        int afound = -1;
        if (!subs.empty()) {
          for (size_t a = subs.size(); a-- > 0;) {
            auto itI = subs[a].find(I);
            if (itI != subs[a].end()) {
              afound = a;
              break;
            }
            auto itJ = subs[a].find(J);
            if (itJ != subs[a].end()) {
              afound = a;
              break;
            }
          }
        }

        if (afound >= 0) {
          subs[afound].insert(I);
          subs[afound].insert(J);
          addedParticles.insert(I);
          addedParticles.insert(J);
        } else {
          std::set<int> s;
          s.insert(I);
          s.insert(J);
          subs.push_back(s);
          addedParticles.insert(I);
          addedParticles.insert(J);
        }

      }  // end 'if same cluster'
    }    // end 'for it'
  }      // end 'for i'

  // Merge sets that have common particle.
  // Remark: the merge of more than 2 sets is not possible because an interface can onle involve 2 particles
  std::vector<std::set<int>> subsFinal;
  for (size_t a = 0; a < subs.size(); a++) {
    bool hasBeenMerged = false;
    for (size_t b = a + 1; b < subs.size(); b++) {
      bool hasCommon = false;
      for (auto it = subs[a].begin(); it != subs[a].end(); ++it) {
        if (subs[b].find(*it) != subs[b].end()) {  // found
          hasCommon = true;
          break;
        }
      }
      if (hasCommon == true) {
        std::set<int> s;
        for (auto ia = subs[a].begin(); ia != subs[a].end(); ++ia) s.insert(*ia);
        for (auto ib = subs[b].begin(); ib != subs[b].end(); ++ib) s.insert(*ib);
        subsFinal.push_back(s);
        subs.erase(subs.begin() + b);
        hasBeenMerged = true;
        break;
      }
    }

    if (hasBeenMerged == false) {
      // add a to subsFinal
      std::set<int> s;
      for (auto ia = subs[a].begin(); ia != subs[a].end(); ++ia) s.insert(*ia);
      subsFinal.push_back(s);
    }
  }

  // Add parts composed of single particle
  for (int n = nDriven; n < (int)Particles.size(); n++) {
    if (addedParticles.find(n) == addedParticles.end()) {
      std::set<int> s;
      s.insert(n);
      subsFinal.push_back(s);
    }
  }

  for (size_t a = 0; a < subsFinal.size(); a++) {
    clusterParticles C;
    C.clusterId = a;
    C.particleId.assign(subsFinal[a].begin(), subsFinal[a].end());
    subParts.push_back(C);
  }

#if 0
  std::cout << " --------- " << std::endl;
  for (size_t c = 0 ; c < subParts.size() ; c++) {
    std::cout << std::endl;
    __SHOW(c);
    __SHOW( subParts[c].particleId.size() );
    for (size_t i = 0 ; i < subParts[c].particleId.size() ; i++) {
      __SHOW( subParts[c].particleId[i] );
    }
  }
#endif
}

// ==============================================================================================================
//  PRE-PROCESSING METHODS
// ==============================================================================================================

/**
   @brief  Apply a transformation (rotation followed by translation)
           to the bodies having number in a given range (boundaries included)

   @param[in]  center  Rotation point
   @param[in]  axis    Rotation axis
   @param[in]  angle   Rotation angle
   @param[in]  trans   Translation vector
*/
/*
void Rockable::rotateAndTranslateBlock(vec3r & center, vec3r & axis, double angle, vec3r & vector) {

}
*/

/**
   @brief Create sticked contacts between vertices (spheres) of particles
          that belong to the same cluster (same cluster-ID)

   @param[in]  epsilonDist  A small distance above which the 'glue' is not added
*/
void Rockable::stickVerticesInClusters(double epsilonDist) {
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].updateObb();
  }

  for (size_t i = nDriven; i < Particles.size(); i++) {

    OBB obbi = Particles[i].obb;
    obbi.enlarge(0.5 * DVerlet);

    for (size_t j = i + 1; j < Particles.size(); j++) {

      if (Particles[i].cluster != Particles[j].cluster) continue;

      OBB obbj = Particles[j].obb;
      obbj.enlarge(0.5 * DVerlet);

      // Check intersection
      if (obbi.intersect(obbj)) {

        for (size_t iv = 0; iv < Particles[i].shape->vertex.size(); iv++) {
          vec3r Vi = Particles[i].GlobVertex(iv);
          double Ri = Particles[i].MinskowskiRadius();
          for (size_t jv = 0; jv < Particles[j].shape->vertex.size(); jv++) {
            vec3r Vj = Particles[j].GlobVertex(jv);
            double Rj = Particles[j].MinskowskiRadius();

            double distSqr = norm2(Vj - Vi);
            double dMaxSqr = epsilonDist + (Ri + Rj);
            dMaxSqr *= dMaxSqr;

            if (distSqr < dMaxSqr) {
              // It is necessarily two free bodies (spheres at vertices) that interact
              double meff = (Particles[i].mass * Particles[j].mass) / (Particles[i].mass + Particles[j].mass);
              double en2 = dataTable.get(idEn2InnerBond, Particles[i].group, Particles[j].group);
              double kn = dataTable.get(idKnInnerBond, Particles[i].group, Particles[j].group);
              double Damp = 0.0;
              if (en2 > 0.0 && en2 < 1.0) {
                double logen = 0.5 * log(en2);
                double dampRate = -logen / sqrt(logen * logen + Mth::piSqr);
                Damp = dampRate * 2.0 * sqrt(kn * meff);
              } else if (en2 <= 0.0)
                Damp = 2.0 * sqrt(kn * meff);
              else
                Damp = 0.0;

              BreakableInterface BI_toInsert(i, j);
              BI_toInsert.isInner = 1;
              BI_toInsert.kn = kn;
              BI_toInsert.kt = dataTable.get(idKtInnerBond, Particles[i].group, Particles[j].group);
              BI_toInsert.kr = 0.0;
              BI_toInsert.fn0 = dataTable.get(idFn0InnerBond, Particles[i].group, Particles[j].group);
              BI_toInsert.ft0 = dataTable.get(idFt0InnerBond, Particles[i].group, Particles[j].group);
              BI_toInsert.mom0 = 1e10;
              BI_toInsert.power = dataTable.get(idPowInnerBond, Particles[i].group, Particles[j].group);
              std::pair<std::set<BreakableInterface>::iterator, bool> ret;
              ret = Interfaces[i].insert(BI_toInsert);
              BreakableInterface* BI = const_cast<BreakableInterface*>(std::addressof(*(ret.first)));
              Interaction I(i, j, vvType, iv, jv, Damp, BI);
              Interaction::UpdateDispatcher[vvType](I, Particles[i], Particles[j]);
              Interaction* Iptr = const_cast<Interaction*>(std::addressof(*((Interactions[i].insert(I)).first)));
              BI->concernedBonds.push_back(Iptr);
              activeInteractions.push_back(Iptr);
            }
          }
        }

      }  // if obb intersect
    }    // j
  }      // i

  std::sort(activeInteractions.begin(), activeInteractions.end(), std::less<Interaction*>());
}

/**
   @brief Create sticked interface between two sub-elements of two distinct particles
          that do not belong to the same cluster (different cluster-ID).
          The idea is to update the neighbor list and then replace the contacts
          between different clusters by a SINGLE sticked link.

   @param[in]  epsilonDist  A small distance below which the 'glue' is not added.
*/
void Rockable::stickClusters(double epsilonDist) {
  // In case the neighbor list has not been yet updated
  UpdateNL();

  // Associate pair(i, j) with vector of interactions
  std::map<std::pair<size_t, size_t>, std::vector<Interaction*>> ctc_packets;

  // Find Packets of contacts
  for (size_t k = 0; k < Interactions.size(); ++k) {
    for (auto it = Interactions[k].begin(); it != Interactions[k].end(); ++it) {
      size_t i = it->i;
      size_t j = it->j;

      if (glue_with_walls == false) {
        if (i < nDriven || j < nDriven) continue;  // no sticked link involving a driven body ('walls')
      }

      // Stick is only possible in-between different clusters
      if (Particles[i].cluster == Particles[j].cluster) continue;

      // update the interaction
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      bool valid = Interaction::UpdateDispatcher[it->type](*I, Particles[i], Particles[j]);

      if (valid && it->dn < epsilonDist) {
        if (j < i) std::swap(i, j);  // NORMALLY, NOT NECESSARY BECAUSE THE NEIGHBOR LIST IS CONSTRUCTED SO THAT i < j
        std::pair<size_t, size_t> duo(i, j);
        auto itf = ctc_packets.find(duo);
        if (itf == ctc_packets.end()) {  // if not found
          std::vector<Interaction*> v;
          auto p = ctc_packets.insert(std::pair<std::pair<size_t, size_t>, std::vector<Interaction*>>(duo, v));
          itf = p.first;
        }
        itf->second.push_back(I);
      }
    }
  }

  if (ctc_packets.empty()) {
    std::cout << msg::info() << "@Rockable::stickClusters, No possible glued points, ctc_packets is empty.\n"
              << msg::normal();
    return;
  }

  // Replace contact packets by single sticked link:
  for (auto p : ctc_packets) {

    // -- barycenter of each packet
    vec3r c;
    size_t nb = p.second.size();
    for (size_t ip = 0; ip < nb; ip++) {
      c += (p.second)[ip]->pos;
    }
    c /= nb;  // remark: we know that nb >= 1

    // -- select the contact point that is the closest (for each packet)
    double d2min = norm2((p.second)[0]->pos - c);
    size_t imin = 0;
    for (size_t ip = 1; ip < nb; ip++) {
      double d2 = norm2((p.second)[ip]->pos - c);
      if (d2 < d2min) {
        d2min = d2;
        imin = ip;
      }
    }

    //  -- Create the interface and plug it to the selected interaction (it will then be 'glued')
    size_t i = (p.second)[imin]->i;
    size_t j = (p.second)[imin]->j;
    if (j < i) std::swap(i, j);
    BreakableInterface BI_toInsert(i, j);
    BI_toInsert.isInner = 0;
    BI_toInsert.dn0 = (p.second)[imin]->dn;
    (p.second)[imin]->mom.reset();
    (p.second)[imin]->ft.reset();
    (p.second)[imin]->fn = 0.0;
    BI_toInsert.kn = dataTable.get(idKnOuterBond, Particles[i].group, Particles[j].group);
    BI_toInsert.kt = dataTable.get(idKtOuterBond, Particles[i].group, Particles[j].group);
    BI_toInsert.kr = dataTable.get(idKrOuterBond, Particles[i].group, Particles[j].group);
    BI_toInsert.fn0 = dataTable.get(idFn0OuterBond, Particles[i].group, Particles[j].group);
    BI_toInsert.ft0 = dataTable.get(idFt0OuterBond, Particles[i].group, Particles[j].group);
    BI_toInsert.mom0 = dataTable.get(idMom0OuterBond, Particles[i].group, Particles[j].group);
    BI_toInsert.power = dataTable.get(idPowOuterBond, Particles[i].group, Particles[j].group);

    std::pair<std::set<BreakableInterface>::iterator, bool> ret;
    ret = Interfaces[i].insert(BI_toInsert);
    BreakableInterface* BI = const_cast<BreakableInterface*>(std::addressof(*(ret.first)));

    Interaction* Iptr = const_cast<Interaction*>(std::addressof(*((p.second)[imin])));

    Iptr->stick = BI;
    BI->concernedBonds.push_back(Iptr);
  }

  // -- remove 'not-glued' contacts between clusters
  for (size_t k = 0; k < Interactions.size(); ++k) {
    for (auto it = Interactions[k].begin(); it != Interactions[k].end();) {
      size_t i = it->i;
      size_t j = it->j;
      if (i < nDriven || j < nDriven) {
        ++it;
        continue;
      }  // contact with walls are not erased
      if (Particles[i].cluster == Particles[j].cluster) {
        ++it;
        continue;
      }  // inner bonds are not erased
      if (it->stick != nullptr) {
        ++it;
        continue;
      }  // bonded links are not erased

      // erase the interaction (when not sticked)
      it = Interactions[k].erase(it);  // 'it' points now to the next iterator
    }
  }

  // activeInteractions needs to be re-set because 'saveConf' uses it
  activeInteractions.clear();
  for (size_t k = 0; k < Interactions.size(); ++k) {
    for (auto it = Interactions[k].begin(); it != Interactions[k].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      if (it->dn < 0.0 || it->stick != nullptr) {
        activeInteractions.push_back(I);
      }
    }
  }
}

/**
   @attention  We suppose that the Damp parameter has already been set
               in the corresponding Interactions
*/
void Rockable::copyParamsToInterfaces(std::string& isInnerStr) {
  int isInner = 0;
  if (isInnerStr == "inner") isInner = 1;

  for (size_t i = 0; i < Interfaces.size(); i++) {
    for (auto it = Interfaces[i].begin(); it != Interfaces[i].end(); ++it) {
      int icluster = Particles[it->i].cluster;
      int jcluster = Particles[it->j].cluster;
      if (isInner == 0 && icluster == jcluster) continue;
      if (isInner == 1 && icluster != jcluster) continue;

      // Since it is not possible to modify an element in a std::set
      // we break this restriction by defining and using the following pointer
      BreakableInterface* I = const_cast<BreakableInterface*>(std::addressof(*it));

      if (isInner) {
        I->kn = dataTable.get(idKnInnerBond, Particles[it->i].group, Particles[it->j].group);
        I->kt = dataTable.get(idKtInnerBond, Particles[it->i].group, Particles[it->j].group);
        I->fn0 = dataTable.get(idFn0InnerBond, Particles[it->i].group, Particles[it->j].group);
        I->ft0 = dataTable.get(idFt0InnerBond, Particles[it->i].group, Particles[it->j].group);
        I->power = dataTable.get(idPowInnerBond, Particles[it->i].group, Particles[it->j].group);
      } else {
        I->kn = dataTable.get(idKnOuterBond, Particles[it->i].group, Particles[it->j].group);
        I->kt = dataTable.get(idKtOuterBond, Particles[it->i].group, Particles[it->j].group);
        I->kr = dataTable.get(idKrOuterBond, Particles[it->i].group, Particles[it->j].group);
        I->fn0 = dataTable.get(idFn0OuterBond, Particles[it->i].group, Particles[it->j].group);
        I->ft0 = dataTable.get(idFt0OuterBond, Particles[it->i].group, Particles[it->j].group);
        I->mom0 = dataTable.get(idMom0OuterBond, Particles[it->i].group, Particles[it->j].group);
        I->power = dataTable.get(idPowOuterBond, Particles[it->i].group, Particles[it->j].group);
      }
    }
  }
  ParamsInInterfaces = 1;  // so that they are stored in the conf files (within interfaces)
}

/**
   @brief  Set the ratio kt/kn in interfaces
           (only for parameters that are stored in Interfaces)
*/
void Rockable::setStiffnessRatioInterfaces(double ratio) {
  for (size_t i = 0; i < Interfaces.size(); i++) {
    for (auto it = Interfaces[i].begin(); it != Interfaces[i].end(); ++it) {
      // Normally, we can't modify a value in a set because the order can be compromised.
      // But in this case, the value of kn and kt will not change the order.
      // To be able to modify kt, we thus make use of a pointer:
      double* Kt_ptr = (double*)&(it->kt);
      *Kt_ptr = ratio * it->kn;
    }
  }
}

/**
   Usage (in an input file):
   variableStickParams paramName inner/outer lambdaValue mValue 0/1
   example: variableStickParams fn0 inner 80 9 1
*/
void Rockable::setVariableStickParams(std::string& paramName, std::string& isInnerStr, double lambda, double m,
                                      bool timeSeeded) {
  std::default_random_engine generator;
  if (timeSeeded == true) {
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }
  std::weibull_distribution<> distribution(m, lambda);

  int shift = 0;
  if (paramName == "kn")
    shift = 0;
  else if (paramName == "kt")
    shift = 1;
  else if (paramName == "kr")
    shift = 2;
  else if (paramName == "fn0")
    shift = 3;
  else if (paramName == "ft0")
    shift = 4;
  else if (paramName == "mom0")
    shift = 5;

  int isInner = 0;
  if (isInnerStr == "inner") isInner = 1;

  for (size_t i = 0; i < Interfaces.size(); i++) {
    for (auto it = Interfaces[i].begin(); it != Interfaces[i].end(); ++it) {
      int icluster = Particles[it->i].cluster;
      int jcluster = Particles[it->j].cluster;
      if (isInner == 0 && icluster == jcluster) continue;
      if (isInner == 1 && icluster != jcluster) continue;

      double* seek = (double*)&(it->kn);
      *(seek + shift) = distribution(generator);
    }
  }

  ParamsInInterfaces = 1;  // so that they are stored in the conf files (within interfaces)
}

/**
   @brief  Set randomly oriented velocity vectors (the magnitudes are all the same).
           The velocities of the driven bodies are not modified
*/
void Rockable::randomlyOrientedVelocities(double velocityMagnitude) {
  quat q;
  q.randomize(true);
  vec3r u(velocityMagnitude, 0.0, 0.0);
  for (size_t i = nDriven; i < Particles.size(); i++) {
    q.randomize();
    Particles[i].vel = q * u;
  }
}

/**
   @brief Set random orientation to the clusters
   @param[in]  velocityMagnitude  Magnitude of all velocities (only orientations change)
   @param[in]  opt                An option. if opt = 1 then all the velocity vectors
                                  will be oriented towards negative y (downward)
*/
void Rockable::randomlyOrientedVelocitiesClusters(double velocityMagnitude, int opt) {
  std::vector<clusterParticles> clusters;
  getClusters(clusters);

  quat q;
  q.randomize(true);
  vec3r u(velocityMagnitude, 0.0, 0.0);
  for (size_t c = 0; c < clusters.size(); c++) {
    q.randomize();
    vec3r v = q * u;
    if (opt == 1) v.y = -fabs(v.y);
    for (size_t i = 0; i < clusters[c].particleId.size(); i++) {
      Particles[clusters[c].particleId[i]].vel = v;
    }
  }
}

/**
   @brief Set Uniform distribution of homothety of a sub-set of particles
   @param[in]  idFirst     ID-Number of the first particle
   @param[in]  idLast      ID-Number of the last particle
   @param[in]  hmin        Minimum of the homothety range
   @param[in]  hmax        Maximum of the homothety range
   @param[in]  timeSeeded  if true, the random generator is seeded with current time

   Usage in input conf-file:
   homothetyRange idFirst idLast hmin hmax timeSeeded(0/1)
*/
void Rockable::homothetyRange(size_t idFirst, size_t idLast, double hmin, double hmax, bool timeSeeded) {
  std::default_random_engine generator;
  if (timeSeeded == true) {
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }
  std::uniform_real_distribution<> distribution(hmin, hmax);

  for (size_t i = idFirst; i <= idLast; i++) {
    double h = distribution(generator);
    Particles[i].homothety = h;
    Particles[i].mass = (h * h * h * Particles[i].shape->volume) * properties.get(idDensity, Particles[i].group);
    Particles[i].inertia = (h * h * Particles[i].shape->inertia_mass) * Particles[i].mass;
  }
}

/**
   Usage in input conf-file:
   particlesClonage idFirst idLast dX dY dZ
*/
void Rockable::particlesClonage(size_t idFirst, size_t idLast, vec3r& translation) {
  if (Particles.empty()) return;

  int numClusterMax = Particles[0].cluster;
  for (size_t i = 1; i < Particles.size(); i++) {
    if (Particles[i].cluster > numClusterMax) numClusterMax = Particles[i].cluster;
  }

  Particle P;
  for (size_t i = idFirst; i <= idLast; i++) {
    P.group = Particles[i].group;
    P.cluster = Particles[i].cluster - idFirst + numClusterMax + 1;
    P.homothety = Particles[i].homothety;
    P.pos = Particles[i].pos + translation;
    P.vel = Particles[i].vel;
    P.acc = Particles[i].acc;
    P.Q = Particles[i].Q;
    P.vrot = Particles[i].vrot;
    P.arot = Particles[i].arot;
    P.shape = Particles[i].shape;
    P.homothety = Particles[i].homothety;
    P.mass = Particles[i].mass;
    P.inertia = Particles[i].inertia;
    P.obb = Particles[i].obb;
    P.force = Particles[i].force;
    P.moment = Particles[i].moment;
    Particles.push_back(P);
  }

  if (Interactions.size() != Particles.size()) Interactions.resize(Particles.size());
  if (Interfaces.size() != Particles.size()) Interfaces.resize(Particles.size());
}
