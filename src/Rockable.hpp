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


#ifndef ROCKABLE_HPP
#define ROCKABLE_HPP

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// common headers
#include "AABB.hpp"
#include "DataTable.hpp"
#include "Mth.hpp"
#include "PerfTimer.hpp"
#include "Properties.hpp"
#include "Tempo.hpp"
#include "factory.hpp"
#include "fileTool.hpp"
#include "kwParser.hpp"
#include "linkCells.hpp"
#include "message.hpp"

// local headers
#include "BreakableInterface.hpp"
#include "ContactPartnership.hpp"
#include "DataExtractor.hpp"
#include "DrivingSystem.hpp"
#include "Interaction.hpp"
#include "Particle.hpp"
#include "clusterParticles.hpp"
#include "BodyForce.hpp"

class Rockable {
 public:
  std::vector<Particle> Particles;                        ///< The particles
  std::vector<std::set<Interaction> > Interactions;       ///< The interactions
                                                          ///< (contacts and potential contacts)
  std::vector<std::set<BreakableInterface> > Interfaces;  ///< The interfaces
  std::vector<Interaction*> activeInteractions;           ///< Hold a pointer to the
                                                          ///< contact interactions that are active
  ContactPartnership ctcPartnership;                      ///< Model to weight the stiffnesses of contacts
                                                          ///< in case of multiple contacts between two particles
  DrivingSystem System;                                   ///< The system driver acting on the nDriven first particles
  size_t nDriven;                                         ///< Number of sphero-polyhedra that are fix at the beginning
                                                          ///< of the list of bodies

  std::vector<DataExtractor*> dataExtractors;  ///< Some data to be saved in
                                               ///< files according to given
                                               ///< processings
  std::vector<Tempo<double> > Tempos;    ///< useful for force ramps for example
  
  BodyForce *bodyForce;                  ///< We can add body-forces and body-moments
                                         ///< in addition to gravity 

  // Shape library
  std::string shapeFile;                  ///< Name of the file that contain the library of shapes
  std::vector<Shape> Shapes;              ///< Loaded library of shapes
  std::map<std::string, size_t> shapeId;  ///< Associate a name of shape with
                                          ///< its id in the vector 'Shapes'

  // Time parameters
  double t;     ///< Current Time
  double tmax;  ///< End time
  double dt;    ///< Time increment

  // Simulation flow
  double interVerlet;  ///< Time interval between each update of the
                       ///< neighbor-list
  double interConf;    ///< Time interval between the CONF files

  // Neighbor list
  double DVerlet;  ///< Distance of Verlet for obb intersections
  double dVerlet;  ///< Distance of Verlet for sub-body proximities

  // Mechanical properties
  Properties properties;  ///< Density of the particles (defined by group)
  DataTable dataTable;    ///< Table that hold the parameters for the force-laws

  // Neighbor list
  int dynamicUpdateNL;   ///< Flag to activate or not the dynamic update of the neighbor list
  double dispUpdateNL;   ///< The moving distance that trigger an update of the neighbor list
  double angleUpdateNL;  ///< The moving rotation-angle that trigger an update of the neighbor list

  // Numerical Damping
  double numericalDampingCoeff;  ///< Coefficient of the so-called Cundall damping

  // Other parameters
  int iconf;                ///< Current configuration ID
  AABB aabb;                ///< Axis Aligned Bounding Box (AABB) of the entire sample
  std::vector<AABB> paabb;  ///< AABB associated to each particle (evolves during computation)
  vec3r gravity;            ///< Gravity acceleration
  int ParamsInInterfaces;   ///< Flag saying if the parameters for sticked blonds
                            ///< are embedded within the interfaces
  bool glue_with_walls;     ///< flag to say if the glue can be also set between free and controlled particles 

  vec3r cellMinSizes;      ///< Sizes of the linked cells in the 3 directions x, y ad z
  int boxForLinkCellsOpt;  ///< Option for defining the master cell when linkCells' stategy in employed


  // Ctor
  Rockable();

  // Initialization methods
  void initOutputFiles();                  ///< Open the output files
  void setInteractive(bool imode = true);  ///< Set in a computation mode
  bool isInteractive() const;              ///< Set in a visualization mode
  void showBanner();                       ///< Display a banner about the code

  void initialChecks();  ///< Checks before runing a computation

  // Core DEM methods
  void velocityVerletStep();               ///< Make a time increment with the
                                           ///< velocity-Verlet scheme
  void integrate();                        ///< Simulation flow (make time increments and check for
                                           ///< updates or saving)
  void accelerations();                    ///< Compute accelerations (both for particles and the
                                           ///< periodic-cell)
  void incrementResultants(Interaction&);  ///< Project force and moment on the
                                           ///< interacting particles

  // Save/Load methods
  void clearMemory();            ///< Clear the Particles and Interactions
  void saveConf(int i);          ///< Save the current configuration in a file named
                                 ///< confX, where X=i
  void loadConf(const char*);    ///< Load a configuration from a conf-file
  void loadShapes(const char*);  ///< Load a shape library from a shape-file
  void readDataExtractors();     ///< Read the file dataExtractors.txt

  void UpdateNL_bruteForce();      ///< Brute-force approach
  void UpdateNL_linkCells();       ///< Link-cells approach
  std::function<void()> UpdateNL;  ///< Pointer funtion to the updateNL (update
                                   ///< the neighbor-list) method

  bool forceLawStickedLinks(Interaction&);        ///< Force law for 'glued' bodies
  bool forceLawAvalanches(Interaction&);          ///< Force law specific for rock avalanches
  std::function<bool(Interaction&)> forceLawPtr;  ///< Pointer funtion to the force-law used

  // Processing methods
  void computeAABB(size_t first = 0, size_t last = 0);  ///< Compute Axis Aligned Bounding
                                                        ///< Box of a part of the sample
  void getCriticalTimeStep(double& dtc);
  void getCurrentCriticalTimeStep(double& dtc);
  void estimateCriticalTimeStep(double& dtc);
  void getMassRange(double& massMin, double& massMax, size_t first = 0, size_t last = 0);
  void getResultantQuickStats(double& Fmax, double& F_fnmax, double& Fmean, double& Fstddev, size_t first = 0,
                              size_t last = 0);
  void getInteractionQuickStats(double& fnMin, double& fnMax, double& fnMean, double& fnStddev);
  void getKineticEnergy(double& Etrans, double& Erot, size_t first = 0, size_t last = 0);
  void getClusters(std::vector<clusterParticles>& clusters);
  void getBrokenSubClusters(std::vector<clusterParticles>& subclusters);
  void getInteractionGroups(std::vector<size_t>& nbInt);

  // Pre-processing methods
  void stickVerticesInClusters(double);  ///< Create bonds between vertices if
                                         ///< they belong to the same cluster
  void stickClusters(double);            ///< Create bonds in-between clusters
  void copyParamsToInterfaces(std::string& isInnerStr);
  void setStiffnessRatioInterfaces(double ratio);

  void setVariableStickParams(std::string& paramName, std::string& isInnerStr, double lambda, double m,
                              bool timeSeeded);
  void randomlyOrientedVelocities(double velocityMagnitude);  ///< Set a randomly oriented
                                                              ///< velocity to all free bodies

  void randomlyOrientedVelocitiesClusters(double velocityMagnitude,
                                          int opt = 0);  ///< Set a randomly oriented velocity
                                                         ///< to all free clusters (same
                                                         ///< velocity for the bodies belonging
                                                         ///< to the same cluster)
  void homothetyRange(size_t idFirst, size_t idLast, double hmin, double hmax, bool timeSeeded);
  void particlesClonage(size_t idFirst, size_t idLast, vec3r& translation);

  // =============================================================================================================
 private:
  void numericalDamping();      ///< The Cundall damping solution
  void dynamicCheckUpdateNL();  ///< Ask for NL reconstruction if the maximum
                                ///< displacement of rotation is more than given
                                ///< values

  int AddOrRemoveInteractions_bruteForce(size_t i, size_t j,
                                         double dmax);                   ///< Brute-force approach
  int AddOrRemoveInteractions_OBBtree(size_t i, size_t j, double dmax);  ///< OBB-tree approach
  std::function<int(size_t, size_t, double)> AddOrRemoveInteractions;    ///< (Pointer function) Add or remove an
                                                                         ///< interaction according to the distance dmax

  void getInteractingGroups(Interaction& I, int& g1, int& g2);

  void readLawData(std::istream&, size_t id);             ///< Helper method to read a law in loadConf
  void writeLawData(std::ostream&, const char* parName);  ///< Helper method to
                                                          ///< write law in
                                                          ///< saveConf

  std::map<std::string, std::string> optionNames;  ///< All the options refered to as a string and set
                                                   ///< as a string also
  bool interactiveMode;                            ///< computation (false) or visualization (true) modes

  // Some predefined identifiers to get (quickly) data from input tables
  size_t idDensity;  ///< Identifier of the density parameter

  size_t idKnContact;   ///< Identifier of normal stiffness for contact
  size_t idEn2Contact;  ///< Identifier of normal energy-restitution coefficient
  size_t idKtContact;   ///< Identifier of tangential stiffness for contact
  size_t idMuContact;   ///< Identifier of friction coefficient for contact
  size_t idKrContact;   ///< Identifier of angular stiffness for contact
  size_t idMurContact;  ///< Identifier of rolling-friction coefficient for contact

  size_t idKnInnerBond;   ///< Identifier of normal stiffness for inner-bonds
  size_t idKtInnerBond;   ///< Identifier of tangential stiffness for inner-bonds
  size_t idEn2InnerBond;  ///< Identifier of normal energy-restitution for inner-bonds
  size_t idFn0InnerBond;  ///< Identifier of normal threshold force for inner-bonds
  size_t idFt0InnerBond;  ///< Identifier of tangential threshold force for inner-bonds
  size_t idPowInnerBond;  ///< Identifier of power in yield function for inner-bonds

  size_t idKnOuterBond;    ///< Identifier of normal stiffness for outer-bonds
  size_t idKtOuterBond;    ///< Identifier of tangential stiffness for outer-bonds
  size_t idKrOuterBond;    ///< Identifier of angular stiffness for outer-bonds
  size_t idEn2OuterBond;   ///< Identifier of normal energy-restitution for outer-bonds
  size_t idFn0OuterBond;   ///< Identifier of normal threshold force for outer-bonds
  size_t idFt0OuterBond;   ///< Identifier of tangential threshold force for outer-bonds
  size_t idMom0OuterBond;  ///< Identifier of threshold bendind-moment for outer-bonds
  size_t idPowOuterBond;   ///< Identifier of power in yield function for outer-bonds

  // Postponed breakage of BreakableInterfaces
  std::set<BreakableInterface*> interfacesToBreak;  ///< Pointers to interfaces
                                                    ///< to be broken

  // dynamic update of the Neighbor-List (NL)
  bool needUpdate;              ///< when true, an updateNL will be done at the next time increment
  double maxDeltaPos;           ///< Overall maximum displacement since the last updateNL
  double maxDeltaRot;           ///< Overall maximum rotation angle since the last updateNL
  std::vector<vec3r> deltaPos;  ///< displacement vector of each particle (since last updateNL)
  std::vector<quat> deltaQ;     ///< rotation quaternion of each particle (since last updateNL)

  // CPU time measurments
  double timeInUpdateNL;          ///< Used to measure elapsed time by the update of NL
  double timeInForceComputation;  ///< Used to measure elapsed time by the force laws

  // Counters for the simulation flow
  double interVerletC;  ///< Counter between Verlet-updates
  double interConfC;    ///< Counter between savings of conf-files

  // time-step constants
  double dt_2;   ///< Half the time-step
  double dt2_2;  ///< Half the squared time-step

  // output files
  std::ofstream perfFile;           ///< to store performance data in the course of a computation
  std::ofstream staticBalanceFile;  ///< to store static balance data in the course of a computation
  std::ofstream kineticEnergyFile;  ///< to store kinetic energy data in the course of a computation
};

#endif /* end of include guard: ROCKABLE_HPP */
