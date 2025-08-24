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
#include <unordered_map>
#include <utility>
#include <vector>

#define FMT_HEADER_ONLY
#include "fmtLogger.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef ROCKABLE_GIT_TAG
#define ROCKABLE_GIT_TAG "unknown"
#endif

// toofus headers
#include "AABB.hpp"
#include "Compliance.hpp"
#include "DataTable.hpp"
#include "Mth.hpp"
#include "PerfTimer.hpp"
#include "Properties.hpp"
#include "Tempo.hpp"
#include "common.hpp"
#include "factory.hpp"
#include "fileTool.hpp"
#include "kwParser.hpp"
#include "linkCells.hpp"
#include "message.hpp"

// local headers
#include "BodyForces/BodyForce.hpp"
#include "BreakableInterface.hpp"
#include "ContactPartnership.hpp"
#include "DataExtractors/DataExtractor.hpp"
#include "DrivingSystem.hpp"
#include "ForceLaws/ForceLaw.hpp"
#include "Interaction.hpp"
#include "Particle.hpp"
#include "PeriodicCell.hpp"
#include "PreproCommands/PreproCommand.hpp"
#include "ProfilingTools.hpp"
#include "SpringJoint.hpp"
#include "clusterParticles.hpp"

#include "MutexVector.hpp"
#include "Traversals.hpp"
#include "preComputeUpdate.hpp"
#ifdef ROCKABLE_ENABLE_BOUNDARY
#include "Boundaries/Boundary.hpp"
#endif

class Rockable {
 public:
  std::vector<Particle> Particles;                           ///< The particles
  std::vector<std::set<Interaction>> Interactions;           ///< The interactions (contacts and potential contacts)
  std::vector<std::set<BreakableInterface>> Interfaces;      ///< The interfaces
  std::vector<std::vector<Interaction*>> m_vecInteractions;  ///< Vector of interaction pointers
  std::vector<Interaction*> activeInteractions;  ///< Hold a pointer to the contact interactions that are active
  std::vector<SpringJoint> joints;               ///< The spring joints

  ContactPartnership ctcPartnership;  ///< Model to weight the stiffnesses of contacts
                                      ///< in case of multiple contacts between two particles

  DrivingSystem System;  ///< The system driver acting on the nDriven first particles
  size_t nDriven;        ///< Number of sphero-polyhedra that are fix at the beginning of the list of bodies

  std::vector<DataExtractor*> dataExtractors;  ///< Some data to be saved in files according to given processings

  std::vector<Tempo<double>> Tempos;  ///< Useful for force ramps for example

  BodyForce* bodyForce;  ///< We can add body-forces and body-moments in addition to gravity

  // Shape library
  std::string shapeFile;                  ///< Name of the file that contain the library of shapes
  std::vector<Shape> Shapes;              ///< Loaded library of shapes
  std::map<std::string, size_t> shapeId;  ///< Associate a name of shape with its id in the vector 'Shapes'

#ifdef ROCKABLE_ENABLE_PERIODIC
  int usePeriodicCell;  ///< Flag indicating if periodic cell is used
  PeriodicCell Cell;    ///< The periodic cell
#endif

  int useSoftParticles;  ///< Flag indicating if soft particles are used
  Compliance Cinv;       ///< The inverse compliance

  // Time parameters
  double t;                  ///< Current Time
  double tmax;               ///< End time
  double dt;                 ///< Time increment
  int computationStopAsked;  ///< Ask to break the integration loop if positive

  // Simulation flow
  double interVerlet;  ///< Time interval between each update of the neighbor-list
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
  double numericalDampingCoeff;    ///< Coefficient of the so-called Cundall damping (0 = disabled)
  double velocityBarrier;          ///< Velocity Barrier by weighting the particle accelerations (0 = disabled)
  double angularVelocityBarrier;   ///< Angular velocity Barrier by weighting the particle accelerations (0 = disabled)
  double velocityBarrierExponent;  ///< The exponent to be used with velocityBarrier
  double angularVelocityBarrierExponent;  ///< The exponent to be used with angularVelocityBarrier

  // Other parameters
  int iconf;                ///< Current configuration ID
  AABB aabb;                ///< Axis Aligned Bounding Box (AABB) that surrounds some free bodies.
  std::vector<AABB> paabb;  ///< AABBs of all particles
  vec3r gravity;            ///< Gravity acceleration
  int paramsInInterfaces;   ///< Flag saying if the parameters for sticked blonds
                            ///< are embedded within the interfaces
  bool glue_with_walls;     ///< flag to say if the glue can be also set between free and controlled particles

  vec3r cellMinSizes;      ///< Sizes of the linked cells in the 3 directions x, y ad z
  int boxForLinkCellsOpt;  ///< Option for defining the master cell when linkCells' stategy in employed

  std::map<std::string, std::string> optionNames;  ///< All the options refered to as a string and set
                                                   ///< as a string also

  ForceLaw* forceLaw;       ///< User defined force law
  linkCells m_linkCells;    ///< Used to dectect contact by using the linkCell method
  Traversals m_traversals;  ///< Used to store different cell traversal.

  double preventCrossingLength;  ///< The length used to prevent crossing

  // Ctor
  Rockable();
  void ExplicitRegistrations();  ///< Performs explicit registrations of various features.

  // Initialization methods
  void setVerboseLevel(int v);                         ///< Defines verbosity with an integer
  void setVerboseLevel(const std::string& levelName);  ///< Defines verbosity with a string label
  void initOutputFiles();                              ///< Open the output files
  void setInteractive(bool imode = true);              ///< Set in a computation mode
  bool isInteractive() const;                          ///< Set in a visualization mode
  void showBanner();                                   ///< Display a banner about the code
  void setOpenMPThreads(int nbThreads = 1);            ///< Set the number of threads

  void initialChecks();  ///< Checks before runing a computation

  // Core DEM methods
  void velocityVerletStep();  ///< Make a time increment with the velocity-Verlet scheme
  void EulerStep();           ///< Make a time increment with the explicit-Euler scheme
  void BeemanStep();          ///< Make a time increment with the Beeman scheme
  void RungeKutta4Step();     ///< Make a time increment with the Runge-Kutta scheme (4th order)
  void initIntegrator();      ///< Set additionally stored data required by some integrator
  void integrate();           ///< Simulation flow (make time increments and check for updates or saving)

  void accelerations();                    ///< Compute accelerations
  void incrementResultants(Interaction&);  ///< Project force and moment on the interacting particles
#ifdef ROCKABLE_ENABLE_PERIODIC
  void incrementPeriodicCellTensorialMoment(Interaction&);  ///< Update the tensorial moment of interacting particles
  void realToReducedKinematics();  ///< Go from real coordinates to reduced ones for positions and velocities
  void reducedToRealKinematics();  ///< Go from reduced coordinates to real ones for positions and velocities
#endif
  std::function<void()> integrationStep;  ///< Pointer function for  tion
  void setIntegrator(std::string& Name);  ///< Select the time-integration scheme

 private:
  // These methods are subparts of the method accelerations()
  void initialise_particle_forces_and_moments();  ///< Initialise the forces and moments
  void update_interactions();                     ///< Update the interactions
  void update_interactions_boundary();            ///< Update the interactions
  void build_activeInteractions();                ///< Build the active interactions
  void compute_forces_and_moments();              ///< Compute forces and moments
  void compute_resultants();                      ///< Compute the resulting forces and moments
  void compute_SpringJoints();                    ///< Compute the spring joints
  void check_breakage_of_interfaces();            ///< ____
  void breakage_of_interfaces();                  ///< Break the interfaces
  void compute_accelerations_from_resultants();   ///< Compute accelerations

  bool computeInterArray();
  void applyForces();
  void buildPreComputeUpdate();
  //  void updateDispatcher();
  void computeForcesAndMoments();
  void updateForcesAndMoments();
  std::vector<std::vector<Interaction*>> testInter;
  MutexVector InteractionsMutex;  ///< mutex to synchornize concurrent writes to Interactions in AddOrRemoveInteractions
  preComputeUpdate m_reorg;       ///< Some stuff to update in parallel forces and moments in OpenMP

 public:
  // Core CD method (TODO)
  // ...

  // Save/Load methods
  kwParser parser;
  void clearMemory();            ///< Clear the Particles and Interactions
  void initParser();             ///< The place where the parsing commands are defined
  void saveConf(int i);          ///< Save the current configuration in a file named confX, where X=i
  void saveConf(const char*);    ///< Save a configuration from a conf-file
  void loadConf(int i);          ///< Load the current configuration in a file named confX, where X=i
  void loadConf(const char*);    ///< Load a configuration from a conf-file
  void loadShapes(const char*);  ///< Load a shape library from a shape-file
  void readDataExtractors();     ///< Read the file dataExtractors.txt

  void console_run(const std::string& confFileName);

  void Interactions_from_set_to_vec();
  void UpdateNL_bruteForce();           ///< Brute-force approach
  void UpdateNL_linkCells();            ///< Link-cells approach
  std::function<void()> UpdateNL;       ///< Pointer funtion to the updateNL
                                        ///< (update the neighbor-list) method
  void setUpdateNL(std::string& Name);  ///< select the Neighbor list updator

#ifdef ROCKABLE_ENABLE_BOUNDARY
  Boundary* boundary{nullptr};
  std::vector<std::set<InteractionBoundary*>> InteractionsBoundary;
  std::vector<InteractionCylinder*> activeInteractionsBoundary;
  std::function<bool(InteractionCylinder&)> forceLawCylinderPtr;
  // void incrementResultantsPiCylinder(InteractionCylinder&);
  // void UpdateNL_bruteForce_Cylinder();
  // bool forceLawAvalanchesCylinder(InteractionCylinder&);
  // bool forceLawHertzCylinder(InteractionCylinder&);
  // bool forceLawCohesiveForceCylinder(InteractionCylinder&);
#endif

  // Processing methods
  void computeAABB(size_t first = 0, size_t last = 0);  ///< Compute Axis Aligned Bounding Box of a part of the sample
  void getCriticalTimeStep(double& dtc);                ///< Get the critical time step
  void getCurrentCriticalTimeStep(double& dtc);         ///< Get the current critical time step
  void estimateCriticalTimeStep(double& dtc);           ///< Estimate the critical time step
  void getResultantQuickStats(double& Fmax, double& F_fnmax, double& Fmean, double& Fstddev, size_t first = 0,
                              size_t last = 0);  ///> Get the quick statistics
  void getInteractionQuickStats(double& fnMin, double& fnMax, double& fnMean,
                                double& fnStddev);                                         ///> Get the quick statistics
  void getKineticEnergy(double& Etrans, double& Erot, size_t first = 0, size_t last = 0);  ///> Get the kinetic energy

  // =============================================================================================================

  void velocityControlledDrive();      ///< Impose the velocities of driven bodies
  void numericalDamping();             ///< The Cundall damping solution
  void applyVelocityBarrier();         ///< The velocity barrier solution
  void applyAngularVelocityBarrier();  ///< The velocity barrier solution for rotations
  void dynamicCheckUpdateNL();         ///< Ask for NL reconstruction if the maximum
                                       ///< displacement of rotation is more than given values

  int AddOrRemoveInteractions_bruteForce(size_t i, size_t j, double dmax);  ///< Brute-force approach
  int AddOrRemoveInteractions_OBBtree(size_t i, size_t j, double dmax);     ///< OBB-tree approach
  std::function<int(size_t, size_t, double)> AddOrRemoveInteractions;       ///< (Pointer function) Add or remove an
                                                                       ///< interaction according to the distance dmax
  void setAddOrRemoveInteractions(std::string& Name);  ///> Select the add or remove interaction

  void readLawData(std::istream&, size_t id);             ///< Helper method to read a law in loadConf
  void writeLawData(std::ostream&, const char* parName);  ///< Helper method to write law in saveConf

  bool interactiveMode;  ///< computation (false) or visualization (true) modes

  // Some predefined identifiers to get (quickly) data from input tables
  size_t idDensity;  ///< Identifier of the density parameter

  size_t idKnContact;     ///< Identifier of normal stiffness for contact
  size_t idEn2Contact;    ///< Identifier of normal energy-restitution coefficient
  size_t idKtContact;     ///< Identifier of tangential stiffness for contact
  size_t idMuContact;     ///< Identifier of friction coefficient for contact
  size_t idKrContact;     ///< Identifier of angular stiffness for contact
  size_t idMurContact;    ///< Identifier of rolling-friction coefficient for contact
  size_t idViscTContact;  ///< Identifier of tangential viscuous parameter for contact

  size_t idKnInnerBond;    ///< Identifier of normal stiffness for inner-bonds
  size_t idKtInnerBond;    ///< Identifier of tangential stiffness for inner-bonds
  size_t idKrInnerBond;    ///< Identifier of angular stiffness for inner-bonds
  size_t idEn2InnerBond;   ///< Identifier of normal energy-restitution for inner-bonds
  size_t idFn0InnerBond;   ///< Identifier of normal threshold force for inner-bonds
  size_t idFt0InnerBond;   ///< Identifier of tangential threshold force for inner-bonds
  size_t idMom0InnerBond;  ///< Identifier of threshold bendind-moment for inner-bonds
  size_t idPowInnerBond;   ///< Identifier of power in yield function for inner-bonds

  size_t idKnOuterBond;    ///< Identifier of normal stiffness for outer-bonds
  size_t idKtOuterBond;    ///< Identifier of tangential stiffness for outer-bonds
  size_t idKrOuterBond;    ///< Identifier of angular stiffness for outer-bonds
  size_t idEn2OuterBond;   ///< Identifier of normal energy-restitution for outer-bonds
  size_t idFn0OuterBond;   ///< Identifier of normal threshold force for outer-bonds
  size_t idFt0OuterBond;   ///< Identifier of tangential threshold force for outer-bonds
  size_t idMom0OuterBond;  ///< Identifier of threshold bendind-moment for outer-bonds
  size_t idPowOuterBond;   ///< Identifier of power in yield function for outer-bonds

  size_t idGcInnerBond;  ///< Identifier of ____
  size_t idGcOuterBond;  ///< Identifier of ____

  // Postponed breakage of BreakableInterfaces
  std::set<BreakableInterface*> interfacesToBreak;  ///< Pointers to interfaces to be broken

  std::set<BreakableInterface*> breakableInterfaces;
  double totalBrokenArea{0.0};  // for tracking area of interface liberated during breakage

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
  double dt_6;   ///< Time-step divided by 6
  double dt2_2;  ///< Half the squared time-step
  double dt2;    ///< Double of time-step
  double dt2_8;  ///< Squared time-step divided by 8
  double dt2_6;  ///< Squared time-step divided by 6

  // input files
  std::string m_path = std::string();  ///< Default path to load our input files is the current directory

  // output files
  std::ofstream perfFile;           ///< to store performance data in the course of a computation
  std::ofstream staticBalanceFile;  ///< to store static balance data in the course of a computation
  std::ofstream kineticEnergyFile;  ///< to store kinetic energy data in the course of a computation
};

#endif /* end of include guard: ROCKABLE_HPP */
