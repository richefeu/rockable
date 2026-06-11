#include "confedit.hpp"

// TODO rename init_snippets()
void init_snippets() {

  comp["Time"] = "t _value_";
  comp["End time"] = "tmax _value_";

  comp["Time step"] = "dt _value_";
  comp["Elapsed time between NL updates"] = "interVerlet _duration_";
  comp["Elapsed time between conf files"] = "interConf _duration_";
  comp["Proximity distance for rigid bodies"] = "DVerlet _distance_";
  comp["Proximity distance for sub-elements"] = "dVerlet _distance_";
  comp["Density of all bodies"] = "density _groupNumber_ _density_";
  comp["Gravity"] = "gravity _value_";
  comp["Parameters stored in the interfaces"] = "ParamsInInterfaces _0/1_";
  comp["Update NL dynamically"] = "dynamicUpdateNL _0/1_\ndispUpdateNL _distance_\nangleUpdateNL _angleDegree_";
  comp["numericalDampingCoeff"] = "numericalDampingCoeff _value_";
  comp["Velocity barriers"] =
      "VelocityBarrier _value_\n"
      "AngularVelocityBarrier _value_\n"
      "VelocityBarrierExponent _value_\n"
      "AngularVelocityBarrierExponent _value_\n";

  comp["Tempo"] = "Tempo NDCoeff _command_ _t1_ _t2_ _val1_ _val2_";

  comp["force law"] = "forceLaw _Default/Avalanche/StickedLinks/GeoVisc/BCM_";

  comp["Contact Parameters"] =
      "knContact _group1_ _group2_ _value_\n"
      "en2Contact _group1_ _group2_ _value_\n"
      "ktContact _group1_ _group2_ _value_\n"
      "muContact _group1_ _group2_ _value_\n"
      "krContact _group1_ _group2_ _value_\n"
      "murContact _group1_ _group2_ _value_\n";

  comp["InnerBond Parameters"] =
      "knInnerBond _group1_ _group2_ _value_\n"
      "en2InnerBond _group1_ _group2_ _value_\n"
      "ktInnerBond _group1_ _group2_ _value_\n"
      "fn0InnerBond _group1_ _group2_ _value_\n"
      "ft0InnerBond _group1_ _group2_ _value_\n"
      "powInnerBond _group1_ _group2_ _value_\n";

  comp["OuterBond Parameters"] =
      "knOuterBond _group1_ _group2_ _value_\n"
      "en2OuterBond _group1_ _group2_ _value_\n"
      "ktOuterBond _group1_ _group2_ _value_\n"
      "krOuterBond _group1_ _group2_ _value_\n"
      "fn0OuterBond _group1_ _group2_ _value_\n"
      "ft0OuterBond _group1_ _group2_ _value_\n"
      "mom0OuterBond _group1_ _group2_ _value_\n"
      "powOuterBond _group1_ _group2_ _value_\n";

  comp["BCM bond fracture energy"] =
      "gcInnerBond _group1_ _group2_ _value_\n"
      "gcOuterBond _group1_ _group2_ _value_\n";

  comp["Periodic cell"] =
      "usePeriodicCell 1\n"
      "h _xx_ _xy_ _xz_ _yx_ _yy_ _yz_ _zx_ _zy_ _zz_\n"
      "vh 0 0 0 0 0 0 0 0 0\n"
      "ah 0 0 0 0 0 0 0 0 0\n"
      "mh _massRatio_\n"
      "dh _damping_\n";

  comp["Soft particles"] = "useSoftParticles _Young_ _Poisson_";

  comp["Spring joint"] = "initSpringJoint _ibody_ _ix_ _iy_ _iz_ _jbody_ _jx_ _jy_ _jz_ _stiffness_";

  comp["Prevent crossing"] = "preventCrossingLength _length_";
}

void init_documentation() {

  docu["Rockable"] =
      "A conf-file always starts with the header: Rockable dd-mm-yyyy (e.g., Rockable 20-02-2017).\n"
      "Each time a noticeable change is made in the format, the date of this change is also changed\n"
      "in the header of the file. It has to be understood as the version of the format.";

  docu["t"] =
      "t [(double) value]\n\n"
      "Current time";
  docu["tmax"] =
      "tmax [(double) value]\n\n"
      "Maximum time (at the end of a simulation)";
  docu["dt"] =
      "dt [(double) value]\n\n"
      "Time step increment";
  docu["interVerlet"] =
      "interVerlet [(double) duration]\n\n"
      "Elapsed time between each update of the neighbor list (NL)";
  docu["interConf"] =
      "interConf [(double) duration]\n\n"
      "Elapsed time between each backup of the configuration";
  docu["DVerlet"] =
      "DVerlet [(double) distance]\n\n"
      "Increment size of the OBBs surrounding the rigid bodies";
  docu["dVerlet"] =
      "dVerlet [(double) distance]\n\n"
      "Increment size of the OBBs surrounding the sub-elements of the rigid bodies";
  docu["density"] =
      "density [(int)groupNumber] [(double) density]\n\n"
      "Set the density of all particles belonging to the group groupNumber";
  docu["gravity"] = "gravity [(vec3r) value] \n\n value = components of the gravity\n";
  docu["ParamsInInterfaces"] =
      "ParamsInInterfaces [0/1]\n\n"
      "States wether the parameters are embedded in the glued interfaces";
  docu["dynamicUpdateNL"] =
      "dynamicUpdateNL [0/1]\n\n"
      "If the dynamic update of the neighbor list is activated (set to 1), then an update will be made\n"
      "if the maximum distance of a body, since the last update, becomes larger than 'dispUpdateNL'.\n"
      "An update will also be made when the maximum rotation becomes larger than 'angleUpdateNL'";
  docu["dispUpdateNL"] =
      "dispUpdateNL [(double) distance]\n\n"
      "The maximum distance that any body is allowed to move without updating the NL";
  docu["angleUpdateNL"] =
      "angleUpdateNL [(double) angleDegree]\n\n"
      "The maximum angle that any body is allowed to rotate without updating the NL";
  docu["numericalDampingCoeff"] =
      "numericalDampingCoeff [(double) value]\n\n"
      "Cundall non-viscous numerical damping coefficient (typically between 0 and 1).\n"
      "At each step the accelerations are scaled by (1 - value) or (1 + value) depending\n"
      "on whether they oppose or follow the velocity, which damps the motion.\n"
      "Set to 0 to disable.";
  docu["VelocityBarrier"] =
      "VelocityBarrier [(double) value]\n\n"
      "Reference translational velocity above which a barrier force progressively slows\n"
      "the bodies down. The opposing force grows like (|v|/VelocityBarrier)^VelocityBarrierExponent.\n"
      "Set to 0 (default) to disable.";
  docu["AngularVelocityBarrier"] =
      "AngularVelocityBarrier [(double) value]\n\n"
      "Same as VelocityBarrier but applied to the angular velocity of the bodies.\n"
      "Set to 0 (default) to disable.";
  docu["VelocityBarrierExponent"] =
      "VelocityBarrierExponent [(double) value]\n\n"
      "Exponent used by the translational velocity barrier (see VelocityBarrier).\n"
      "The larger the exponent, the sharper the barrier.";
  docu["AngularVelocityBarrierExponent"] =
      "AngularVelocityBarrierExponent [(double) value]\n\n"
      "Exponent used by the angular velocity barrier (see AngularVelocityBarrier).";
  docu["Tempo"] =
      "Tempo [(string) target] ...\n\n"
      "Defines a time-evolving ('tempo') parameter interpolated between two times.\n"
      "Three forms are available:\n\n"
      "  Tempo NDCoeff [command] [t1] [t2] [val1] [val2]\n"
      "      drives the numericalDampingCoeff\n"
      "  Tempo Inter [paramName] [group1] [group2] [command] [t1] [t2] [val1] [val2]\n"
      "      drives an interaction parameter (e.g. knContact) between two groups\n"
      "  Tempo Body [paramName] [group] [command] [t1] [t2] [val1] [val2]\n"
      "      drives a body property (e.g. density) of a group\n\n"
      "'command' selects the interpolation profile between (t1,val1) and (t2,val2).";
  docu["forceLaw"] =
      "forceLaw [(string)Name]\n\n"
      "Select a model for the computation of forces.\n"
      "Possible Names are: 'Default', 'Avalanche', 'StickedLinks', 'GeoVisc' or 'BCM'.";

  docu["AddOrRemoveInteractions"] =
      "AddOrRemoveInteractions [(string) Name]\n\n"
      "Strategy used to detect which interactions must be created or deleted.\n"
      "Possible Names are: 'bruteForce' (default) or 'OBBtree'.";
  docu["UpdateNL"] =
      "UpdateNL [(string) Name]\n\n"
      "Strategy used to rebuild the neighbor list (NL).\n"
      "Possible Names are: 'bruteForce' (default) or 'linkCells'.";
  docu["Integrator"] =
      "Integrator [(string) Name]\n\n"
      "Time-integration scheme.\n"
      "Possible Names are: 'velocityVerlet' (default), 'Euler', 'Beeman' or 'RungeKutta4'.";
  docu["cellMinSizes"] =
      "cellMinSizes [(double) value]\n\n"
      "Minimum size of the cells used by the 'linkCells' neighbor-list strategy.";
  docu["boxForLinkCellsOpt"] =
      "boxForLinkCellsOpt [(int) 0/1]\n\n"
      "Selects the bounding box used to build the link-cells grid.\n"
      "0 (default): box surrounding all particles;\n"
      "1: box surrounding only the free (non-driven) bodies.";
  docu["ContactPartnership"] =
      "ContactPartnership [(string) Name]\n\n"
      "Model used to weight the contact stiffnesses when several sub-elements of two\n"
      "bodies interact at once (used by some force laws).\n"
      "Possible Names are: 'None' (default), 'NumberWeight', 'OverlapWeight' or 'SurfaceWeight'.";

  docu["preventCrossingLength"] =
      "preventCrossingLength [(double) value]\n\n"
      "If positive, enables a safeguard that prevents bodies from crossing/tunneling\n"
      "through each other when the overlap would exceed this characteristic length.";
  docu["initSpringJoint"] =
      "initSpringJoint [(int) ibody] [(vec3r) ipos0] [(int) jbody] [(vec3r) jpos0] [(double) stiffness]\n\n"
      "Create a linear spring joint between bodies ibody and jbody, anchored at the\n"
      "local positions ipos0 and jpos0, with the given stiffness.";
  docu["PostProcessor"] =
      "PostProcessor [(string) Name] ...\n\n"
      "Adds a post-processor (used by the 'postpro' tool).\n"
      "Possible Names are: 'ClusterGranulo' or 'ParticleStress'.\n"
      "The remaining parameters depend on the chosen post-processor.";

  // periodic cell
  docu["usePeriodicCell"] =
      "usePeriodicCell [0/1]\n\n"
      "Enable (1) or disable (0) the periodic boundary cell.\n"
      "Requires Rockable compiled with ROCKABLE_ENABLE_PERIODIC.";
  docu["h"] =
      "h [(mat9r) matrix]\n\n"
      "The 3x3 matrix (9 components) defining the periodic cell (column vectors are the\n"
      "periodicity vectors).";
  docu["vh"] = "vh [(mat9r) matrix]\n\nTime derivative (velocity) of the periodic cell matrix 'h'.";
  docu["ah"] = "ah [(mat9r) matrix]\n\nSecond time derivative (acceleration) of the periodic cell matrix 'h'.";
  docu["mh"] = "mh [(double) value]\n\nMass ratio assigned to the periodic cell degrees of freedom.";
  docu["dh"] = "dh [(double) value]\n\nNumerical damping applied to the periodic cell.";
  docu["cellMomentumCorrection"] = "cellMomentumCorrection [0/1]\n\nEnable correction of the total momentum within the periodic cell.";
  docu["cellVelocityCorrection"] = "cellVelocityCorrection [0/1]\n\nEnable correction of the velocities within the periodic cell.";
  docu["useKineticStress"] = "useKineticStress [0/1]\n\nInclude the kinetic (velocity fluctuation) contribution in the stress of the periodic cell.";

  // soft particles
  docu["useSoftParticles"] =
      "useSoftParticles [(double) Young] [(double) Poisson]\n\n"
      "Enable deformable (soft) particles with the given Young modulus and Poisson ratio.\n"
      "Requires Rockable compiled with ROCKABLE_ENABLE_SOFT_PARTICLES.";

  docu["knContact"] = "knContact [(int) group1] [(int) group2] [(double) value]";
  docu["en2Contact"] = "en2Contact [(int) group1] [(int) group2] [(double) value]";
  docu["en2ContactFromViscRate"] = "en2ContactFromViscRate [(int) group1] [(int) group2] [(double) value]";
  docu["ktContact"] = "ktContact [(int) group1] [(int) group2] [(double) value]";
  docu["muContact"] = "muContact [(int) group1] [(int) group2] [(double) value]";
  docu["krContact"] = "krContact [(int) group1] [(int) group2] [(double) value]";
  docu["murContact"] = "murContact [(int) group1] [(int) group2] [(double) value]";

  docu["knInnerBond"] = "knInnerBond [(int) group1] [(int) group2] [(double) value]";
  docu["en2InnerBond"] = "en2InnerBond [(int) group1] [(int) group2] [(double) value]";
  docu["ktInnerBond"] = "ktInnerBond [(int) group1] [(int) group2] [(double) value]";
  docu["fn0InnerBond"] = "fn0InnerBond [(int) group1] [(int) group2] [(double) value]";
  docu["ft0InnerBond"] = "ft0InnerBond [(int) group1] [(int) group2] [(double) value]";
  docu["powInnerBond"] = "powInnerBond [(int) group1] [(int) group2] [(double) value]";

  docu["knOuterBond"] = "knOuterBond [(int) group1] [(int) group2] [(double) value]";
  docu["ktOuterBond"] = "ktOuterBond [(int) group1] [(int) group2] [(double) value]";
  docu["krOuterBond"] = "krOuterBond [(int) group1] [(int) group2] [(double) value]";
  docu["fn0OuterBond"] = "fn0OuterBond [(int) group1] [(int) group2] [(double) value]";
  docu["ft0OuterBond"] = "ft0OuterBond [(int) group1] [(int) group2] [(double) value]";
  docu["mom0OuterBond"] = "mom0OuterBond [(int) group1] [(int) group2] [(double) value]";
  docu["powOuterBond"] = "powOuterBond [(int) group1] [(int) group2] [(double) value]";
  docu["en2OuterBond"] = "en2OuterBond [(int) group1] [(int) group2] [(double) value]";

  docu["krInnerBond"] = "krInnerBond [(int) group1] [(int) group2] [(double) value]";
  docu["mom0InnerBond"] = "mom0InnerBond [(int) group1] [(int) group2] [(double) value]";
  docu["gcInnerBond"] =
      "gcInnerBond [(int) group1] [(int) group2] [(double) value]\n\n"
      "Critical fracture energy (energy release rate Gc) of the inner bonds,\n"
      "used by the BCM (Bonded Cell Model) force law.";
  docu["gcOuterBond"] =
      "gcOuterBond [(int) group1] [(int) group2] [(double) value]\n\n"
      "Critical fracture energy (energy release rate Gc) of the outer bonds,\n"
      "used by the BCM (Bonded Cell Model) force law.";

  docu["iconf"] =
      "iconf [(int) value]"
      "Number of the current configuration";
  docu["nDriven"] =
      "nDriven [(int) value]\n\n"
      "The value is the number of bodies, at the beginning of the list, that nDriven first bodies\n"
      "are fixed (all velocities imposed to zero), but if we want to set a velocity or a force/moment,\n"
      "some commands have to be added in a file named drivingSystem.txt";
  docu["shapeFile"] =
      "shapeFile [(string) path]\n\n"
      "Paths of the file that defines the shapes used.";

  docu["separator"] =
      "By default, the separator character in a conf-file is the 'space'\n"
      "but it can be changed to 'tab' or 'semicolon'";
  docu["glue_with_walls"] =
      "glue_with_walls [0/1]\n\n"
      "State wether the glue will be activated between driven and not-driven bodies\n"
      "when the command 'stickClusters' will be invoked";
  docu["precision"] =
      "precision [(int) value]\n\n"
      "Set the number of digits to be used in the conf-files";

  docu["Particles"] =
      "Particles [(int) numberOfParticles]\n\n"
      "REPEAT FOR EACH PARTICLE {\n"
      "   [(string) shapeName] [(int) group] [(int) cluster] [(double) homothety] \n"
      "   [(vec3r) position] [(vec3r) velocity] [(vec3r) acceleration] \n"
      "   [(quat) angularPosition] [(vec3r) angularVelocity] [(vec3r) angularAcceleration]\n"
      "}";

  docu["Interactions"] =
      "Interactions [(int) numberOfInteractions]\n\n"
      "REPEAT FOR EACH INTERACTION {\n"
      "   [(int) i] [(int) j] [(int) type] [(int) isub] [(int) jsub] [(vec3r) n] [(double) dn] [(vec3r) position]\n"
      "   [(vec3r) relativeVelocity] [(double) fn] [(vec3r) ft] [(vec3r) mom] [(double) viscousDampingValue]\n"
      "}";
  docu["Interfaces"] =
      "Interfaces [(int)numberOfInterfaces]\n\n"
      "REPEAT FOR EACH INTERFACE {\n"
      "   [(int) i] [(int) j] [(int) nbBonds]\n"
      "   IF (ParamsInInterfaces is not 0) {\n"
      "      [(double) kn] [(double) kt] [(double) kr] [(double) fn0] [(double) ft0] [(double) mom0] [(double) power]\n"
      "   }\n"
      "   REPEAT FOR EACH BOND { [(int) type] [(int) isub] [(int) jsub] }\n"
      "}";
  docu["DataExtractor"] =
      "DataExtractor [(string) Name] ...\n\n"
      "Adds a data extractor that records quantities during the simulation.\n"
      "Possible Names are: 'ClusterAABB', 'dnStat', 'DuoBalance', 'MeanVelocity' or 'TrackBody'.\n"
      "The remaining parameters depend on the chosen extractor.\n"
      "Kept for compatibility: prefer defining them in a file named 'dataExtractors.txt'.";
  docu["BodyForce"] =
      "BodyForce [(string) Name] ...\n\n"
      "Activates a body force applied to the particles.\n"
      "Possible Names are: 'AttractingPoint', 'PreferredDirection' or 'ViscousFluid'.\n"
      "The remaining parameters depend on the chosen body force.";

  docu["stickVerticesInClusters"] =
      "stickVerticesInClusters [(double) Epsilon]\n\n"
      "This command will add glued interfaces between bodies having\n"
      "the same cluster identifier. Only bonds between vertices (spheres)\n"
      "are added when the distance is less than Epsilon.";
  docu["stickClusters"] =
      "stickClusters [(double) Epsilon]\n\n"
      "This command will add glued interfaces between bodies having different cluster identifiers.\n"
      "Bonds are added when the distance is less than Epsilon.";
  docu["stickVerticesInClustersMoments"] =
      "stickVerticesInClustersMoments [(double) Epsilon]\n\n"
      "Like stickVerticesInClusters, but the created bonds also transmit moments.\n"
      "Bonds between vertices (spheres) of a same cluster are added when the distance\n"
      "is less than Epsilon.";
  docu["stickBCM"] =
      "stickBCM [(double) Epsilon]\n\n"
      "Create the bonded interfaces required by the BCM (Bonded Cell Model) force law\n"
      "between neighboring sub-elements closer than Epsilon.";

  docu["copyParamsToInterfaces"] = "copyParamsToInterfaces [(string) isInnerStr]";
  docu["setStiffnessRatioInterfaces"] =
      "setStiffnessRatioInterfaces [(double) ratio]\n\n"
      "Set the tangential-to-normal stiffness ratio kt/kn for every bond stored in the\n"
      "glued interfaces (only for parameters that are embedded in the Interfaces).";
  docu["setVariableStickParams"] =
      "setVariableStickParams [(string) paramName] [(string) isInnerStr] [(double) lambda] [(double) m] [(int) "
      "timeSeeded]";
  docu["randomlyOrientedVelocities"] = "randomlyOrientedVelocities [(double) velocity]";
  docu["randomlyOrientedVelocitiesClusters"] = "randomlyOrientedVelocitiesClusters [(double) velocity] [(int) option]";
  docu["setAllVelocities"] =
      "setAllVelocities [(vec3r) velocity]\n\n"
      "Set the same velocity vector to ALL free (non-driven) bodies.";
  docu["homothetyRange"] =
      "homothetyRange [(int) idFirst] [(int) idLast] [(double) hmin] [(double) hmax] [(int) timeSeeded]\n\n"
      "Assign to the bodies idFirst..idLast a homothety drawn from a uniform distribution\n"
      "in [hmin, hmax]. If timeSeeded is 1, the random generator is seeded with the current time.";
  docu["particlesClonage"] =
      "particlesClonage [(int) idFirst] [(int) idLast] [(vec3r) translation]\n\n"
      "Duplicate the bodies idFirst..idLast, shifting the clones by the translation vector.\n"
      "The clones receive new cluster identifiers.";

  // driving system
  docu["Control"] =
      "Control [(string) mode] [(int) bodyNumber] [(double) value]\n\n"
      "where mode is either:\n"
      "                     _x_Vel_, _y_Vel_, _z_Vel_,\n"
      "                     _xrot_Vel_, _yrot_Vel_, _zrot_Vel_,\n"
      "                     _x_For_, _y_For_, _z_For_,\n"
      "                     _xrot_Mom_, _yrot_Mom_, _zrot_Mom_,\n"
      "                     _xyzrot_Vel_, or _xyzrot_Mom_ \n\n"
      "This is a simple solution to impose a constant velocity or force/moment to a body.\n"
      "Notice that bodyNumber needs to be lower than 'nDriven'";
  docu["Servo"] =
      "Servo [(string) Name] ...\n\n"
      "Activates a predefined servo-controller that drives a set of walls/bodies.\n"
      "The 'tritri*' servos (e.g. tritriIsostaticCompression, tritriBiaxialCompression,\n"
      "tritriCustom, tritriLodeAngle) expect the 6 wall body numbers\n"
      "(idXmin idXmax idYmin idYmax idZmin idZmax) followed by servo-specific parameters.\n"
      "Other servos: shaker, triangle_shaker, sawtooth_shaker, ramp.\n"
      "Servo commands are defined in the file 'drivingSystem.txt'.";
  docu["PeriodicLoading"] =
      "PeriodicLoading [(string) Name] ...\n\n"
      "Drive the periodic cell with a predefined loading.\n"
      "Example: 'PeriodicLoading IsotropicCompression [(double) pressure]'.\n"
      "Defined in the file 'drivingSystem.txt'.";

  // post-processing (postpro tool)
  docu["PostProcessor"] =
      "PostProcessor [(string) Name] ...\n\n"
      "Adds a post-processor (used by the 'postpro' tool).\n"
      "Possible Names are: 'ClusterGranulo' or 'ParticleStress'.";
  docu["firstConf"] = "firstConf [(int) value]\n\nNumber of the first conf-file to post-process.";
  docu["lastConf"] = "lastConf [(int) value]\n\nNumber of the last conf-file to post-process.";
  docu["stepConf"] = "stepConf [(int) value]\n\nStep between two processed conf-files.";
  docu["SievingSizes"] =
      "SievingSizes [(int) nb] [(double) size1] [(double) size2] ...\n\n"
      "List of nb sieving sizes used by the ClusterGranulo post-processor.";
  docu["ConfVolumes"] =
      "ConfVolumes [(int) nb]\n"
      "REPEAT nb TIMES { [(int) iconf] [(double) volume] }\n\n"
      "Per-conf reference volumes used by the ParticleStress post-processor.";
  docu["Volume"] = "Volume [(double) value]\n\nReference volume used by the ParticleStress post-processor.";

  // shape
  docu["name"] =
      "name [(string) shapeName]\n\n"
      "Set the name of the shape";
  docu["radius"] =
      "radius [(double) value]\n\n"
      "set the radius of the rounded edges of the shape";
  docu["preCompDone"] =
      "preCompDone [y/n]\n\n"
      "If yes, the computation of mass properties will not be run\n"
      "at the beginning of a computation in Rockable (or in the program shapeSurvey)";
  docu["volume"] = "volume [(double) V]";
  docu["I/m"] = "I/m [(double) I1/m] [(double) I2/m] [(double) I3/m]";
  docu["obb.extent"] =
      "obb.extent [(double)extent1] [(double)extent2] [(double)extent3]\n\n"
      "The 3 components of the extents in the 3 directions defined by obb.e1, obb.e2 and obb.e3";
  docu["obb.center"] =
      "obb.center [(vec3r) center]\n\n"
      "Position of the OBB center relative to the mass center of the shape,\n"
      "given in the framework (principal direction of inertia) of the shape.";
  docu["obb.e1"] = "obb.e1 [(vec3r) e1]\n\nExtent unit vector in the first direction";
  docu["obb.e2"] = "obb.e2 [(vec3r) e2]\n\nExtent unit vector in the second direction";
  docu["obb.e3"] = "obb.e3 [(vec3r) e3]\n\nExtent unit vector in the third direction";
  docu["fibObbOption"] =
      "fibObbOption [(int) option]\n\n"
      "Algorithm used to build the OBB of the shape:\n"
      "  0: covariance-based orientation\n"
      "  1: minimum-volume strategy\n"
      "  2: OBB aligned with the axes (AABB) [default]\n"
      "  3: imposed axis";
  docu["OBBtreeLevel"] = "DEPRECATED!!";
  docu["position"] =
      "position [(vec3r) center]\n\n"
      "Position of the mass center of the shape, expressed in the local frame of the shape.\n"
      "Usually computed automatically (Monte-Carlo integration) when preCompDone is 'n'.";
  docu["orientation"] =
      "orientation [(quat) Q]\n\n"
      "Angular position (quaternion) of the principal-inertia frame of the shape.\n"
      "Usually computed automatically when preCompDone is 'n'.";
  docu["isSurface"] =
      "isSurface\n\n"
      "Flag (no value) marking the shape as an open surface rather than a closed volume.";
  docu["MCnstep"] =
      "MCnstep [(int) value]\n\n"
      "Number of steps used in the Monte-Carlo integration of the mass properties\n"
      "(volume, mass center, inertia). Clamped internally to the range [1000, 10000000].";
  docu["nv"] =
      "nv [(int) numberOfVertices]\n\n"
      "REPEAT FOR EACH VERTEX { [(vec3r) position] }\n"
      "List of the vertices (sphero-points) of the shape, in the local frame.";
  docu["ne"] =
      "ne [(int) numberOfEdges]\n\n"
      "REPEAT FOR EACH EDGE { [(int) fromVertex] [(int) toVertex] }\n"
      "List of the edges (sphero-segments) connecting two vertices.";
  docu["nf"] =
      "nf [(int) numberOfFaces]\n\n"
      "REPEAT FOR EACH FACE { [(int) nbNodes] [(int) node1] [(int) node2] ... }\n"
      "List of the faces (sphero-polygons), each given by its number of nodes\n"
      "followed by the corresponding vertex indices.";
}

void init_keywords() {
  // code_keywords.insert(std::string("____"));

  code_keywords.insert(std::string("Rockable"));
  code_keywords.insert(std::string("t"));
  code_keywords.insert(std::string("tmax"));
  code_keywords.insert(std::string("dt"));
  code_keywords.insert(std::string("interVerlet"));
  code_keywords.insert(std::string("interConf"));
  code_keywords.insert(std::string("DVerlet"));
  code_keywords.insert(std::string("dVerlet"));
  code_keywords.insert(std::string("density"));
  code_keywords.insert(std::string("gravity"));
  code_keywords.insert(std::string("ParamsInInterfaces"));
  code_keywords.insert(std::string("dynamicUpdateNL"));
  code_keywords.insert(std::string("dispUpdateNL"));
  code_keywords.insert(std::string("angleUpdateNL"));
  code_keywords.insert(std::string("numericalDampingCoeff"));
  code_keywords.insert(std::string("VelocityBarrier"));
  code_keywords.insert(std::string("AngularVelocityBarrier"));
  code_keywords.insert(std::string("VelocityBarrierExponent"));
  code_keywords.insert(std::string("AngularVelocityBarrierExponent"));
  code_keywords.insert(std::string("Tempo"));
  code_keywords.insert(std::string("forceLaw"));
  code_keywords.insert(std::string("AddOrRemoveInteractions"));
  code_keywords.insert(std::string("UpdateNL"));
  code_keywords.insert(std::string("Integrator"));
  code_keywords.insert(std::string("cellMinSizes"));
  code_keywords.insert(std::string("boxForLinkCellsOpt"));
  code_keywords.insert(std::string("ContactPartnership"));
  code_keywords.insert(std::string("preventCrossingLength"));
  code_keywords.insert(std::string("initSpringJoint"));
  code_keywords.insert(std::string("PostProcessor"));

  // periodic cell (requires Rockable compiled with ROCKABLE_ENABLE_PERIODIC)
  code_keywords.insert(std::string("usePeriodicCell"));
  code_keywords.insert(std::string("h"));
  code_keywords.insert(std::string("vh"));
  code_keywords.insert(std::string("ah"));
  code_keywords.insert(std::string("mh"));
  code_keywords.insert(std::string("dh"));
  code_keywords.insert(std::string("cellMomentumCorrection"));
  code_keywords.insert(std::string("cellVelocityCorrection"));
  code_keywords.insert(std::string("useKineticStress"));

  // soft particles (requires Rockable compiled with ROCKABLE_ENABLE_SOFT_PARTICLES)
  code_keywords.insert(std::string("useSoftParticles"));

  code_keywords.insert(std::string("knContact"));
  code_keywords.insert(std::string("en2Contact"));
  code_keywords.insert(std::string("en2ContactFromViscRate"));
  code_keywords.insert(std::string("ktContact"));
  code_keywords.insert(std::string("muContact"));
  code_keywords.insert(std::string("krContact"));
  code_keywords.insert(std::string("murContact"));
  code_keywords.insert(std::string("viscTContact"));

  code_keywords.insert(std::string("knInnerBond"));
  code_keywords.insert(std::string("en2InnerBond"));
  code_keywords.insert(std::string("ktInnerBond"));
  code_keywords.insert(std::string("krInnerBond"));
  code_keywords.insert(std::string("fn0InnerBond"));
  code_keywords.insert(std::string("ft0InnerBond"));
  code_keywords.insert(std::string("mom0InnerBond"));
  code_keywords.insert(std::string("powInnerBond"));
  code_keywords.insert(std::string("gcInnerBond"));

  code_keywords.insert(std::string("knOuterBond"));
  code_keywords.insert(std::string("ktOuterBond"));
  code_keywords.insert(std::string("krOuterBond"));
  code_keywords.insert(std::string("fn0OuterBond"));
  code_keywords.insert(std::string("ft0OuterBond"));
  code_keywords.insert(std::string("mom0OuterBond"));
  code_keywords.insert(std::string("powOuterBond"));
  code_keywords.insert(std::string("en2OuterBond"));
  code_keywords.insert(std::string("gcOuterBond"));

  code_keywords.insert(std::string("iconf"));
  code_keywords.insert(std::string("nDriven"));
  code_keywords.insert(std::string("shapeFile"));

  code_keywords.insert(std::string("separator"));
  code_keywords.insert(std::string("glue_with_walls"));
  code_keywords.insert(std::string("precision"));

  code_keywords.insert(std::string("Particles"));
  code_keywords.insert(std::string("Interactions"));
  code_keywords.insert(std::string("Interfaces"));
  code_keywords.insert(std::string("DataExtractor"));
  code_keywords.insert(std::string("BodyForce"));

  // processing
  code_keywords.insert(std::string("stickVerticesInClusters"));
  code_keywords.insert(std::string("stickVerticesInClustersMoments"));
  code_keywords.insert(std::string("stickClusters"));
  code_keywords.insert(std::string("stickBCM"));
  code_keywords.insert(std::string("copyParamsToInterfaces"));
  code_keywords.insert(std::string("setStiffnessRatioInterfaces"));
  code_keywords.insert(std::string("setVariableStickParams"));
  code_keywords.insert(std::string("randomlyOrientedVelocities"));
  code_keywords.insert(std::string("randomlyOrientedVelocitiesClusters"));
  code_keywords.insert(std::string("setAllVelocities"));
  code_keywords.insert(std::string("homothetyRange"));
  code_keywords.insert(std::string("particlesClonage"));

  // driving system
  code_keywords.insert(std::string("Control"));
  code_keywords.insert(std::string("Servo"));
  code_keywords.insert(std::string("PeriodicLoading"));

  // post-processing (postpro tool input files)
  code_keywords.insert(std::string("PostProcessor"));
  code_keywords.insert(std::string("firstConf"));
  code_keywords.insert(std::string("lastConf"));
  code_keywords.insert(std::string("stepConf"));
  code_keywords.insert(std::string("ConfVolumes"));
  code_keywords.insert(std::string("SievingSizes"));
  code_keywords.insert(std::string("Volume"));

  // shapes
  code_keywords.insert(std::string("name"));
  code_keywords.insert(std::string("radius"));
  code_keywords.insert(std::string("preCompDone"));
  code_keywords.insert(std::string("volume"));
  code_keywords.insert(std::string("I/m"));
  code_keywords.insert(std::string("obb.extent"));
  code_keywords.insert(std::string("obb.center"));
  code_keywords.insert(std::string("obb.e1"));
  code_keywords.insert(std::string("obb.e2"));
  code_keywords.insert(std::string("obb.e3"));
  code_keywords.insert(std::string("fibObbOption"));
  code_keywords.insert(std::string("OBBtreeLevel"));
  code_keywords.insert(std::string("position"));
  code_keywords.insert(std::string("orientation"));
  code_keywords.insert(std::string("isSurface"));
  code_keywords.insert(std::string("MCnstep"));
  code_keywords.insert(std::string("nv"));
  code_keywords.insert(std::string("ne"));
  code_keywords.insert(std::string("nf"));
}

void init_types() {
  // code_types.insert(std::string("____"));

  // yes/no
  code_types.insert(std::string("y"));
  code_types.insert(std::string("n"));

  // tempo
  code_types.insert(std::string("NDCoeff"));
  code_types.insert(std::string("Inter"));
  code_types.insert(std::string("Body"));

  // force laws
  code_types.insert(std::string("Default"));
  code_types.insert(std::string("Avalanche"));
  code_types.insert(std::string("StickedLinks"));
  code_types.insert(std::string("GeoVisc"));
  code_types.insert(std::string("BCM"));

  // add or remove interaction
  code_types.insert(std::string("bruteForce"));
  code_types.insert(std::string("OBBtree"));

  // updateNL
  code_types.insert(std::string("linkCells"));

  // contact partnership models
  code_types.insert(std::string("None"));
  code_types.insert(std::string("NumberWeight"));
  code_types.insert(std::string("OverlapWeight"));
  code_types.insert(std::string("SurfaceWeight"));

  // integrator
  code_types.insert(std::string("velocityVerlet"));
  code_types.insert(std::string("Euler"));
  code_types.insert(std::string("Beeman"));
  code_types.insert(std::string("RungeKutta4"));

  // body forces
  code_types.insert(std::string("AttractingPoint"));
  code_types.insert(std::string("PreferredDirection"));
  code_types.insert(std::string("ViscousFluid"));

  // data extractor
  code_types.insert(std::string("ClusterAABB"));
  code_types.insert(std::string("dnStat"));
  code_types.insert(std::string("DuoBalance"));
  code_types.insert(std::string("MeanVelocity"));
  code_types.insert(std::string("TrackBody"));
  code_types.insert(std::string("TrackDamage"));
  code_types.insert(std::string("TrackRockfall"));

  // post-processors
  code_types.insert(std::string("ClusterGranulo"));
  code_types.insert(std::string("ParticleStress"));

  // controls (driving system)
  code_types.insert(std::string("_x_Vel_"));
  code_types.insert(std::string("_y_Vel_"));
  code_types.insert(std::string("_z_Vel_"));
  code_types.insert(std::string("_xrot_Vel_"));
  code_types.insert(std::string("_yrot_Vel_"));
  code_types.insert(std::string("_zrot_Vel_"));
  code_types.insert(std::string("_x_For_"));
  code_types.insert(std::string("_y_For_"));
  code_types.insert(std::string("_z_For_"));
  code_types.insert(std::string("_xrot_Mom_"));
  code_types.insert(std::string("_yrot_Mom_"));
  code_types.insert(std::string("_zrot_Mom_"));
  code_types.insert(std::string("_xyzrot_Vel_"));
  code_types.insert(std::string("_xyzrot_Mom_"));

  // servo names (driving system)
  code_types.insert(std::string("tritriIsostaticCompression"));
  code_types.insert(std::string("tritriBiaxialCompression"));
  code_types.insert(std::string("tritriCustom"));
  code_types.insert(std::string("tritriLodeAngle"));
  code_types.insert(std::string("shaker"));
  code_types.insert(std::string("triangle_shaker"));
  code_types.insert(std::string("sawtooth_shaker"));
  code_types.insert(std::string("ramp"));
}

// Returns the stream's current offset, clamped to [0, length]. std::istream
// returns -1 from tellg() once EOF has been reached (e.g. after extracting the
// last token of a buffer without a trailing separator); without this clamp the
// computed offsets go negative and we write style[] out of bounds.
static int stream_pos(std::istream& is, int length) {
  std::streampos p = is.tellg();
  if (p == std::streampos(-1)) return length;
  int pos = static_cast<int>(p);
  if (pos < 0) return 0;
  if (pos > length) return length;
  return pos;
}

void style_parse(const char* text, char* style, int length) {
  for (int i = 0; i < length; i++) {
    if (text[i] != '\n') {
      style[i] = 'A';
    }
  }
  std::istringstream is(std::string(text, length));

  std::string token;
  while (is >> token) {
    int to = stream_pos(is, length);
    int from = to - static_cast<int>(token.size());
    if (from < 0) from = 0;

    if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
      getline(is, token);  // ignore the rest of the current line
      int lineEnd = stream_pos(is, length);
      for (int pos = from; pos < lineEnd; pos++) {
        style[pos] = 'B';
      }
      continue;
    }

    if (code_keywords.find(token) != code_keywords.end()) {
      for (int pos = from; pos < to; pos++) {
        style[pos] = 'D';
      }
    } else if (code_types.find(token) != code_types.end()) {
      for (int pos = from; pos < to; pos++) {
        style[pos] = 'C';
      }
    }
  }  // while
}

void style_init(void) {
  char* style = new char[textbuf->length() + 1];
  char* text = textbuf->text();

  memset(style, 'A', textbuf->length());
  style[textbuf->length()] = '\0';

  if (!stylebuf) {
    stylebuf = new Fl_Text_Buffer(textbuf->length());
  }

  style_parse(text, style, textbuf->length());

  stylebuf->text(style);
  delete[] style;
  free(text);
}

void style_unfinished_cb(int, void*) {}

void style_update(int pos,                      // I - Position of update
                  int nInserted,                // I - Number of inserted chars
                  int nDeleted,                 // I - Number of deleted chars
                  int /*nRestyled*/,            // I - Number of restyled chars
                  const char* /*deletedText*/,  // I - Text that was deleted
                  void* cbArg) {                // I - Callback data
  int start,                                    // Start of text
      end;                                      // End of text
  char *style,                                  // Style data
      *text;                                    // Text data

  // If this is just a selection change, just unselect the style buffer...
  if (nInserted == 0 && nDeleted == 0) {
    stylebuf->unselect();
    return;
  }

  // Track changes in the text buffer...
  if (nInserted > 0) {
    // Insert characters into the style buffer...
    style = new char[nInserted + 1];
    memset(style, 'A', nInserted);
    style[nInserted] = '\0';

    stylebuf->replace(pos, pos + nDeleted, style);
    delete[] style;
  } else {
    // Just delete characters in the style buffer...
    stylebuf->remove(pos, pos + nDeleted);
  }

  // Select the area that was just updated to avoid unnecessary
  // callbacks...
  stylebuf->select(pos, pos + nInserted - nDeleted);

  // Re-parse only the line(s) touched by the change. The Rockable conf format
  // has no multi-line constructs (comments end at the line), so a line's style
  // depends solely on its own content; there is no need to reparse the rest of
  // the buffer, which would be very slow on large files.
  start = textbuf->line_start(pos);
  end = textbuf->line_end(pos + nInserted);
  text = textbuf->text_range(start, end);
  style = stylebuf->text_range(start, end);

  style_parse(text, style, end - start);

  stylebuf->replace(start, end, style);
  ((Fl_Text_Editor*)cbArg)->redisplay_range(start, end);

  free(text);
  free(style);
}

EditorWindow::EditorWindow(int w, int h, const char* t) : Fl_Double_Window(w, h, t) {

  // replace dialog
  replace_dlg = new Fl_Window(300, 105, "Replace");
  replace_find = new Fl_Input(80, 10, 210, 25, "Find:");
  replace_find->align(FL_ALIGN_LEFT);

  replace_with = new Fl_Input(80, 40, 210, 25, "Replace:");
  replace_with->align(FL_ALIGN_LEFT);

  replace_all = new Fl_Button(10, 70, 90, 25, "Replace All");
  replace_all->callback((Fl_Callback*)replall_cb, this);

  replace_next = new Fl_Return_Button(105, 70, 120, 25, "Replace Next");
  replace_next->callback((Fl_Callback*)replace2_cb, this);

  replace_cancel = new Fl_Button(230, 70, 60, 25, "Cancel");
  replace_cancel->callback((Fl_Callback*)replcan_cb, this);
  replace_dlg->end();
  replace_dlg->set_non_modal();

  // keyword dialog
  keyword_dlg = new Fl_Window(300, 400, "Keywords");

  keyword_browser = new Fl_Browser(10, 10, 280, 380 - 25 - 5);
  keyword_browser->type(FL_HOLD_BROWSER);

  for (std::map<std::string, std::string>::iterator entry = comp.begin(); entry != comp.end(); ++entry) {
    keyword_browser->add(entry->first.c_str());
  }

  keyword_insert = new Fl_Button(230 - 60 - 10, 365, 60, 25, "Insert");
  keyword_insert->callback((Fl_Callback*)keyword_insert_cb, this);

  keyword_close = new Fl_Button(230, 365, 60, 25, "Close");
  keyword_close->callback((Fl_Callback*)keyword_close_cb, this);

  keyword_dlg->end();
  keyword_dlg->set_non_modal();

  // add particle dialog
  add_particle_dlg = new addParticleDialog(400, 150, "Add a particle");

  //  ========
  editor = 0;
  *search = (char)0;
  wrap_mode = 0;
  line_numbers = 0;
}

EditorWindow::~EditorWindow() {
  delete replace_dlg;
  delete keyword_dlg;
  delete add_particle_dlg;
}

void doc_selection_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Buffer* buf = e->editor->buffer();
  if (buf->selected()) {
    std::map<std::string, std::string>::iterator it = docu.find(buf->selection_text());
    if (it != docu.end()) {
      fl_message("%s", it->second.c_str());
    }
  }
}

void deg_to_rad_cb(Fl_Widget*, void* v) {
  const char* deg = fl_input("Value in degree: ", "0");
  char radtxt[100];
  snprintf(radtxt, 100, "%f", atof(deg) * M_PI / 180.0);

  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Buffer* buf = e->editor->buffer();
  buf->insert(e->editor->insert_position(), radtxt);
}

void quat_axis_angle_cb(Fl_Widget*, void* v) {
  const char* input = fl_input("Enter axis vector followed by an angle (degree): ", "1 0 0 90");
  std::stringstream ss(input);
  double axis_x, axis_y, axis_z, angleDeg;
  ss >> axis_x >> axis_y >> axis_z >> angleDeg;
  double s, vx, vy, vz;
  double half_angle = 0.5 * angleDeg * M_PI / 180.;
  double len = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
  s = cos(half_angle);
  vx = sin(half_angle) * axis_x / len;
  vy = sin(half_angle) * axis_y / len;
  vz = sin(half_angle) * axis_z / len;
  char txt[100];
  snprintf(txt, 100, "%f  %f %f %f", s, vx, vy, vz);

  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Buffer* buf = e->editor->buffer();
  buf->insert(e->editor->insert_position(), txt);
}

void add_particle_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->add_particle_dlg->show();
}

int check_save(void) {
  if (!changed) return 1;

  int r = fl_choice(
      "The current file has not been saved.\n"
      "Would you like to save it now?",
      "Cancel", "Save", "Don't Save");

  if (r == 1) {
    save_cb();  // Save the file...
    return !changed;
  }

  return (r == 2) ? 1 : 0;
}

void load_file(const char* newfile, int ipos) {
  loading = 1;
  int insert = (ipos != -1);
  changed = insert;
  if (!insert) {
    filename[0] = '\0';
  }
  int r;
  if (!insert) {
    r = textbuf->loadfile(newfile);
  } else {
    r = textbuf->insertfile(newfile, ipos);
  }
  changed = changed || textbuf->input_file_was_transcoded;
  if (r) {
    fl_alert("Error reading from file \'%s\':\n%s.", newfile, strerror(errno));
  } else if (!insert) {
    snprintf(filename, sizeof(filename), "%s", newfile);
  }
  loading = 0;
  textbuf->call_modify_callbacks();
}

void save_file(const char* newfile) {
  if (textbuf->savefile(newfile)) {
    fl_alert("Error writing to file \'%s\':\n%s.", newfile, strerror(errno));
  } else {
    snprintf(filename, sizeof(filename), "%s", newfile);
  }
  changed = 0;
  textbuf->call_modify_callbacks();
}

void copy_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Editor::kf_copy(0, e->editor);
}

void cut_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Editor::kf_cut(0, e->editor);
}

void delete_cb(Fl_Widget*, void*) { textbuf->remove_selection(); }

void linenumbers_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Menu_Bar* m = (Fl_Menu_Bar*)w;
  const Fl_Menu_Item* i = m->mvalue();
  if (i->value()) {
    e->editor->linenumber_width(line_num_width);  // enable
    e->editor->linenumber_size(e->editor->textsize());
  } else {
    e->editor->linenumber_width(0);  // disable
  }
  e->line_numbers = (i->value() ? 1 : 0);
  e->redraw();
}

void wordwrap_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Menu_Bar* m = (Fl_Menu_Bar*)w;
  const Fl_Menu_Item* i = m->mvalue();
  if (i->value()) {
    e->editor->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
  } else {
    e->editor->wrap_mode(Fl_Text_Display::WRAP_NONE, 0);
  }
  e->wrap_mode = (i->value() ? 1 : 0);
  e->redraw();
}

void find_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* val;

  val = fl_input("Search String:", e->search);
  if (val != NULL) {
    // User entered a string - go find it!
    snprintf(e->search, sizeof(e->search), "%s", val);
    find2_cb(w, v);
  }
}

void find2_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  if (e->search[0] == '\0') {
    // Search string is blank; get a new one...
    find_cb(w, v);
    return;
  }

  int pos = e->editor->insert_position();
  int found = textbuf->search_forward(pos, e->search, &pos);
  if (found) {
    // Found a match; select and update the position...
    textbuf->select(pos, pos + strlen(e->search));
    e->editor->insert_position(pos + strlen(e->search));
    e->editor->show_insert_position();
  } else {
    fl_alert("No occurrences of \'%s\' found!", e->search);
  }
}

void set_title(Fl_Window* w) {
  if (filename[0] == '\0') {
    snprintf(title, sizeof(title), "Untitled");
  } else {
    char* slash;
    slash = strrchr(filename, '/');
#ifdef WIN32
    if (slash == NULL) {
      slash = strrchr(filename, '\\');
    }
#endif
    if (slash != NULL) {
      snprintf(title, sizeof(title), "%s", slash + 1);
    } else {
      snprintf(title, sizeof(title), "%s", filename);
    }
  }

  if (changed) {
    size_t len = strlen(title);
    snprintf(title + len, sizeof(title) - len, " (modified)");
  }

  w->label(title);
}

void changed_cb(int, int nInserted, int nDeleted, int, const char*, void* v) {
  if ((nInserted || nDeleted) && !loading) {
    changed = 1;
  }
  EditorWindow* w = (EditorWindow*)v;
  set_title(w);
  update_status(w);
  if (loading) {
    w->editor->show_insert_position();
  }
}

void new_cb(Fl_Widget*, void*) {
  if (!check_save()) {
    return;
  }

  filename[0] = '\0';
  textbuf->select(0, textbuf->length());
  textbuf->remove_selection();
  changed = 0;
  textbuf->call_modify_callbacks();
}

void open_cb(Fl_Widget*, void*) {
  if (!check_save()) {
    return;
  }
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Open file");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_FILE);
  if (fnfc.show()) {
    return;
  }
  load_file(fnfc.filename(), -1);
}

void insert_cb(Fl_Widget*, void* v) {
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Insert file");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_FILE);
  if (fnfc.show()) {
    return;
  }
  EditorWindow* w = (EditorWindow*)v;
  load_file(fnfc.filename(), w->editor->insert_position());
}

void paste_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Editor::kf_paste(0, e->editor);
}

void close_cb(Fl_Widget*, void* v) {
  EditorWindow* w = (EditorWindow*)v;

  if (num_windows == 1) {
    if (!check_save()) {
      return;
    }
  }

  w->hide();
  w->editor->buffer(0);
  textbuf->remove_modify_callback(style_update, w->editor);
  textbuf->remove_modify_callback(changed_cb, w);
  Fl::delete_widget(w);

  num_windows--;
  if (!num_windows) {
    exit(0);
  }
}

void quit_cb(Fl_Widget*, void*) {
  if (changed && !check_save()) {
    return;
  }

  exit(0);
}

void replace_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->replace_dlg->show();
}

void keyword_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->keyword_dlg->show();
}

void replace2_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* find = e->replace_find->value();
  const char* replace = e->replace_with->value();

  if (find[0] == '\0') {
    // Search string is blank; get a new one...
    e->replace_dlg->show();
    return;
  }

  e->replace_dlg->hide();

  int pos = e->editor->insert_position();
  int found = textbuf->search_forward(pos, find, &pos);

  if (found) {
    // Found a match; update the position and replace text...
    textbuf->select(pos, pos + strlen(find));
    textbuf->remove_selection();
    textbuf->insert(pos, replace);
    textbuf->select(pos, pos + strlen(replace));
    e->editor->insert_position(pos + strlen(replace));
    e->editor->show_insert_position();
  } else {
    fl_alert("No occurrences of \'%s\' found!", find);
  }
}

void replall_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* find = e->replace_find->value();
  const char* replace = e->replace_with->value();

  find = e->replace_find->value();
  if (find[0] == '\0') {
    // Search string is blank; get a new one...
    e->replace_dlg->show();
    return;
  }

  e->replace_dlg->hide();

  e->editor->insert_position(0);
  int times = 0;

  // Loop through the whole string
  for (int found = 1; found;) {
    int pos = e->editor->insert_position();
    found = textbuf->search_forward(pos, find, &pos);

    if (found) {
      // Found a match; update the position and replace text...
      textbuf->select(pos, pos + strlen(find));
      textbuf->remove_selection();
      textbuf->insert(pos, replace);
      e->editor->insert_position(pos + strlen(replace));
      e->editor->show_insert_position();
      times++;
    }
  }

  if (times) {
    fl_message("Replaced %d occurrences.", times);
  } else {
    fl_alert("No occurrences of \'%s\' found!", find);
  }
}

void replcan_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->replace_dlg->hide();
}

void keyword_close_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->keyword_dlg->hide();
}

void keyword_insert_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  int line = e->keyword_browser->value();
  if (line > 0) {
    std::string toInsert = comp[e->keyword_browser->text(line)];
    textbuf->insert(e->editor->insert_position(), toInsert.c_str());
  }
}

void save_cb() {
  if (filename[0] == '\0') {
    // No filename - get one!
    saveas_cb();
    return;
  } else {
    save_file(filename);
  }
}

void saveas_cb() {
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Save File As?");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
  if (fnfc.show()) {
    return;
  }
  save_file(fnfc.filename());
}

void view_cb(Fl_Widget*, void*) {
  Fl_Window* w = new_view();
  w->show();
}

void goto_start_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->editor->insert_position(0);
  e->editor->show_insert_position();
}

void goto_end_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->editor->insert_position(textbuf->length());
  e->editor->show_insert_position();
}

void goto_line_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* val = fl_input("Go to line:", "1");
  if (val == NULL) {
    return;
  }
  int line = atoi(val);
  if (line < 1) {
    line = 1;
  }
  // skip_lines(0, line - 1) returns the buffer position at the start of 'line'
  int pos = textbuf->skip_lines(0, line - 1);
  e->editor->insert_position(pos);
  e->editor->show_insert_position();
}

// ----- Status bar -----------------------------------------------------------

void update_status(EditorWindow* w) {
  if (!w || !w->status || !w->editor) {
    return;
  }
  int pos = w->editor->insert_position();
  int line = textbuf->count_lines(0, pos) + 1;
  int col = pos - textbuf->line_start(pos) + 1;
  char buf[128];
  snprintf(buf, sizeof(buf), "  Ln %d, Col %d  |  %d chars", line, col, textbuf->length());
  w->status->copy_label(buf);
}

int ConfTextEditor::handle(int e) {
  int r = Fl_Text_Editor::handle(e);
  update_status(win);
  return r;
}

// ----- Theme ----------------------------------------------------------------

int dark_mode = 1;  // dark is the default theme

void apply_theme_to_editor(Fl_Text_Editor* ed) {
  if (dark_mode) {
    ed->color(fl_rgb_color(30, 30, 30));
    ed->textcolor(fl_rgb_color(220, 220, 220));
    ed->selection_color(fl_rgb_color(60, 90, 140));
    ed->cursor_color(fl_rgb_color(80, 170, 255));
    ed->linenumber_bgcolor(fl_rgb_color(58, 58, 58));
    ed->linenumber_fgcolor(fl_rgb_color(140, 140, 140));
  } else {
    ed->color(FL_WHITE);
    ed->textcolor(fl_rgb_color(30, 30, 30));
    ed->selection_color(fl_rgb_color(180, 210, 255));
    ed->cursor_color(fl_rgb_color(0, 120, 215));
    ed->linenumber_bgcolor(fl_rgb_color(240, 240, 240));
    ed->linenumber_fgcolor(fl_rgb_color(150, 150, 150));
  }
}

// Theme the window "chrome" (menu bar + status bar) plus the editor itself.
void apply_theme_to_window(EditorWindow* e) {
  if (!e) {
    return;
  }
  Fl_Color bg, fg;
  if (dark_mode) {
    bg = fl_rgb_color(45, 45, 45);
    fg = fl_rgb_color(220, 220, 220);
  } else {
    bg = fl_rgb_color(238, 238, 238);
    fg = fl_rgb_color(40, 40, 40);
  }
  if (e->editor) {
    apply_theme_to_editor(e->editor);
  }
  if (e->menubar) {
    e->menubar->color(bg);
    e->menubar->textcolor(fg);
  }
  if (e->status) {
    e->status->color(bg);
    e->status->labelcolor(fg);
  }
  e->redraw();
}

void apply_theme() {
  if (dark_mode) {
    // Syntax colors tuned for a dark background
    styletable[0].color = fl_rgb_color(220, 220, 220);  // A - Plain
    styletable[1].color = fl_rgb_color(120, 190, 120);  // B - Comments
    styletable[2].color = fl_rgb_color(90, 200, 200);   // C - Types
    styletable[3].color = fl_rgb_color(110, 160, 255);  // D - Keywords
    Fl::background(45, 45, 45);
    Fl::background2(30, 30, 30);
    Fl::foreground(220, 220, 220);
  } else {
    styletable[0].color = FL_BLACK;       // A - Plain
    styletable[1].color = FL_DARK_GREEN;  // B - Comments
    styletable[2].color = FL_DARK_CYAN;   // C - Types
    styletable[3].color = FL_BLUE;        // D - Keywords
    Fl::background(238, 238, 238);
    Fl::background2(255, 255, 255);
    Fl::foreground(30, 30, 30);
  }
  Fl::reload_scheme();

  // Push the new look to every open editor view
  for (Fl_Window* win = Fl::first_window(); win; win = Fl::next_window(win)) {
    apply_theme_to_window(dynamic_cast<EditorWindow*>(win));
  }
}

void theme_light_cb(Fl_Widget*, void*) {
  dark_mode = 0;
  apply_theme();
}

void theme_dark_cb(Fl_Widget*, void*) {
  dark_mode = 1;
  apply_theme();
}

void ConfTextEditor::draw() {
  Fl_Text_Editor::draw();
  int lnw = linenumber_width();
  if (lnw <= 0) {
    return;  // gutter hidden -> no separator
  }
  // Vertical separator at the right edge of the line-number gutter
  fl_color(dark_mode ? fl_rgb_color(90, 90, 90) : fl_rgb_color(190, 190, 190));
  int sepx = x() + lnw - 1;
  fl_push_clip(x(), y(), w(), h());
  fl_line(sepx, y() + 1, sepx, y() + h() - 2);
  fl_pop_clip();
}

Fl_Window* new_view() {
  EditorWindow* w = new EditorWindow(with0, height0, title);

  const int menuH = 30;
  const int statusH = 24;

  w->begin();
  Fl_Menu_Bar* m = new Fl_Menu_Bar(0, 0, with0, menuH);
  m->copy(menuitems, w);
  m->box(FL_FLAT_BOX);  // honor color() instead of the scheme's gradient
  w->menubar = m;

  ConfTextEditor* ed = new ConfTextEditor(0, menuH, with0, height0 - menuH - statusH);
  ed->win = w;
  w->editor = ed;
  w->editor->textfont(TF);
  w->editor->textsize(TS);
  w->editor->cursor_style(Fl_Text_Display::SIMPLE_CURSOR);
  // w->editor->wrap_mode(Fl_Text_Editor::WRAP_AT_BOUNDS, 250);
  w->editor->buffer(textbuf);

  // Line numbers shown by default (matches the FL_MENU_VALUE flag on the menu toggle)
  w->editor->linenumber_width(line_num_width);
  w->editor->linenumber_size(TS);
  w->editor->linenumber_font(TF);
  w->editor->linenumber_align(FL_ALIGN_RIGHT);
  w->line_numbers = 1;

  w->editor->highlight_data(stylebuf, styletable, sizeof(styletable) / sizeof(styletable[0]), 'A', style_unfinished_cb,
                            0);
  w->plugDialogsWithEditor();

  // Status bar pinned to the bottom (editor is the resizable widget above it)
  w->status = new Fl_Box(0, height0 - statusH, with0, statusH, "Ln 1, Col 1");
  w->status->box(FL_FLAT_BOX);
  w->status->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
  w->status->labelsize(12);
  w->status->labelfont(FL_HELVETICA);

  apply_theme_to_window(w);  // colors for editor + menu bar + status (current theme)

  w->end();
  w->resizable(w->editor);
  w->size_range(300, 200);
  w->callback((Fl_Callback*)close_cb, w);

  textbuf->add_modify_callback(style_update, w->editor);
  textbuf->add_modify_callback(changed_cb, w);
  textbuf->call_modify_callbacks();
  num_windows++;
  return w;
}

void cb(const char* fname) { load_file(fname, -1); }

int main(int argc, char** argv) {
  Fl::scheme("gleam");  // modern, flat theme that honors background colors

  textbuf = new Fl_Text_Buffer;
  // textbuf->transcoding_warning_action = NULL;
  init_keywords();
  init_types();
  init_documentation();
  init_snippets();
  style_init();
  fl_open_callback(cb);
  fl_message_hotspot(1);

  apply_theme();  // set syntax/background colors for the default (dark) theme

  Fl_Window* window = new_view();

  window->show(1, argv);

  // Re-apply after show(): showing the window re-initializes some scheme colors,
  // so we set the theme colors again here to make sure the dark theme wins.
  apply_theme();

  // #ifndef __APPLE__
  if (argc > 1) {
    load_file(argv[1], -1);
  }
  // #endif

  return Fl::run();
}
