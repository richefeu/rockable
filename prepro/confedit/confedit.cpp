#include "confedit.hpp"

void init_completion() {
  // comp["____"] = "____";

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
  comp["numericalDampingCoeff"] = "_____";
  comp["Velocity barriers"] =
      "VelocityBarrier _value_\nAngularVelocityBarrier _value_\n"
      "VelocityBarrierExponent _value_\nAngularVelocityBarrierExponent _value_\n";

  comp["Tempo"] = "_____";

  comp["force law"] = "forceLaw _Default/Avalanches/StickedLinks_";

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
}

void init_documentation() {
  // docu["___"] = "____";

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
      "Time step";
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
  docu["numericalDampingCoeff"] = "_____";
  docu["VelocityBarrier"] = "_____";
  docu["AngularVelocityBarrier"] = "_____";
  docu["VelocityBarrierExponent"] = "_____";
  docu["AngularVelocityBarrierExponent"] = "_____";
  docu["Tempo"] = "_____";
  docu["forceLaw"] =
      "forceLaw [(string)Name]\n\n"
      "Select a model for the computation of forces.\n"
      "Possible Names are: 'Default', 'Avalanches' or 'StickedLinks'";

  docu["AddOrRemoveInteractions"] = "_____";
  docu["UpdateNL"] = "_____";
  docu["Integrator"] = "_____";
  docu["cellMinSizes"] = "_____";
  docu["boxForLinkCellsOpt"] = "_____";
  docu["ContactPartnership"] = "_____";

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
  docu["DataExtractor"] = "____";
  docu["BodyForce"] = "____";

  docu["stickVerticesInClusters"] =
      "stickVerticesInClusters [(double) Epsilon]\n\n"
      "This command will add glued interfaces between bodies having\n"
      "the same cluster identifier. Only bonds between vertices (spheres)\n"
      "are added when the distance is less than Epsilon.";
  docu["stickClusters"] =
      "stickClusters [(double) Epsilon]\n\n"
      "This command will add glued interfaces between bodies having different cluster identifiers.\n"
      "Bonds are added when the distance is less than Epsilon.";

  docu["copyParamsToInterfaces"] = "copyParamsToInterfaces [(string) isInnerStr]";
  docu["setStiffnessRatioInterfaces"] = "____";
  docu["setVariableStickParams"] =
      "setVariableStickParams [(string) paramName] [(string) isInnerStr] [(double) lambda] [(double) m] [(int) "
      "timeSeeded]";
  docu["randomlyOrientedVelocities"] = "randomlyOrientedVelocities [(double) velocity]";
  docu["randomlyOrientedVelocitiesClusters"] = "randomlyOrientedVelocitiesClusters [(double) velocity] [(int) option]";
  docu["setAllVelocities"] = "____";
  docu["homothetyRange"] = "____";
  docu["particlesClonage"] = "____";

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
  docu["Servo"] = "____";

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
  docu["OBBtreeLevel"] = "DEPRECATED!!";
  docu["position"] = "____";
  docu["orientation"] = "____";
  docu["isSurface"] = "____";
  docu["MCnstep"] = "____";
  docu["nv"] = "____";
  docu["ne"] = "____";
  docu["nf"] = "____";
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
  code_keywords.insert(std::string("paramsInInterfaces"));
  code_keywords.insert(std::string("dynamicUpdateNL"));
  code_keywords.insert(std::string("dispUpdateNL"));
  code_keywords.insert(std::string("angleUpdateNL"));
  code_keywords.insert(std::string("numericalDampingCoeff"));
  code_keywords.insert(std::string("velocityBarrier"));
  code_keywords.insert(std::string("angularVelocityBarrier"));
  code_keywords.insert(std::string("velocityBarrierExponent"));
  code_keywords.insert(std::string("angularVelocityBarrierExponent"));
  code_keywords.insert(std::string("Tempo"));
  code_keywords.insert(std::string("forceLaw"));
  code_keywords.insert(std::string("AddOrRemoveInteractions"));
  code_keywords.insert(std::string("UpdateNL"));
  code_keywords.insert(std::string("Integrator"));
  code_keywords.insert(std::string("cellMinSizes"));
  code_keywords.insert(std::string("boxForLinkCellsOpt"));
  code_keywords.insert(std::string("ContactPartnership"));

  code_keywords.insert(std::string("knContact"));
  code_keywords.insert(std::string("en2Contact"));
  code_keywords.insert(std::string("en2ContactFromViscRate"));
  code_keywords.insert(std::string("ktContact"));
  code_keywords.insert(std::string("muContact"));
  code_keywords.insert(std::string("krContact"));
  code_keywords.insert(std::string("murContact"));

  code_keywords.insert(std::string("knInnerBond"));
  code_keywords.insert(std::string("en2InnerBond"));
  code_keywords.insert(std::string("ktInnerBond"));
  code_keywords.insert(std::string("fn0InnerBond"));
  code_keywords.insert(std::string("ft0InnerBond"));
  code_keywords.insert(std::string("powInnerBond"));

  code_keywords.insert(std::string("knOuterBond"));
  code_keywords.insert(std::string("ktOuterBond"));
  code_keywords.insert(std::string("krOuterBond"));
  code_keywords.insert(std::string("fn0OuterBond"));
  code_keywords.insert(std::string("ft0OuterBond"));
  code_keywords.insert(std::string("mom0OuterBond"));
  code_keywords.insert(std::string("powOuterBond"));
  code_keywords.insert(std::string("en2OuterBond"));

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
  code_keywords.insert(std::string("stickClusters"));
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

  code_types.insert(std::string("y"));
  code_types.insert(std::string("n"));

  // tempo
  code_types.insert(std::string("NDCoeff"));
  code_types.insert(std::string("Inter"));
  code_types.insert(std::string("Body"));

  // force laws
  code_types.insert(std::string("Default"));
  code_types.insert(std::string("Avalanches"));
  code_types.insert(std::string("StickedLinks"));

  // add or remove interaction
  code_types.insert(std::string("bruteForce"));
  code_types.insert(std::string("OBBtree"));

  // updateNL
  code_types.insert(std::string("linkCells"));

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

void style_parse(const char* text, char* style, int length) {
  for (int i = 0; i < length; i++) {
    if (text[i] != '\n') style[i] = 'A';
  }
  std::istringstream is((char*)text);

  std::string token;
  is >> token;

  while (is.good()) {

    if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
      int from = is.tellg();
      from -= token.size();
      getline(is, token);  // ignore the rest of the current line
      int to = is.tellg();
      for (int pos = from; pos < to; pos++) style[pos] = 'B';
      is >> token;  // next token before continuing
      continue;
    } else {
      std::set<std::string>::iterator it = code_keywords.find(token);
      if (it != code_keywords.end()) {
        int to = is.tellg();
        int from = to - token.size();
        for (int pos = from; pos < to; pos++) style[pos] = 'D';
      } else {
        std::set<std::string>::iterator ittype = code_types.find(token);
        if (ittype != code_types.end()) {
          int to = is.tellg();
          int from = to - token.size();
          for (int pos = from; pos < to; pos++) style[pos] = 'C';
        }
      }
    }

    is >> token;
  }  // while
}

void style_init(void) {
  char* style = new char[textbuf->length() + 1];
  char* text = textbuf->text();

  memset(style, 'A', textbuf->length());
  style[textbuf->length()] = '\0';

  if (!stylebuf) stylebuf = new Fl_Text_Buffer(textbuf->length());

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
  char last,                                    // Last style on line
      *style,                                   // Style data
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

  // Re-parse the changed region; we do this by parsing from the
  // beginning of the previous line of the changed region to the end of
  // the line of the changed region...  Then we check the last
  // style character and keep updating if we have a multi-line
  // comment character...
  start = textbuf->line_start(pos);
  //  if (start > 0) start = textbuf->line_start(start - 1);
  end = textbuf->line_end(pos + nInserted);
  text = textbuf->text_range(start, end);
  style = stylebuf->text_range(start, end);
  if (start == end)
    last = 0;
  else
    last = style[end - start - 1];

  style_parse(text, style, end - start);

  stylebuf->replace(start, end, style);
  ((Fl_Text_Editor*)cbArg)->redisplay_range(start, end);

  if (start == end || last != style[end - start - 1]) {
    // Either the user deleted some text, or the last character
    // on the line changed styles, so reparse the
    // remainder of the buffer...
    free(text);
    free(style);

    end = textbuf->length();
    text = textbuf->text_range(start, end);
    style = stylebuf->text_range(start, end);

    style_parse(text, style, end - start);

    stylebuf->replace(start, end, style);
    ((Fl_Text_Editor*)cbArg)->redisplay_range(start, end);
  }

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
  if (!insert) strcpy(filename, "");
  int r;
  if (!insert)
    r = textbuf->loadfile(newfile);
  else
    r = textbuf->insertfile(newfile, ipos);
  changed = changed || textbuf->input_file_was_transcoded;
  if (r)
    fl_alert("Error reading from file \'%s\':\n%s.", newfile, strerror(errno));
  else if (!insert)
    strcpy(filename, newfile);
  loading = 0;
  textbuf->call_modify_callbacks();
}

void save_file(const char* newfile) {
  if (textbuf->savefile(newfile))
    fl_alert("Error writing to file \'%s\':\n%s.", newfile, strerror(errno));
  else
    strcpy(filename, newfile);
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
  if (i->value())
    e->editor->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
  else
    e->editor->wrap_mode(Fl_Text_Display::WRAP_NONE, 0);
  e->wrap_mode = (i->value() ? 1 : 0);
  e->redraw();
}

void find_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* val;

  val = fl_input("Search String:", e->search);
  if (val != NULL) {
    // User entered a string - go find it!
    strcpy(e->search, val);
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
  } else
    fl_alert("No occurrences of \'%s\' found!", e->search);
}

void set_title(Fl_Window* w) {
  if (filename[0] == '\0')
    strcpy(title, "Untitled");
  else {
    char* slash;
    slash = strrchr(filename, '/');
#ifdef WIN32
    if (slash == NULL) slash = strrchr(filename, '\\');
#endif
    if (slash != NULL)
      strcpy(title, slash + 1);
    else
      strcpy(title, filename);
  }

  if (changed) strcat(title, " (modified)");

  w->label(title);
}

void changed_cb(int, int nInserted, int nDeleted, int, const char*, void* v) {
  if ((nInserted || nDeleted) && !loading) changed = 1;
  EditorWindow* w = (EditorWindow*)v;
  set_title(w);
  if (loading) w->editor->show_insert_position();
}

void new_cb(Fl_Widget*, void*) {
  if (!check_save()) return;

  filename[0] = '\0';
  textbuf->select(0, textbuf->length());
  textbuf->remove_selection();
  changed = 0;
  textbuf->call_modify_callbacks();
}

void open_cb(Fl_Widget*, void*) {
  if (!check_save()) return;
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Open file");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_FILE);
  if (fnfc.show()) return;
  load_file(fnfc.filename(), -1);
}

void insert_cb(Fl_Widget*, void* v) {
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Insert file");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_FILE);
  if (fnfc.show()) return;
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
    if (!check_save()) return;
  }

  w->hide();
  w->editor->buffer(0);
  textbuf->remove_modify_callback(style_update, w->editor);
  textbuf->remove_modify_callback(changed_cb, w);
  Fl::delete_widget(w);

  num_windows--;
  if (!num_windows) exit(0);
}

void quit_cb(Fl_Widget*, void*) {
  if (changed && !check_save()) return;

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
  } else
    fl_alert("No occurrences of \'%s\' found!", find);
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

  if (times)
    fl_message("Replaced %d occurrences.", times);
  else
    fl_alert("No occurrences of \'%s\' found!", find);
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
  } else
    save_file(filename);
}

void saveas_cb() {
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Save File As?");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
  if (fnfc.show()) return;
  save_file(fnfc.filename());
}

void view_cb(Fl_Widget*, void*) {
  Fl_Window* w = new_view();
  w->show();
}

Fl_Window* new_view() {
  EditorWindow* w = new EditorWindow(with0, height0, title);

  w->begin();
  Fl_Menu_Bar* m = new Fl_Menu_Bar(0, 0, with0, 30);
  m->copy(menuitems, w);
  w->editor = new Fl_Text_Editor(0, 30, with0, height0 - 30);
  w->editor->textfont(FL_SCREEN);
  w->editor->textsize(TS);
  w->editor->cursor_style(Fl_Text_Display::HEAVY_CURSOR);
  w->editor->cursor_color(FL_RED);
  // w->editor->wrap_mode(Fl_Text_Editor::WRAP_AT_BOUNDS, 250);
  w->editor->buffer(textbuf);
  w->editor->highlight_data(stylebuf, styletable, sizeof(styletable) / sizeof(styletable[0]), 'A', style_unfinished_cb,
                            0);
  w->plugDialogsWithEditor();

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
  textbuf = new Fl_Text_Buffer;
  // textbuf->transcoding_warning_action = NULL;
  init_keywords();
  init_types();
  init_documentation();
  init_completion();
  style_init();
  fl_open_callback(cb);
  fl_message_hotspot(1);

  Fl_Window* window = new_view();

  window->show(1, argv);

  //#ifndef __APPLE__
  if (argc > 1) {
    load_file(argv[1], -1);
  }
  //#endif

  return Fl::run();
}
