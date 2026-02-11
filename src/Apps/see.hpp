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

#ifndef SEE_HPP
#define SEE_HPP

#include "toofus-gate/tclap/CmdLine.h"

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#include <GL/freeglut.h>

#ifdef WITH_LIBPNG
#include <png.h>
#endif

#include <nlohmann/json.hpp>

#include "ColorTable.hpp"
#include "fileTool.hpp"
#include "geoTool.hpp"
#include "glTools.hpp"
#include "message.hpp"

#include "BlockRelease.hpp"
#include "Core/Rockable.hpp"
#include "ProcessingTools/processingTool_probeSolidFraction.hpp"

Rockable box;
int confNum = 0;

std::vector<BlockRelease> releases;

int main_window;

nlohmann::json params = {{"colorMode", 0},
                         {"rescaleColorRange", 1},
                         {"colorRangeMin", 0.0},
                         {"colorRangeMax", 1.0},
                         {"show_background", 0},
                         {"show_particles", 1},
                         {"show_driven", 1},
                         {"show_velocities", 0},
                         {"show_colorBar", 1},
                         {"show_forces", 0},
                         {"show_normal_forces", 0},
                         {"show_obb", 0},
                         {"enlarged_obb", 0},
                         {"show_interFrames", 0},
                         {"show_interTypes", 0},
                         {"show_keybinds", 0},
                         {"show_traj", 0},
                         {"show_probe", 0},
                         {"show_periodicBox", 1},
                         {"ParticleColor", {207, 174, 85}},
                         {"alpha_particles", 1.0f},
                         {"alpha_fixparticles", 0.15f},
                         {"window", {{"width", 768}, {"height", 768}}},
                         {"camera",
                          {{"view_angle", 45.0f},
                           {"znear", 0.01f},
                           {"zfar", 10.0f},
                           {"eye", {0.0, 0.0, 1.0}},
                           {"center", {0.0, 0.0, 0.0}},
                           {"up", {0.0, 1.0, 0.0}}}}};

// brick color rgb is: 178, 34, 34

AABB probe;
size_t probe_MCnsteps = 10000;

ColorTable CT;
std::vector<colorRGBA> pcolors;

int shapeWithoutThickness = 0;
size_t totalNumberOfVertices = 0;  // it says how the sample will be long to display

int selectedParticle = -1;

// double arrowSize = 0.0005;
// double arrowAngle = 0.7;
// double radiusMin;
// double radiusMax;
// double radiusMean;

int width = 700;
int height = 700;
float wh_ratio = (float)width / (float)height;

glTextZone textZone(3, &width, &height);

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int mouse_start[2];
float view_angle;
float znear;
float zfar;
GLfloat Rot_Matrix[16];
GLfloat max_length;

vec3r eye;
vec3r center;
vec3r up;

void drawShape(Shape* s, double homothety = 1.0, const mat9r& T = mat9r::unit());
void drawForces();
void drawC2CNormalForce();
void drawVelocities();
void drawParticles();
void drawTrajectories();
void drawInteractionTypes();
void drawInteractionFrames();
void drawGlobalFrame();
void drawGlobalAABB();
void drawOBBs();
void drawProbe();
void drawColorBar();

#ifdef ROCKABLE_ENABLE_PERIODIC
void drawPeriodicCell();
#endif

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void local_to_json();
void json_to_local();
void load_see_json();
void save_see_json();
void buildMenu();
vec3r rotatePoint(vec3r const& p, vec3r const& center, vec3r const& axis, double theta);
void computePerspective();
void adjustClippingPlans();
void fitView();
bool fileExists(const char* fileName);
void readTraj(const char* name);
bool tryToReadConf(int num);
int screenshot(const char* filename);
void selection(int x, int y);
void editSelection();
void resetColors(int mode, int rescale = 1);


#endif /* end of include guard: SEE_HPP */
