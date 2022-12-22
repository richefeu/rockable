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

#ifndef SEE3_HPP
#define SEE3_HPP
#include <tclap/CmdLine.h>
#include <unistd.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#ifdef PNG_H
#include <png.h>
#endif

#include "ColorTable.hpp"
#include "fileTool.hpp"
#include "geoTool.hpp"
#include "glTools.hpp"
#include "message.hpp"

#include "BlockRelease.hpp"
#include "Rockable.hpp"

Rockable box;
int confNum = 0;
int lastConfNumOK = confNum;

std::vector<BlockRelease> releases;

int main_window;

char cwd[256];
// char textBuf[256];
bool ImGui_window_focused = true;
bool show_demo_window = false;
bool show_another_window = false;
bool request_quit = false;

float clear_bottom_color[3] = {135.0f / 255.0f, 206.0f / 255.0f, 250.0f / 255.0f};
float clear_top_color[3] = {1.0f, 1.0f, 1.0f};

// flags with default values
int show_background = 1;
int show_globalFrame = 1;
int show_particles = 1;
int show_driven = 1;
int show_velocities = 0;
int show_slice = 0;
int show_forces = 0;
int show_obb = 0;
int enlarged_obb = 0;
int link_obb = 1;
int show_interFrames = 0;
int show_interTypes = 0;
int show_keybinds = 0;
int show_traj = 0;
int show_probe = 0;

bool show_window_conf_data = false;
bool show_window_info = false;

double forceFactor = 1.0;

AABB probe;
size_t probe_MCnsteps = 10000;

ColorTable CT;
std::vector<colorRGBA> pcolors;
int colorMode = 0;
int rescaleColorRange = 1;
double colorRangeMin = 0.0;
double colorRangeMax = 1.0;

int complexMode = 0;
size_t complexityNumber = 0;  // it says how the sample will be long to display

GLfloat alpha_particles = 1.0f;
GLfloat alpha_fixparticles = 0.1f;

int selectedParticle = -1;

double arrowSize = 0.0005;
double arrowAngle = 0.7;

double forceTubeFactor = 1.0;

double radiusMin;
double radiusMax;
double radiusMean;

int width = 1100;
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

void drawShape(Shape* s, double homothety = 1.0);
void drawForces();
void drawVelocities();
void drawParticles();
void drawTrajectories();
void drawInteractionTypes();
void drawInteractionFrames();
void drawGlobalFrame();
void drawOBBs();
void drawProbe();

// Callback functions
void display();
void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse(GLFWwindow* window, int button, int action, int mods);
void motion(GLFWwindow* window, double x, double y);
void reshape(GLFWwindow* window, int x, int y);

// Helper functions
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

#endif /* end of include guard: SEE3_HPP */
