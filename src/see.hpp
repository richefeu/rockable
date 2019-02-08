//        Rockable, 3D-DEM with sphero-polyhedra
//        Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>
//        
//        This program is free software: you can redistribute it and/or modify
//        it under the terms of the GNU General Public License as published by
//        the Free Software Foundation, either version 3 of the License, or
//        (at your option) any later version.
//        
//        This program is distributed in the hope that it will be useful,
//        but WITHOUT ANY WARRANTY; without even the implied warranty of
//        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//        GNU General Public License for more details.
//        
//        You should have received a copy of the GNU General Public License
//        along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef SEE_HPP_E29BD15E
#define SEE_HPP_E29BD15E

#include <tclap/CmdLine.h>

# ifdef __APPLE__
#  include <GLUT/glut.h>
# else
#  include <GL/glut.h>
# endif

# ifdef PNG_SCREENSHOTS
#  include <png.h>
# endif

#include "message.hpp"
#include "glutTools.hpp"
#include "geoTool.hpp"
#include "fileTool.hpp"

#include "Rockable.hpp"
#include "BlockRelease.hpp"

Rockable box;
int confNum = 0;

std::vector<BlockRelease> releases;

int main_window;

// flags with default values
int show_background  = 1;
int show_particles   = 1;
int show_velocities  = 0;
int show_slice       = 0;
int show_forces      = 0;
int show_obb         = 0;
int enlarged_obb     = 0;
int link_obb         = 1;
int show_interFrames = 0;
int show_interTypes  = 0;
int show_keybinds    = 0;
int show_traj        = 0;

int complexMode = 0;
size_t complexityNumber = 0; // it says how the sample will be long to display

GLfloat alpha_particles = 1.0f;
GLfloat alpha_fixparticles = 0.1f;

int selectedParticle = -1;

double arrowSize   = 0.0005;
double arrowAngle  = 0.7;

double forceTubeFactor = 1.0;

double radiusMin;
double radiusMax;
double radiusMean;

int width  = 800;
int height = 800;
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

void drawShape(Shape * s, double homothety = 1.0);
void drawForces();
void drawVelocities();
void drawParticles();
void drawTrajectories();
void drawInteractionTypes();
void drawInteractionFrames();
void drawOBBs();

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion (int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void buildMenu();
vec3r rotatePoint(vec3r const &p, vec3r const &center, vec3r const &axis, double theta);
void computePerspective();
void adjustClippingPlans();
void fitView();
bool fileExists(const char * fileName);
void readTraj(const char * name);
bool tryToReadConf(int num);
int screenshot(const char *filename);
void selection(int x, int y);
void editSelection();

#endif /* end of include guard: SEE_HPP_E29BD15E */
