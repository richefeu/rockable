#ifndef SEE_HPP
#define SEE_HPP

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#include <GL/glut.h>

#include "Core/Shape.hpp"
#include "Core/ShapeMesh.hpp"
#include "stackTracer.hpp"

#include <fstream>
#include <functional>
#include <map>

std::vector<Shape> Shapes;
size_t ishape = 0;  // id of the current shape beeing shown
std::string shapeFileName = "shapes";

// ---- CLI screenshot (headless one-shot rendering) ----
// When outputFile is non-empty, shapeSurvey renders a single clean frame of the
// selected shape (no axes, no OBB, no text) and writes it to that PNG, then
// exits, instead of opening the interactive window.
std::string outputFile = "";              // empty => interactive mode
std::string skinFile = "";                // optional .rmsh: draw the skin mesh instead of primitives
std::map<std::string, ShapeMesh> skinMeshes;  // keyed by shape name
double shot_elev = 22.0;                   // camera elevation (deg), z is up
double shot_azim = -55.0;                  // camera azimuth (deg)
int shot_transparent = 0;                  // 1 => transparent background (alpha 0)
int show_wire = 0;                         // 1 => overlay the mesh as black wireframe
float shapeColor[3] = {0.761f, 0.733f, 0.976f};  // fill colour (RGB in [0,1])
double radiusOverride = -1.0;              // >=0 => force the Minkowski radius R (0 => raw triangles)

int main_window;

// flags
int show_background = 0;
int show_help = 0;

int width = 700;
int height = 700;
float wh_ratio = (float)width / (float)height;

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int mouse_start[2];
float view_angle;
float znear;
float zfar;
GLfloat Rot_Matrix[16];
GLfloat max_length;
GLfloat alpha = 1.0;

int maxOBBLevel = 0;

vec3r eye;
vec3r center;
vec3r up;

void drawShape(size_t ishp);
void drawSkin(size_t ishp);
void set_view(double elev_deg, double azim_deg);
void screenshotDisplay();
void saveScreenshot(const char* filename);
void recursiveDrawOBB(OBBnode<subBox> * node, int wantedLevel, int level = 0);
void drawObbLevel(size_t ishp, size_t wantedLevel);
void drawFrame();
void drawInfo();

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void buildMenu();
void printHelp();
vec3r rotatePoint(vec3r const& p, vec3r const& center, vec3r const& axis, double theta);
void adjust_clipping_plans();
void fit_view();
int readShapeLib(const char* fileName);
void saveShapeLib(const char* fileName);
void exportSample();

#endif /* end of include guard: SEE_HPP */
