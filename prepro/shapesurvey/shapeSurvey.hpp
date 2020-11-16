#ifndef SEE_HPP_E29BD15E
#define SEE_HPP_E29BD15E

#include <GL/freeglut.h>

#include "../Rockable/Shape.hpp"
#include <fstream>
#include <functional>
#include <map>

std::vector<Shape> Shapes;
size_t ishape = 0;  // id of the current shape beeing shown
std::string shapeFileName = "shapes";

int main_window;

// flags
int show_background = 0;
int show_help = 0;

int width = 800;
int height = 800;
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
void drawObbLevel(size_t ishp, size_t level);
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
void exportImageTerrain(const char* fileName);
void exportPositionning();

#endif /* end of include guard: SEE_HPP_E29BD15E */
