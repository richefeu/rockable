#ifndef CONFEDIT_HPP
#define CONFEDIT_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctype.h>
#include <errno.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#ifdef __MWERKS__
#define FL_DLL
#endif

#include <FL/Fl.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Native_File_Chooser.H>
#include <FL/Fl_Return_Button.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Text_Editor.H>
#include <FL/filename.H>
#include <FL/fl_ask.H>
#include <FL/x.H>  // for fl_open_callback

#include "quat.hpp"
#include "transformation.hpp"
Transformation<double> globalTransformation;
quat individualParticleRotation;
#include "addParticle.hpp"

const int with0 = 1000;
const int height0 = 800;

int changed = 0;
char filename[FL_PATH_MAX] = "";
char title[FL_PATH_MAX];
Fl_Text_Buffer* textbuf = 0;

// width of line number display, if enabled
const int line_num_width = 60;

// Syntax highlighting stuff...
#define TS 14  // default editor textsize
Fl_Text_Buffer* stylebuf = 0;
Fl_Text_Display::Style_Table_Entry styletable[] = {
    // Style table
    {FL_BLACK, FL_SCREEN, TS},        // A - Plain
    {FL_DARK_YELLOW, FL_SCREEN, TS},  // B - Line comments
    {FL_DARK_CYAN, FL_SCREEN, TS},    // C - Types
    {FL_BLUE, FL_SCREEN, TS},         // D - Keywords
};

std::set<std::string> code_keywords;
std::set<std::string> code_types;
std::map<std::string, std::string> docu;  // documentation
std::map<std::string, std::string> comp;  // completion


//void cb(const char* fname);
void save_cb();
void saveas_cb();

void new_cb(Fl_Widget*, void*);
void open_cb(Fl_Widget*, void*);
void insert_cb(Fl_Widget*, void* v);
void paste_cb(Fl_Widget*, void* v);
void close_cb(Fl_Widget*, void* v);
void quit_cb(Fl_Widget*, void*);
void cut_cb(Fl_Widget*, void* v);
void copy_cb(Fl_Widget*, void* v);
void linenumbers_cb(Fl_Widget* w, void* v);
void wordwrap_cb(Fl_Widget* w, void* v);

void find_cb(Fl_Widget* w, void* v);
void find2_cb(Fl_Widget*, void*);
void replall_cb(Fl_Widget*, void*);
void replace_cb(Fl_Widget*, void*);
void replace2_cb(Fl_Widget*, void*);
void replcan_cb(Fl_Widget*, void*);
void view_cb(Fl_Widget*, void*);

void doc_selection_cb(Fl_Widget*, void* v);
void deg_to_rad_cb(Fl_Widget*, void* v);
void quat_axis_angle_cb(Fl_Widget*, void* v);
void add_particle_cb(Fl_Widget*, void* v);

void keyword_cb(Fl_Widget*, void* v);
void keyword_close_cb(Fl_Widget*, void* v);
void keyword_insert_cb(Fl_Widget*, void* v);

class addParticleDialog : public Fl_Window {
 public:
  Fl_Input* position;
  Fl_Input* rot_axis;
  Fl_Input* rot_angle_deg;
  Fl_Button* insert;
  Fl_Button* close_but;

  Fl_Text_Editor* editor;

  addParticleDialog(int w, int h, const char* t) : Fl_Window(w, h, t) {
    position = new Fl_Input(160, 10, 210, 25, "Position:");
    rot_axis = new Fl_Input(160, 40, 210, 25, "Rotation axis:");
    rot_angle_deg = new Fl_Input(160, 70, 210, 25, "Rotation angle (deg):");

    insert = new Fl_Button(w - 70 - 60 - 10, h - 40, 60, 25, "Insert");
    insert->callback((Fl_Callback*)addParticleDialog::insert_cb, this);

    close_but = new Fl_Button(w - 70, h - 40, 60, 25, "Close");
    close_but->callback((Fl_Callback*)addParticleDialog::close_cb, this);
  }
  
  ~addParticleDialog() {
    
  }

  static void close_cb(Fl_Widget*, void* v) {
    addParticleDialog* w = (addParticleDialog*)v;
    w->hide();
  }

  static void insert_cb(Fl_Widget*, void* v) {
    addParticleDialog* w = (addParticleDialog*)v;

    std::stringstream ss(w->rot_axis->value());
    double axis_x, axis_y, axis_z;
    ss >> axis_x >> axis_y >> axis_z;

    std::stringstream ssangle(w->rot_angle_deg->value());
    double angleDeg;
    ssangle >> angleDeg;

    double s, vx, vy, vz;
    double half_angle = 0.5 * angleDeg * M_PI / 180.;
    double len = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
    s = cos(half_angle);
    vx = sin(half_angle) * axis_x / len;
    vy = sin(half_angle) * axis_y / len;
    vz = sin(half_angle) * axis_z / len;

    char txt[512];
    snprintf(txt, 512, "_shapeName_ _group_ _cluster_ _homothety_  %s  0 0 0  0 0 0  %f %f %f %f  0 0 0  0 0 0\n",
            w->position->value(), s, vx, vy, vz);
    w->editor->buffer()->insert(w->editor->insert_position(), txt);
    w->hide();
  }
};
// Fl_Text_Editor* addParticleDialog::editor = 0;

class EditorWindow : public Fl_Double_Window {
 public:
  EditorWindow(int w, int h, const char* t);
  ~EditorWindow();

  Fl_Window* replace_dlg;
  Fl_Input* replace_find;
  Fl_Input* replace_with;
  Fl_Button* replace_all;
  Fl_Return_Button* replace_next;
  Fl_Button* replace_cancel;

  Fl_Window* keyword_dlg;
  Fl_Browser* keyword_browser;
  Fl_Button* keyword_insert;
  Fl_Button* keyword_close;

  addParticleDialog* add_particle_dlg;

  int wrap_mode;
  int line_numbers;

  Fl_Text_Editor* editor;
  char search[256];

  void plugDialogsWithEditor() { add_particle_dlg->editor = editor; }
};

int loading = 0;
int num_windows = 0;
Fl_Window* new_view();

Fl_Menu_Item menuitems[] = {
    {"&File", 0, 0, 0, FL_SUBMENU},
    {"&New File", 0, (Fl_Callback*)new_cb},
    {"&Open File...", FL_COMMAND + 'o', (Fl_Callback*)open_cb},
    {"&Insert File...", FL_COMMAND + 'i', (Fl_Callback*)insert_cb, 0, FL_MENU_DIVIDER},
    {"&Save File", FL_COMMAND + 's', (Fl_Callback*)save_cb},
    {"Save File &As...", FL_COMMAND + FL_SHIFT + 's', (Fl_Callback*)saveas_cb, 0, FL_MENU_DIVIDER},
    {"New &View",
     FL_ALT
#ifdef __APPLE__
         + FL_COMMAND
#endif
         + 'v',
     (Fl_Callback*)view_cb, 0},
    {"&Close View", FL_COMMAND + 'w', (Fl_Callback*)close_cb, 0, FL_MENU_DIVIDER},
    {"E&xit", FL_COMMAND + 'q', (Fl_Callback*)quit_cb, 0},
    {0},

    {"&Edit", 0, 0, 0, FL_SUBMENU},
    {"Cu&t", FL_COMMAND + 'x', (Fl_Callback*)cut_cb},
    {"&Copy", FL_COMMAND + 'c', (Fl_Callback*)copy_cb},
    {"&Paste", FL_COMMAND + 'v', (Fl_Callback*)paste_cb},
    //{"&Delete", 0, (Fl_Callback*)delete_cb},
    {"&Documentation", FL_COMMAND + 'd', (Fl_Callback*)doc_selection_cb},
    {"Preferences", 0, 0, 0, FL_SUBMENU},
    {"Line Numbers", FL_COMMAND + 'l', (Fl_Callback*)linenumbers_cb, 0, FL_MENU_TOGGLE},
    {"Word Wrap", 0, (Fl_Callback*)wordwrap_cb, 0, FL_MENU_TOGGLE},
    {0},
    {0},

    {"&Search", 0, 0, 0, FL_SUBMENU},
    {"&Find...", FL_COMMAND + 'f', (Fl_Callback*)find_cb},
    {"F&ind Again", FL_COMMAND + 'g', (Fl_Callback*)find2_cb},
    {"&Replace...", FL_COMMAND + 'r', (Fl_Callback*)replace_cb},
    {"Re&place Again", FL_COMMAND + 't', (Fl_Callback*)replace2_cb},
    {0},

    {"T&ools", 0, 0, 0, FL_SUBMENU},
    {"Sni&ppets...", FL_COMMAND + 'p', (Fl_Callback*)keyword_cb},
    {"addParticle...", 0, (Fl_Callback*)add_particle_cb},
    {"Degree to Radian...", 0, (Fl_Callback*)deg_to_rad_cb},
    {"(Axis + angle) to quaternion...", 0, (Fl_Callback*)quat_axis_angle_cb},
    {0},

    {0}};

#endif /* end of include guard: CONFEDIT_HPP */
