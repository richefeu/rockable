#include "shapeSurvey.hpp"

#include "fileTool.hpp"
#include "glTools.hpp"

void printHelp() {
  switch2D::go(width, height);

  glColor4f(1.0f, 1.0f, 1.0f, 0.6f);
  glBegin(GL_QUADS);
  int nbLines = 17;  // update this value when a line is added
  int by = height - nbLines * 15 - 3;
  glVertex2i(0, height);
  glVertex2i(width, height);
  glVertex2i(width, by);
  glVertex2i(0, by);
  glEnd();

  glColor3i(0, 0, 0);
  int dhline = -15;
  int hline = height;
#define _nextLine_ (hline += dhline)
  glText::print(15, _nextLine_, "[A][a]   Tune alpha (transparency)");
  glText::print(15, _nextLine_, "[b]      Background (color gradient) on/off");
  glText::print(15, _nextLine_,
                "[c]      Compute mass properties of the current shape (only "
                "if preCompDone = n)");
  glText::print( 15, _nextLine_,
                "[C]      Compute mass properties of ALL shapes (only if "
                "preCompDone = n)");
  glText::print(15, _nextLine_, "[d]      delete duplicated edges in all shapes");
  glText::print(15, _nextLine_, "[e]      print extents of the current shape");
  glText::print(15, _nextLine_, "[h]      Show this help");
  glText::print(15, _nextLine_, "[K][k]   Tune the level of displayed OBB-tree");
  glText::print(15, _nextLine_, "[N][n]   Tune number of Monte-Carlo steps");
  glText::print(15, _nextLine_, "[p]      Export as particles (Rockable sample)");
  glText::print(15, _nextLine_, "[q]      Quit");
  glText::print(15, _nextLine_, "[s]      Save the shape library");
  glText::print(15, _nextLine_, "[t]      Compute the OBB-tree of the current shape");
  glText::print(15, _nextLine_, "[W][w]   Rotate arround the view axis");
  glText::print(15, _nextLine_, "[+][-]   Navigate through the shapes");
  glText::print(15, _nextLine_, "[*]      reset preCompDone to 'n'");
#undef _nextLine_

  switch2D::back();
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

    case 'A': {
      if (alpha < 1.0f) {
        alpha += 0.05f;
        if (alpha > 1.0f) {
          alpha = 1.0f;
        }
      }
    } break;
    case 'a': {
      if (alpha >= 0.1f) alpha -= 0.05f;
    } break;

    case 'b': {
      show_background = 1 - show_background;
    } break;

    case 'c': {
      if (Shapes[ishape].preCompDone == 'n') {
        Shapes[ishape].massProperties();
        Shapes[ishape].preCompDone = 'y';
        fit_view();
      }
    } break;

    case 'C': {
      for (size_t i = 0; i < Shapes.size(); i++) {
        if (Shapes[i].preCompDone == 'n') {
          Shapes[i].massProperties();
          Shapes[i].preCompDone = 'y';
        }
      }
      fit_view();
    } break;

    case 'd': {
      for (size_t i = 0; i < Shapes.size(); i++) {
        Shapes[i].clean();
      }
    } break;

    case 'e': {
      std::cout << "extents: " << Shapes[ishape].obb.extent[0] << " " << Shapes[ishape].obb.extent[1] << " "
                << Shapes[ishape].obb.extent[2] << std::endl;
    } break;

    case 'h': {
      show_help = 1 - show_help;
    } break;

    case 'K': {
      if (maxOBBLevel < /*Shapes[ishape].OBBtreeLevel*/ 10) {
        maxOBBLevel += 1;
      }
    } break;
    case 'k': {
      if (maxOBBLevel > 0) {
        maxOBBLevel -= 1;
      }
    } break;

    case 'N': {
      if (Shapes[ishape].MCnstep < 100000000) {
        Shapes[ishape].MCnstep *= 10;
      }
    } break;
    case 'n': {
      if (Shapes[ishape].MCnstep > 1000) {
        Shapes[ishape].MCnstep = (size_t)floor(Shapes[ishape].MCnstep * 0.1);
      }
    } break;
    
    case 'o': {
      Shapes[ishape].fibObbOption += 1;
      if (Shapes[ishape].fibObbOption > 3) { Shapes[ishape].fibObbOption = 0; }
      if (Shapes[ishape].fibObbOption == 0) {std::cout << "fit OBB strategy = COVARIANCE" << std::endl; }
      if (Shapes[ishape].fibObbOption == 1) {std::cout << "fit OBB strategy = MIN_VOLUME" << std::endl; }
      if (Shapes[ishape].fibObbOption == 2) {std::cout << "fit OBB strategy = IS_AABB" << std::endl; }
      if (Shapes[ishape].fibObbOption == 3) {std::cout << "fit OBB strategy = IMPOSED_AXIS" << std::endl; }
    } break;

    case 'p': {
      exportSample();
    } break;

    case 'q': {
      exit(0);
    } break;

    case 's': {
      saveShapeLib(shapeFileName.c_str());
    } break;

    case 't': {
      maxOBBLevel = 0;
      Shapes[ishape].buildOBBtree();
    } break;

    case 'w': {
      vec3r axis = center - eye;
      axis.normalize();

      up = rotatePoint(up, eye, axis, -0.005 * M_PI);
      up.normalize();
    } break;

    case 'W': {
      vec3r axis = center - eye;
      axis.normalize();

      up = rotatePoint(up, eye, axis, 0.005 * M_PI);
      up.normalize();
    } break;

    case '-': {
      if (ishape > 0) ishape--;
      if (Shapes[ishape].preCompDone == 'n') Shapes[ishape].fitObb();
      fit_view();
    } break;

    case '+': {
      ishape++;
      if (ishape >= Shapes.size()) ishape = Shapes.size() - 1;
      if (Shapes[ishape].preCompDone == 'n') Shapes[ishape].fitObb();
      fit_view();
    } break;

    case '*': {
      if (Shapes[ishape].preCompDone == 'y') Shapes[ishape].preCompDone = 'n';
    } break;
  };

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    //display();
    glutPostRedisplay();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON: {
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
          mouse_mode = PAN;
        } else {
          mouse_mode = ROTATION;
        }
      } break;
      case GLUT_MIDDLE_BUTTON: {
        mouse_mode = ZOOM;
      } break;
    }
  }
}

vec3r rotatePoint(vec3r const& p, vec3r const& center_, vec3r const& axis, double theta) {
  double const c = cos(theta), s = sin(theta);
  double const C = 1.0 - c;
  vec3r tmp = p - center_;
  return center_ + vec3r(tmp[0] * (axis[0] * axis[0] * C + c) + tmp[1] * (axis[0] * axis[1] * C - axis[2] * s) +
                             tmp[2] * (axis[0] * axis[2] * C + axis[1] * s),
                         tmp[0] * (axis[1] * axis[0] * C + axis[2] * s) + tmp[1] * (axis[1] * axis[1] * C + c) +
                             tmp[2] * (axis[1] * axis[2] * C - axis[0] * s),
                         tmp[0] * (axis[2] * axis[0] * C - axis[1] * s) +
                             tmp[1] * (axis[2] * axis[1] * C + axis[0] * s) + tmp[2] * (axis[2] * axis[2] * C + c));
}

void motion(int x, int y) {
  if (mouse_mode == NOTHING) {
    return;
  }

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;
  double length;
  vec3r axis;

  switch (mouse_mode) {

    case ROTATION: {
      axis = (cross(up, center - eye));
      axis.normalize();
      eye = rotatePoint(eye, center, up, -dx * M_PI);
      eye = rotatePoint(eye, center, axis, dy * M_PI);
      up = (rotatePoint((center + up), center, axis, dy * M_PI) - center);
      up.normalize();
    } break;

    case ZOOM: {
      eye = center + (eye - center) * (dy + 1.0);
    } break;

    case PAN: {
      length = (eye - center).length() * tan(view_angle * M_PI / 360.0) * 2.0;
      axis = cross(up, center - eye);
      axis.normalize();
      center = center + axis * dx * length * 0.8;
      center = center + up * dy * length;
    } break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  //display();
  glutPostRedisplay();
}

void drawInfo() {
  switch2D::go(width, height);
  glColor3f(1.0f, 0.388f, 0.278f);  // dark-orange

  glText::print(10, 10, "Shape %lu/%lu, named %s", ishape + 1, Shapes.size(),
                Shapes[ishape].name.c_str());
  glText::print(10, 25, "Radius = %g, OBBtreeLevel = %d", Shapes[ishape].radius, maxOBBLevel);
  glText::print(10, 40, "nb vertex = %lu, nb edge = %lu, nb face = %lu",
                Shapes[ishape].vertex.size(), Shapes[ishape].edge.size(), Shapes[ishape].face.size());
  glText::print(10, 55, "preCompDone %c", Shapes[ishape].preCompDone);

  if (Shapes[ishape].preCompDone == 'y') glColor3f(0.153f, 0.486f, 0.22f);  // green

  glText::print(10, 70, "MCnstep = %lu, Volume = %g", Shapes[ishape].MCnstep,
                Shapes[ishape].volume);
  std::string v = "Solid";
  if (Shapes[ishape].isSurface) v = "Surface";

  glText::print(10, 85, "I/m %g %g %g (%s)", Shapes[ishape].inertia_mass[0],
                Shapes[ishape].inertia_mass[1], Shapes[ishape].inertia_mass[2], v.c_str());

  switch2D::back();
}

void display() {
  glTools::clearBackground(show_background);
  adjust_clipping_plans();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

  drawFrame();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);

  drawShape(ishape);

  glColor3f(0.8f, 0.11f, 0.78f);
  glShape::obb(Shapes[ishape].obb);

  drawObbLevel(ishape, maxOBBLevel);

  drawInfo();
  if (show_help) {
    printHelp();
  }

  glFlush();
  glutSwapBuffers();
}

void adjust_clipping_plans() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  wh_ratio = (float)width / (float)height;
  double zf = (eye - center).normalize();
  OBB& obb = Shapes[ishape].obb;
  vec3r mx = 2 * (obb.extent[0] * obb.e[0] + obb.extent[1] * obb.e[1] + obb.extent[2] * obb.e[2]);
  max_length = (GLfloat)(2 * norm(mx));
  znear = zf - 0.5 * max_length;
  double close_dst = 0.1 * zf;
  if (znear < close_dst) znear = close_dst;
  zfar = zf + max_length;
  gluPerspective(view_angle, wh_ratio, znear, zfar);
  glMatrixMode(GL_MODELVIEW);
}

void fit_view() {
  vec3r dir = (eye - center);
  OBB& obb = Shapes[ishape].obb;
  vec3r diag = 2.0 * (obb.extent[0] * obb.e[0] + obb.extent[1] * obb.e[1] + obb.extent[2] * obb.e[2]);
  dir.normalize();
  center = obb.center;
  GLfloat d = 0.5 * diag.length() / (atan(view_angle * M_PI / 360.0));
  eye = center + d * dir;
}

void reshape(int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, width, height);

  adjust_clipping_plans();
  glutPostRedisplay();
}

void recursiveDrawOBB(OBBnode<subBox>* node, int wantedLevel, int level) {
  if (node == nullptr) {
    return;
  }

  if (level == wantedLevel) {
    glShape::obb(node->boundary);
  }

  if (node->first != nullptr) {
    recursiveDrawOBB(node->first, wantedLevel, level + 1);
  }
  if (node->second != nullptr) {
    recursiveDrawOBB(node->second, wantedLevel, level + 1);
  }
  return;
}

void drawObbLevel(size_t ishp, size_t wantedLevel) {
  glDisable(GL_LIGHTING);
  glColor3f(0.2f, 0.2f, 0.8f);

  recursiveDrawOBB(Shapes[ishp].tree.root, wantedLevel);
}

void drawFrame() {
  OBB& obb = Shapes[ishape].obb;
  vec3r diag = 2 * (obb.extent[0] * obb.e[0] + obb.extent[1] * obb.e[1] + obb.extent[2] * obb.e[2]);
  double len = diag.length() * 0.333;

  glColor3f(1.0f, 0.0f, 0.0f);
  glShape::arrow(vec3r::zero(), len * vec3r::unit_x());
  glColor3f(0.0f, 1.0f, 0.0f);
  glShape::arrow(vec3r::zero(), len * vec3r::unit_y());
  glColor3f(0.0f, 0.0f, 1.0f);
  glShape::arrow(vec3r::zero(), len * vec3r::unit_z());
}

void drawShape(size_t ishp) {
  if (ishp >= Shapes.size()) {
    return;
  }
  if (mouse_mode != NOTHING) {
    return;
  }

  double R = Shapes[ishp].radius;
  // glColor4f(0.666f, 0.729f, 0.09f, alpha); // yellow
  glColor4f(.761f, .733f, .976f, alpha);

  size_t nv = Shapes[ishp].vertex.size();
  for (size_t v = 0; v < nv; ++v) {
    vec3r pos = Shapes[ishp].vertex[v];
    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);
    glShape::sphere(R, 3);
    glPopMatrix();
  }

  size_t ne = Shapes[ishp].edge.size();
  for (size_t e = 0; e < ne; ++e) {
    size_t ideb = Shapes[ishp].edge[e].first;
    size_t iend = Shapes[ishp].edge[e].second;
    vec3r orig = Shapes[ishp].vertex[ideb];
    vec3r arrow = Shapes[ishp].vertex[iend] - orig;
    glShape::tube(orig, arrow, 2.0 * R);
  }

  size_t nf = Shapes[ishp].face.size();
  for (size_t f = 0; f < nf; ++f) {
    if (Shapes[ishp].face[f].size() < 3) {
      continue;
    }  // At least 3 pts!
    vec3r N = cross(Shapes[ishp].vertex[Shapes[ishp].face[f][1]] - Shapes[ishp].vertex[Shapes[ishp].face[f][0]],
                    Shapes[ishp].vertex[Shapes[ishp].face[f][2]] - Shapes[ishp].vertex[Shapes[ishp].face[f][0]]);
    N.normalize();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(N.x, N.y, N.z);
    for (size_t v = 0; v < Shapes[ishp].face[f].size(); ++v) {
      glVertex3f(Shapes[ishp].vertex[Shapes[ishp].face[f][v]].x + N.x * R,
                 Shapes[ishp].vertex[Shapes[ishp].face[f][v]].y + N.y * R,
                 Shapes[ishp].vertex[Shapes[ishp].face[f][v]].z + N.z * R);
    }
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(-N.x, -N.y, -N.z);
    for (size_t v = 0; v < Shapes[ishp].face[f].size(); ++v) {
      glVertex3f(Shapes[ishp].vertex[Shapes[ishp].face[f][v]].x - N.x * R,
                 Shapes[ishp].vertex[Shapes[ishp].face[f][v]].y - N.y * R,
                 Shapes[ishp].vertex[Shapes[ishp].face[f][v]].z - N.z * R);
    }
    glEnd();
  }
}

int readShapeLib(const char* fileName) {
  if (!fileTool::fileExists(fileName)) {
    std::cout << "Shape Library named '" << fileName << "' has not been found." << std::endl;
    return 0;
  }
  shapeFileName = std::string(fileName);

  std::ifstream is(fileName);

  std::string token;
  is >> token;
  while (is) {
    if (token == "<") {
      Shape S;
      S.read(is);
      Shapes.push_back(S);
    }
    is >> token;
  }

  std::cout << "Number of Shapes found: " << Shapes.size() << std::endl;

  ishape = 0;
  if (Shapes[ishape].preCompDone == 'n') {
    Shapes[ishape].fitObb();
  }
  OBB& obb = Shapes[ishape].obb;
  center.set(obb.center.x, obb.center.y, obb.center.y);  // where we look at
  eye.set(obb.center.x + obb.extent.x, obb.center.y,
          obb.center.y);  // from where we look
  up.set(0.0, 1.0, 0.0);  // direction (normalized)

  return 1;
}

void saveShapeLib(const char* fileName) {
  std::ofstream os;

  if (fileTool::fileExists(fileName)) {
    std::string newFileName = fileTool::GetFilePath(fileName) + "/mod_" + fileTool::GetFileName(fileName) + "." +
                              fileTool::GetFileExt(fileName);
    std::cout << "save " << newFileName << std::endl;
    os.open(newFileName);
  } else {
    std::cout << "save " << fileName << std::endl;
    os.open(fileName);
  }

  for (size_t s = 0; s < Shapes.size(); s++) {
    Shapes[s].write(os);
  }
}

void exportSample() {
  std::cout << "Find sourrounding box " << std::endl;
  vec3r minBox(1e12, 1e12, 1e12);
  vec3r maxBox(-1e12, -1e12, -1e12);

  for (size_t i = 0; i < Shapes.size(); i++) {
    for (size_t isub = 0; isub < Shapes[i].vertex.size(); isub++) {
      vec3r pos = Shapes[i].position + Shapes[i].orientation * Shapes[i].vertex[isub];
      vec3r rad(Shapes[i].radius, Shapes[i].radius, Shapes[i].radius);
      vec3r posMin = pos - rad;
      vec3r posMax = pos + rad;
      if (minBox.x > posMin.x) minBox.x = posMin.x;
      if (maxBox.x < posMax.x) maxBox.x = posMax.x;
      if (minBox.y > posMin.y) minBox.y = posMin.y;
      if (maxBox.y < posMax.y) maxBox.y = posMax.y;
      if (minBox.z > posMin.z) minBox.z = posMin.z;
      if (maxBox.z < posMax.z) maxBox.z = posMax.z;
    }
  }

  std::cout << "Export sample " << std::endl;

  std::ofstream file("exportedSample.txt");
  file << "periodicity " << maxBox.x - minBox.x << ' ' << maxBox.y - minBox.y << ' ' << maxBox.z - minBox.z << '\n';
  file << "Particles " << Shapes.size() << '\n';
  for (size_t i = 0; i < Shapes.size(); i++) {
    file << Shapes[i].name << " 1 0 1  " << Shapes[i].position - minBox << "  0 0 0  0 0 0  " << Shapes[i].orientation
         << "  0 0 0  0 0 0\n";
  }
}

void menu(int num) {
  switch (num) {
    case 0:
      exit(0);
      break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  glutCreateMenu(menu);  // Main menu

  glutAddMenuEntry("Export release config. for DEMbox", 1);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {

  StackTracer::initSignals();

  if (argc == 1) {
    if (readShapeLib("shapes") == 0) {
      return 0;
    }
  } else if (argc == 2) {
    if (readShapeLib(argv[1]) == 0) {
      return 0;
    }
  }

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutSetOption(GLUT_MULTISAMPLE, 8);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH | GLUT_MULTISAMPLE);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("ShapeSurvey");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  // glutSpecialFunc(processSpecialKeys);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  glText::init();

  // ==== Init the visualizer
  center.set(0.0, 0.0, 0.0);  // where we look at
  eye.set(0.0, 0.0, 1.0);     // from where we look
  up.set(0.0, 1.0, 0.0);      // direction (normalized)

  mouse_mode = NOTHING;
  view_angle = 45.0;
  znear = 0.01;
  zfar = 10.0;

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_COLOR_MATERIAL);

  // Create light components
  GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
  GLfloat diffuseLight[] = {0.8f, 0.8f, 0.8, 1.0f};
  GLfloat specularLight[] = {0.5f, 0.5f, 0.5f, 1.0f};
  GLfloat positionLight0[] = {1000000.0f, 1000000.0f, 1000000.0f, 1.0f};
  GLfloat positionLight1[] = {-1000000.0f, -1000000.0f, -1000000.0f, 1.0f};

  // Assign created components to GL_LIGHT0
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, positionLight0);

  // Assign created components to GL_LIGHT1
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT1, GL_POSITION, positionLight1);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  // ==== Enter GLUT event processing cycle
  adjust_clipping_plans();
  fit_view();
  glutMainLoop();
  return 0;
}
