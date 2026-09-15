#include "shapeSurvey.hpp"

#include "fileTool.hpp"
#include "glTools.hpp"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <vector>

#ifdef __APPLE__
#include <OpenGL/glext.h>
#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "toofus-gate/stb/stb_image_write.h"

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

// Draw the given geometry again as a thin black wireframe sitting on the surface
// (a small polygon offset keeps the lines from z-fighting with the fill).
static void wireOverlay(const std::function<void()>& emit) {
  GLboolean lit = glIsEnabled(GL_LIGHTING);
  glDisable(GL_LIGHTING);
  glColor3f(0.0f, 0.0f, 0.0f);
  glLineWidth(wire_width);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glEnable(GL_POLYGON_OFFSET_LINE);
  glPolygonOffset(-1.0f, -1.0f);
  emit();
  glDisable(GL_POLYGON_OFFSET_LINE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  if (lit) glEnable(GL_LIGHTING);
}

// Area-weighted centroid of a planar face, in 3D. A plain vertex average can
// land on a reentrant corner (or outside the polygon) for a non-convex face,
// which would make the inside() probe below unreliable.
static vec3r faceCentroid(size_t ishp, const std::vector<size_t>& F, const vec3r& N) {
  const std::vector<vec3r>& V = Shapes[ishp].vertex;
  vec3r o = V[F[0]];
  vec3r bx = V[F[1]] - o; bx.normalize();
  vec3r by = cross(N, bx);
  double A = 0.0, cx = 0.0, cy = 0.0;
  size_t n = F.size();
  for (size_t k = 0; k < n; ++k) {
    vec3r p0 = V[F[k]] - o, p1 = V[F[(k + 1) % n]] - o;
    double x0 = p0 * bx, y0 = p0 * by, x1 = p1 * bx, y1 = p1 * by;
    double cr = x0 * y1 - x1 * y0;
    A += cr; cx += (x0 + x1) * cr; cy += (y0 + y1) * cr;
  }
  if (std::fabs(A) < 1e-14) {                 // degenerate: fall back to the average
    vec3r c(0, 0, 0);
    for (size_t k = 0; k < n; ++k) c += V[F[k]];
    return c / (double)n;
  }
  return o + (cx / (3.0 * A)) * bx + (cy / (3.0 * A)) * by;
}

// Emit the overlapping-primitive geometry of a shape (no colour set here, so the
// same call serves both the filled pass and the wireframe pass).
static void emitPrimitives(size_t ishp) {
  double R = Shapes[ishp].radius;

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

  bool solid = !Shapes[ishp].isSurface;
  size_t nf = Shapes[ishp].face.size();
  for (size_t f = 0; f < nf; ++f) {
    const std::vector<size_t>& F = Shapes[ishp].face[f];
    if (F.size() < 3) {
      continue;
    }  // At least 3 pts!
    vec3r N = cross(Shapes[ishp].vertex[F[1]] - Shapes[ishp].vertex[F[0]],
                    Shapes[ishp].vertex[F[2]] - Shapes[ishp].vertex[F[0]]);
    N.normalize();

    // Which offset polygons are actually on the surface?
    //   open surface -> both (+R and -R): the shape is a slab of thickness 2R.
    //   solid        -> only the outward one; the inward offset is an interior
    //                   face (hidden when opaque, but it shows through in
    //                   transparency, so we must not draw it).
    int sLo = -1, sHi = +1;                    // default: draw both (surface)
    if (solid) {
      vec3r cf = faceCentroid(ishp, F, N);
      double m = 1.25 * R;                      // probe just beyond the surface
      bool inPlus  = Shapes[ishp].inside(cf + N * m);
      bool inMinus = Shapes[ishp].inside(cf - N * m);
      if (inPlus != inMinus) {                  // outward = the side that is outside
        sLo = sHi = (!inPlus) ? +1 : -1;
      }                                         // else ambiguous: keep both (safe)
    }

    for (int s = sLo; s <= sHi; s += 2) {
      vec3r Ns = (s > 0) ? N : -N;
      glBegin(GL_TRIANGLE_FAN);
      glNormal3f(Ns.x, Ns.y, Ns.z);
      for (size_t v = 0; v < F.size(); ++v) {
        glVertex3f(Shapes[ishp].vertex[F[v]].x + Ns.x * R,
                   Shapes[ishp].vertex[F[v]].y + Ns.y * R,
                   Shapes[ishp].vertex[F[v]].z + Ns.z * R);
      }
      glEnd();
    }
  }
}

void drawShape(size_t ishp) {
  if (ishp >= Shapes.size()) {
    return;
  }
  if (mouse_mode != NOTHING) {
    return;
  }

  glColor4f(shapeColor[0], shapeColor[1], shapeColor[2], alpha);
  emitPrimitives(ishp);
  if (show_wire) {
    wireOverlay([ishp]() { emitPrimitives(ishp); });
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

// Draw the skin mesh of a shape (the .rmsh companion), as an alternative to the
// overlapping-primitive rendering of drawShape(). Positions and exact normals
// are stored per vertex, in the body frame, exactly like drawShape() uses.
void drawSkin(size_t ishp) {
  if (ishp >= Shapes.size()) return;
  auto it = skinMeshes.find(Shapes[ishp].name);
  if (it == skinMeshes.end()) {
    std::cout << "No skin mesh named '" << Shapes[ishp].name << "' in " << skinFile << std::endl;
    return;
  }
  const ShapeMesh& m = it->second;

  // Triangle draw order. Marching cubes emits its triangles grid-slab by
  // grid-slab, and blending a transparent mesh in that order against a written
  // depth buffer makes the far side appear and vanish in slab-wide bands, which
  // read as internal planes. Sorting the triangles back to front (the same cure
  // as in `see`) blends them correctly and still leaves a valid depth buffer
  // for the wireframe overlay. Opaque rendering needs no order.
  static std::vector<size_t> order;
  order.resize(m.tri.size());
  for (size_t t = 0; t < order.size(); ++t) order[t] = t;
  if (alpha < 0.999f) {
    std::vector<double> d2(m.tri.size());
    for (size_t t = 0; t < m.tri.size(); ++t) {
      vec3r c = (1.0 / 3.0) * (m.P[m.tri[t].a] + m.P[m.tri[t].b] + m.P[m.tri[t].c]);
      d2[t] = norm2(c - eye);
    }
    std::sort(order.begin(), order.end(), [&d2](size_t a, size_t b) { return d2[a] > d2[b]; });
  }

  auto emit = [&m]() {
    glBegin(GL_TRIANGLES);
    for (size_t o = 0; o < order.size(); ++o) {
      size_t t = order[o];
      size_t idx[3] = {m.tri[t].a, m.tri[t].b, m.tri[t].c};
      for (int k = 0; k < 3; ++k) {
        const vec3r& n = m.N[idx[k]];
        const vec3r& p = m.P[idx[k]];
        glNormal3f(n.x, n.y, n.z);
        glVertex3f(p.x, p.y, p.z);
      }
    }
    glEnd();
  };
  glColor4f(shapeColor[0], shapeColor[1], shapeColor[2], alpha);
  emit();
  if (show_wire) {
    wireOverlay(emit);
  }
}

// Place the camera on a sphere around the current shape, using an elevation and
// azimuth (degrees, z up) like matplotlib's mplot3d, at the fit-view distance.
void set_view(double elev_deg, double azim_deg) {
  OBB& obb = Shapes[ishape].obb;
  center = obb.center;
  vec3r diag = 2.0 * (obb.extent[0] * obb.e[0] + obb.extent[1] * obb.e[1] + obb.extent[2] * obb.e[2]);
  double d = 0.5 * diag.length() / (atan(view_angle * M_PI / 360.0));
  double e = elev_deg * M_PI / 180.0, a = azim_deg * M_PI / 180.0;
  vec3r dir(cos(e) * cos(a), cos(e) * sin(a), sin(e));
  eye = center + d * dir;
  vec3r worldUp(0.0, 0.0, 1.0);
  if (fabs(dir * worldUp) > 0.98) worldUp.set(0.0, 1.0, 0.0);
  up = worldUp;
}

// Read the rendered frame (W x H pixels of the current read buffer) back and
// write it as a PNG (row-flipped).
void saveScreenshot(const char* filename, int W, int H) {
  std::vector<unsigned char> pix((size_t)W * H * 4);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, W, H, GL_RGBA, GL_UNSIGNED_BYTE, pix.data());
  std::vector<unsigned char> flip((size_t)W * H * 4);
  for (int y = 0; y < H; ++y) {
    std::memcpy(&flip[(size_t)(H - 1 - y) * W * 4], &pix[(size_t)y * W * 4], (size_t)W * 4);
  }
  if (stbi_write_png(filename, W, H, 4, flip.data(), W * 4) == 0) {
    std::cout << "Could not write " << filename << std::endl;
  } else {
    std::cout << "wrote " << filename << " (" << W << "x" << H << ")" << std::endl;
  }
}

#if defined(GL_EXT_framebuffer_object) && defined(GL_EXT_framebuffer_multisample) && \
    defined(GL_EXT_framebuffer_blit)
#define SHOT_OFFSCREEN 1

// Off-screen render target of the screenshot: the frame is drawn into a
// multisampled framebuffer object of the requested size, then resolved into a
// plain one that is read back. The image size is thus not limited by the size of
// the window, which the system clamps to the screen.
struct OffscreenTarget {
  GLuint msFbo{0}, msColor{0}, msDepth{0};  // multisampled, drawn into
  GLuint fbo{0}, color{0};                  // single-sampled, read back
};

static GLuint attachRenderbuffer(GLenum attachment, GLenum format, GLsizei samples, int W, int H) {
  GLuint rb = 0;
  glGenRenderbuffersEXT(1, &rb);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rb);
  if (samples > 0) {
    glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, samples, format, W, H);
  } else {
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, format, W, H);
  }
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, attachment, GL_RENDERBUFFER_EXT, rb);
  return rb;
}

// Returns false when the driver cannot provide the buffers; the caller then
// draws into the window instead
static bool createOffscreen(OffscreenTarget& t, int W, int H) {
  GLint maxSamples = 0;
  glGetIntegerv(GL_MAX_SAMPLES_EXT, &maxSamples);
  GLsizei samples = std::min(8, (int)maxSamples);

  glGenFramebuffersEXT(1, &t.msFbo);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, t.msFbo);
  t.msColor = attachRenderbuffer(GL_COLOR_ATTACHMENT0_EXT, GL_RGBA8, samples, W, H);
  t.msDepth = attachRenderbuffer(GL_DEPTH_ATTACHMENT_EXT, GL_DEPTH_COMPONENT24, samples, W, H);
  bool ok = (glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT) == GL_FRAMEBUFFER_COMPLETE_EXT);

  glGenFramebuffersEXT(1, &t.fbo);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, t.fbo);
  t.color = attachRenderbuffer(GL_COLOR_ATTACHMENT0_EXT, GL_RGBA8, 0, W, H);
  ok = ok && (glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT) == GL_FRAMEBUFFER_COMPLETE_EXT);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, ok ? t.msFbo : 0);
  return ok;
}

// Copy the multisampled frame into the plain buffer and make it the read buffer
static void resolveOffscreen(const OffscreenTarget& t, int W, int H) {
  glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, t.msFbo);
  glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, t.fbo);
  glBlitFramebufferEXT(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_NEAREST);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, t.fbo);
  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
}
#endif

// One-shot display callback used in screenshot mode: render a clean frame (just
// the shape, no axes/OBB/text), save it, and quit.
void screenshotDisplay() {
#ifdef SHOT_OFFSCREEN
  OffscreenTarget target;
  bool offscreen = createOffscreen(target, shot_width, shot_height);
#else
  bool offscreen = false;
#endif
  if (offscreen) {
    // Same aspect ratio and viewport as a window of the requested size
    width = shot_width;
    height = shot_height;
    glViewport(0, 0, width, height);
    glEnable(GL_MULTISAMPLE);
  } else {
    std::cout << "No off-screen buffer: the image is limited to the window size" << std::endl;
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    width = vp[2];
    height = vp[3];
  }

  if (shot_transparent) {
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
  } else {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  }
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  adjust_clipping_plans();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);

  if (!skinFile.empty()) {
    drawSkin(ishape);
  } else {
    drawShape(ishape);
  }

  glFinish();
#ifdef SHOT_OFFSCREEN
  if (offscreen) {
    resolveOffscreen(target, width, height);
  } else {
    glReadBuffer(GL_BACK);
  }
#else
  glReadBuffer(GL_BACK);
#endif
  saveScreenshot(outputFile.c_str(), width, height);
  exit(0);
}

void printCLIhelp() {
  std::cout << "shapeSurvey - inspect r-shapes, or render one to an image (CLI)\n\n"
            << "Usage: shapeSurvey <shapeFile> [options]\n\n"
            << "  -o, --output <file.png>  render one clean frame to <file.png> and exit\n"
            << "  -i, --index  <n>         shape to render (0-based, default 0)\n"
            << "  -N, --name   <name>      shape to render, by name\n"
            << "  -W, --width  <px>        image width  (default 800)\n"
            << "  -H, --height <px>        image height (default 800)\n"
            << "  -A, --alpha  <a>         opacity in [0,1] (default 1)\n"
            << "      --elev   <deg>       camera elevation, z up (default 22)\n"
            << "      --azim   <deg>       camera azimuth        (default -55)\n"
            << "      --skin   <file.rmsh> draw the skin mesh instead of the primitives\n"
            << "      --radius <R>         override the Minkowski radius (0 = raw triangles)\n"
            << "      --wire               overlay the mesh as a thin black wireframe\n"
            << "      --line-width <px>    line width of that wireframe (default 0.7)\n"
            << "      --color  <r,g,b>     fill colour in [0,1] (default 0.761,0.733,0.976)\n"
            << "      --transparent        transparent background\n"
            << "  -h, --help               this help\n\n"
            << "With no --output, shapeSurvey opens the interactive viewer.\n";
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {

  StackTracer::initSignals();

  // ---- Parse the command line (positional shapeFile + options) ----
  std::string inputFile = "shapes";
  int shapeIndex = 0;
  std::string shapeName = "";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&](const char* def) -> std::string {
      return (i + 1 < argc) ? std::string(argv[++i]) : std::string(def);
    };
    if (a == "-o" || a == "--output") outputFile = next("out.png");
    else if (a == "-i" || a == "--index") shapeIndex = std::stoi(next("0"));
    else if (a == "-N" || a == "--name") shapeName = next("");
    else if (a == "-W" || a == "--width") width = std::stoi(next("800"));
    else if (a == "-H" || a == "--height") height = std::stoi(next("800"));
    else if (a == "-A" || a == "--alpha") alpha = std::stof(next("1"));
    else if (a == "--elev") shot_elev = std::stod(next("22"));
    else if (a == "--azim") shot_azim = std::stod(next("-55"));
    else if (a == "--skin") skinFile = next("");
    else if (a == "--radius") radiusOverride = std::stod(next("-1"));
    else if (a == "--transparent") shot_transparent = 1;
    else if (a == "--wire") show_wire = 1;
    else if (a == "--line-width") wire_width = std::stof(next("0.7"));
    else if (a == "--color") {
      std::string c = next("0.761,0.733,0.976");
      std::replace(c.begin(), c.end(), ',', ' ');
      std::istringstream iss(c);
      iss >> shapeColor[0] >> shapeColor[1] >> shapeColor[2];
    }
    else if (a == "--bg") show_background = 1;
    else if (a == "-h" || a == "--help") { printCLIhelp(); return 0; }
    else if (!a.empty() && a[0] != '-') inputFile = a;
    else { std::cout << "Unknown option: " << a << std::endl; printCLIhelp(); return 0; }
  }

  // The window can be clamped to the screen (and its size is then overwritten by
  // reshape), so the screenshot keeps the size that was asked for
  shot_width = width;
  shot_height = height;

  if (readShapeLib(inputFile.c_str()) == 0) {
    return 0;
  }

  // ---- Select the shape to show/render ----
  if (!shapeName.empty()) {
    bool found = false;
    for (size_t i = 0; i < Shapes.size(); ++i) {
      if (Shapes[i].name == shapeName) { ishape = i; found = true; break; }
    }
    if (!found) std::cout << "No shape named '" << shapeName << "', using shape 0" << std::endl;
  } else if (shapeIndex >= 0 && (size_t)shapeIndex < Shapes.size()) {
    ishape = (size_t)shapeIndex;
  }
  if (radiusOverride >= 0.0) {  // e.g. --radius 0 shows the imported triangles, no rounding
    Shapes[ishape].radius = radiusOverride;
  }
  if (Shapes[ishape].preCompDone == 'n') Shapes[ishape].fitObb();
  if (!skinFile.empty()) {
    skinMeshes = loadRmsh(skinFile);
    std::cout << "Loaded " << skinMeshes.size() << " skin mesh(es) from " << skinFile << std::endl;
  }

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  //glutSetOption(GLUT_MULTISAMPLE, 8);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH | GLUT_MULTISAMPLE);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("ShapeSurvey");

  bool screenshotMode = !outputFile.empty();

  // ==== Register callbacks
  if (screenshotMode) {
    glutDisplayFunc(screenshotDisplay);  // renders one frame, saves, exits
    glutReshapeFunc(reshape);
  } else {
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    // glutSpecialFunc(processSpecialKeys);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);

    // ==== Menu
    buildMenu();
    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }

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
  if (screenshotMode) {
    set_view(shot_elev, shot_azim);
  } else {
    fit_view();
  }
  adjust_clipping_plans();
  glutMainLoop();
  return 0;
}
