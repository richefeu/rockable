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

#include "see3.hpp"

void showKeybinds() {
  switch2D::go(width, height);

  glColor4f(1.0f, 1.0f, 1.0f, 0.6f);
  glBegin(GL_QUADS);
  int nbLines = 21;  // !!!!! increment this when a line is added !!!!!  <=====
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
  glText::print(15, _nextLine_, "[ESC]   quit");
  glText::print(15, _nextLine_, "[RIGHT] load next configuration file");
  glText::print(15, _nextLine_, "[LEFT]  load previous configuration file");
  glText::print(15, _nextLine_, "[=]     fit the view");
  // glText::print(15, _nextLine_, "[a][A]  decrease/increase alpha (transparence) of driven bodies");
  glText::print(15, _nextLine_, "[b]     switch ON/OFF the background color");
  // glText::print(15, _nextLine_, "[c]     run 5000 steps of computation (used for debugging)");
  // glText::print(15, _nextLine_, "[e][E]  decrease/increase alpha (transparence) of free bodies");
  // glText::print(15, _nextLine_, "[g]     open another file");
  glText::print(15, _nextLine_, "[k]     print this help");
  glText::print(15, _nextLine_, "[l]     switch ON/OFF the links (normal vector at contact)");
  glText::print(15, _nextLine_, "[m]     switch ON/OFF the links colored by type");
  // glText::print(15, _nextLine_, "[n]     switch ON/OFF the body displays");
  // glText::print(15, _nextLine_, "[o]     switch ON/OFF the OBB displays");
  // glText::print(15, _nextLine_, "[O]     switch ON/OFF the enlargement of OBBs by the Verlet distance");
  // glText::print(15, _nextLine_, "[w]     set the view so that gravity appears vertical");
  // glText::print(15, _nextLine_, "[x]     print the space limits of the current scene");
  // glText::print(15, _nextLine_, "[p]     edit selected body");

  glText::print(15, _nextLine_, "[y]     make the display faster (and less nice)");
#ifdef PNG_H
  glText::print(15, _nextLine_, "[z]     make a screenshot (oneshot.png)");
  glText::print(15, _nextLine_, "[Z]     make a series of screenshots (shotX.png)");
#else
  glText::print(15, _nextLine_, "[z]     make a screenshot (oneshot.tga)");
  glText::print(15, _nextLine_, "[Z]     make a series of screenshots (shotX.tga)");
#endif

#undef _nextLine_
  switch2D::back();
}

void keyboard(GLFWwindow* /*window*/, int key, int /*scancode*/, int action, int mods) {
  if (ImGui_window_focused == true) return;
  if (action == GLFW_PRESS) {
    // std::cout << "> " << glfwGetKeyName(key, 0) << '\n';
    switch (key) {
      case GLFW_KEY_ESCAPE:
        request_quit = true;
        break;

        /*
      case '!': {
        show_probe = 1 - show_probe;
        if (show_probe == 1 && fileTool::fileExists("probe.txt")) {
          std::ifstream f("probe.txt");
          f >> probe.min >> probe.max >> probe_MCnsteps;
          std::cout << "Probe: min = " << probe.min << "\n       max = " << probe.max << "\n";
        }
      } break;

      case '@': {
        std::cout << "Probing with " << probe_MCnsteps << " Monte-Carlo steps...\n" << std::flush;
        double SF = box.probeSolidFraction(probe, probe_MCnsteps);
        std::cout << "Solid Fraction in probe: " << SF << '\n';
      } break;
      */

      case GLFW_KEY_SPACE:
        textZone.reset();
        break;

      case GLFW_KEY_B:
        show_background = 1 - show_background;
        break;

        /*
      case GLFW_KEY_C: {  // compute a few steps (for debugging)
        box.UpdateNL();
        for (int i = 0; i < 5000; i++) box.velocityVerletStep();
        textZone.addLine("5000 time steps have been done.");
      } break;
*/

      case GLFW_KEY_D:
        show_driven = 1 - show_driven;
        break;

        /*
      case GLFW_KEY_E: {
        if (mods == GLFW_MOD_SHIFT && alpha_particles <= 0.95)
          alpha_particles += 0.05;
        else if (alpha_particles > 0.1)
          alpha_particles -= 0.05;
      }
        */

        break;

      case GLFW_KEY_F:
        show_forces = 1 - show_forces;
        break;

        /*
      case GLFW_KEY_G: {
        std::cout << ">>>>>>>>>> CONF NUMBER: ";
        int Num;
        std::cin >> Num;
        tryToReadConf(Num);
      } break;
*/

      case GLFW_KEY_H: {
        textZone.addLine("");
        textZone.addLine("MOUSE:");
        textZone.addLine("        Left button + move = ROTATE");
        textZone.addLine("SHIFT + Left button + move = PAN");
        textZone.addLine("CTRL  + Left button + move = ZOOM (or middle button)");
        textZone.addLine("ALT   + Left button        = SELECT");
        textZone.addLine("The KEYBOARD shortcuts -> key [k]");
      } break;

      case GLFW_KEY_K:
        show_keybinds = 1 - show_keybinds;
        break;

      case GLFW_KEY_L:
        show_interFrames = 1 - show_interFrames;
        break;

      case GLFW_KEY_M:
        show_interTypes = 1 - show_interTypes;
        break;

      case GLFW_KEY_R: {
        rescaleColorRange = 1 - rescaleColorRange;
        if (rescaleColorRange == 0) {
          std::cout << "rescaleColorRange = 0\n";
          std::cout << "colorRangeMin = " << colorRangeMin << '\n';
          std::cout << "colorRangeMax = " << colorRangeMax << '\n';
        } else if (rescaleColorRange == 1) {
          std::cout << "rescaleColorRange = 1\n";
        }
      } break;

      case GLFW_KEY_S:
        show_slice = 1 - show_slice;
        break;

      case GLFW_KEY_T: {
        if (mods == GLFW_MOD_SHIFT) {
          if (forceTubeFactor <= 0.9) forceTubeFactor += 0.05;
        } else {
          if (forceTubeFactor > 0.2) forceTubeFactor -= 0.05;
        }
      } break;

      case GLFW_KEY_Z: {
        if (mods == GLFW_MOD_SHIFT) {
          // be carreful there's no way to stop this loop
          // if the process if too long
          while (tryToReadConf(confNum + 1)) {
            char name[256];
#ifdef PNG_H
            snprintf(name, 256, "shot%d.png", confNum);
#else
            snprintf(name, 256, "shot%d.tga", confNum);
#endif
            display();
            screenshot(name);
          }
        } else {
#ifdef PNG_H
          screenshot("oneshot.png");
#else
          screenshot("oneshot.tga");
#endif
        }
      } break;

      case GLFW_KEY_UP:
        textZone.increase_nbLine();
        break;

      case GLFW_KEY_DOWN:
        textZone.decrease_nbLine();
        break;

      case GLFW_KEY_LEFT:
        if (confNum > 0) tryToReadConf(confNum - 1);
        break;

      case GLFW_KEY_RIGHT:
        tryToReadConf(confNum + 1);
        break;

      case GLFW_KEY_0: {
        colorMode = 0;
        pcolors.clear();
      } break;

      case GLFW_KEY_1: {
        colorMode = 1;
        pcolors.clear();
        colorRGBA col;
        for (size_t i = 0; i < box.Particles.size(); ++i) {
          CT.getCyclicRGB8(&col);
          pcolors.push_back(col);
        }
      } break;

      case GLFW_KEY_2: {
        colorMode = 2;
        resetColors(colorMode, 1);
      } break;
    }
  }
}

void resetColors(int mode, int rescale) {
  switch (mode) {
    case 2: {  // velocity magnitude
      if (rescale == 1) {
        colorRangeMin = colorRangeMax = 0.0;
        for (size_t i = box.nDriven; i < box.Particles.size(); ++i) {
          double v = norm(box.Particles[i].vel);
          if (v > colorRangeMax) colorRangeMax = v;
        }
        std::cout << "colorRangeMin = " << colorRangeMin << '\n';
        std::cout << "colorRangeMax = " << colorRangeMax << '\n';
      }

      CT.setTableID(MATLAB_HOT);
      CT.setMinMax(colorRangeMin, colorRangeMax);
      CT.Rebuild();
      pcolors.clear();
      colorRGBA col;
      for (size_t i = 0; i < box.Particles.size(); ++i) {
        double v = norm(box.Particles[i].vel);
        CT.getRGB(v, &col);
        pcolors.push_back(col);
      }
    } break;
  }
}

#define BUFSIZE 1024
void selection(int x, int y) {
  GLuint selectBuf[BUFSIZE];
  GLint hits;
  GLint viewport[4];

  selectedParticle = -1;

  glGetIntegerv(GL_VIEWPORT, viewport);  // Get viewport data
  glSelectBuffer(BUFSIZE, selectBuf);    // select the buffer
  glRenderMode(GL_SELECT);

  // The following is quasi the same as adjustClippingPlans
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPickMatrix(x, viewport[3] - y, 2.0, 2.0, viewport);

  wh_ratio = (float)width / (float)height;
  computePerspective();
  gluPerspective(view_angle, wh_ratio, znear, zfar);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

  glEnable(GL_DEPTH_TEST);
  glInitNames();
  glPushName(0);

  if (show_particles) {
    glColor3f(0.5f, 0.5f, 0.5f);
    size_t i0 = 0;
    if (show_driven == 0) i0 = box.nDriven;
    for (size_t i = i0; i < box.Particles.size(); ++i) {
      glLoadName(i);
      vec3r pos = box.Particles[i].pos;

      glPushMatrix();
      glTranslatef(pos.x, pos.y, pos.z);
      quat2GLMatrix<GLfloat>(box.Particles[i].Q, Rot_Matrix);
      glMultMatrixf(Rot_Matrix);
      drawShape(box.Particles[i].shape, box.Particles[i].homothety);
      glPopMatrix();
    }
  }

  glFlush();

  hits = glRenderMode(GL_RENDER);
  // processHits:
  if (hits == 0) {
    selectedParticle = -1;
    return;
  }

  GLuint* ptr = (GLuint*)selectBuf;

  GLuint names = *ptr;
  ptr++;
  GLuint minZ = *ptr;
  GLuint* ptrName = ptr + 2;
  ptr += names + 2;

  for (GLint i = 1; i < hits; i++) {
    names = *ptr;
    ptr++;
    if (*ptr < minZ) {
      minZ = *ptr;
      ptrName = ptr + 2;
    }

    ptr += names + 2;
  }

  selectedParticle = (int)*ptrName;  // Only 1 name since we used only 1
                                     // 'pushName' followed by 'loadName' for
                                     // each master-body
}
#undef BUFSIZE

void mouse(GLFWwindow* window, int button, int action, int mods) {
  if (ImGui_window_focused == true) return;

  if (action == GLFW_RELEASE) {
    mouse_mode = NOTHING;
  } else if (action == GLFW_PRESS) {
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    mouse_start[0] = (int)floor(xpos);
    mouse_start[1] = (int)floor(ypos);
    switch (button) {
      case GLFW_MOUSE_BUTTON_RIGHT:
        mouse_mode = PAN;
        break;
      case GLFW_MOUSE_BUTTON_LEFT:
        if (mods == GLFW_MOD_SHIFT)
          mouse_mode = PAN;
        else if (mods == GLFW_MOD_ALT)
          selection(mouse_start[0], mouse_start[1]);
        else if (mods == GLFW_MOD_CONTROL)
          mouse_mode = ZOOM;
        else
          mouse_mode = ROTATION;
        break;
      case GLFW_MOUSE_BUTTON_MIDDLE:
        mouse_mode = ZOOM;
        break;
    }
  }
}

void motion(GLFWwindow* /*window*/, double x, double y) {
  if (ImGui_window_focused == true) return;
  if (mouse_mode == NOTHING) return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;
  double length;
  vec3r axis;

  switch (mouse_mode) {

    case ROTATION:
      axis = (cross(up, center - eye));
      axis.normalize();
      eye = geoTool::rotatePoint(eye, center, up, -dx * M_PI);
      eye = geoTool::rotatePoint(eye, center, axis, dy * M_PI);
      up = (geoTool::rotatePoint((center + up), center, axis, dy * M_PI) - center);
      up.normalize();
      break;

    case ZOOM:
      eye = center + (eye - center) * (dy + 1.0);
      break;

    case PAN:
      length = (eye - center).length() * tan(view_angle * M_PI / 360.0) * 2.0;
      axis = cross(up, center - eye);
      axis.normalize();
      center = center + axis * dx * length * 0.8;
      center = center + up * dy * length;
      break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;
}

void display() {
  sleep(0);  // it is supposed to accelerate the display
  glTools::clearBackground(show_background, (int)floor(255 * clear_bottom_color[0]),
                           (int)floor(255 * clear_bottom_color[1]), (int)floor(255 * clear_bottom_color[2]),
                           (int)floor(255 * clear_top_color[0]), (int)floor(255 * clear_top_color[1]),
                           (int)floor(255 * clear_top_color[2]));
  adjustClippingPlans();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  if (show_globalFrame) drawGlobalFrame();
  if (show_traj) drawTrajectories();
  if (show_forces) drawForces();
  if (show_interFrames) drawInteractionFrames();
  if (show_interTypes) drawInteractionTypes();
  if (show_obb) drawOBBs();
  if (show_particles) drawParticles();
  if (show_probe) drawProbe();

  if (show_keybinds) showKeybinds();
  textZone.draw();
  glFlush();
}

void computePerspective() {
  double zf = (eye - center).normalize();

  vec3r mx = (box.aabb.max - box.aabb.min);
  max_length = (GLfloat)(2.0 * norm(mx));

  znear = zf - 0.5 * max_length;
  double close_dst = 0.1 * zf;
  if (znear < close_dst) znear = close_dst;
  zfar = zf + 0.5 * max_length;
}

void adjustClippingPlans() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  wh_ratio = (float)width / (float)height;
  computePerspective();
  gluPerspective(view_angle, wh_ratio, znear, zfar);
  glMatrixMode(GL_MODELVIEW);
}

void fitView() {
  vec3r dir = (eye - center);
  vec3r diag = (box.aabb.max - box.aabb.min);
  dir.normalize();
  center = 0.5 * (box.aabb.max + box.aabb.min);
  GLfloat d = 0.5 * diag.length() / (atan(view_angle * M_PI / 360.0));
  eye = center + d * dir;
}

void reshape(GLFWwindow* /*window*/, int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, width, height);

  adjustClippingPlans();
}

// Draw the shape of the sphero-polyhedron in its own framework
void drawShape(Shape* s, double homothety) {
  double R = homothety * s->radius;
  int nbLevelSphere = 2;
  if (totalNumberOfVertices > 10000) {
    nbLevelSphere = 1;
  }

  if (shapeWithoutThickness == 0) {
    // vertixes (spheres)
    for (size_t v = 0; v < s->vertex.size(); ++v) {
      glPushMatrix();
      glTranslatef(homothety * s->vertex[v].x, homothety * s->vertex[v].y, homothety * s->vertex[v].z);
      facetSphere::draw(nbLevelSphere, R);
      glPopMatrix();
    }

    // edges (tubes)
    for (size_t e = 0; e < s->edge.size(); ++e) {
      vec3r orig = homothety * s->vertex[s->edge[e].first];
      vec3r arrow = homothety * s->vertex[s->edge[e].second];
      arrow -= orig;
      glShape::tube(orig, arrow, 2.0 * R);
    }
  }

  // faces (3D polygones)
  for (size_t f = 0; f < s->face.size(); ++f) {
    if (s->face[f].size() < 3) continue;  // At least 3 pts!
    vec3r N =
        cross(s->vertex[s->face[f][1]] - s->vertex[s->face[f][0]], s->vertex[s->face[f][2]] - s->vertex[s->face[f][0]]);
    N.normalize();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(N.x, N.y, N.z);
    for (size_t v = 0; v < s->face[f].size(); ++v) {
      glVertex3f(homothety * s->vertex[s->face[f][v]].x + N.x * R, homothety * s->vertex[s->face[f][v]].y + N.y * R,
                 homothety * s->vertex[s->face[f][v]].z + N.z * R);
    }
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(-N.x, -N.y, -N.z);
    for (size_t v = 0; v < s->face[f].size(); ++v) {
      glVertex3f(homothety * s->vertex[s->face[f][v]].x - N.x * R, homothety * s->vertex[s->face[f][v]].y - N.y * R,
                 homothety * s->vertex[s->face[f][v]].z - N.z * R);
    }
    glEnd();
  }
}

void drawParticles() {
  if (mouse_mode != NOTHING && totalNumberOfVertices > 100) {
    drawOBBs();
    return;
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  for (size_t i = box.nDriven; i < box.Particles.size(); ++i) {
    if (selectedParticle >= 0 && i == (size_t)selectedParticle) {
      glColor4ub(234, 255, 0, (int)floor(alpha_particles * 255));  // yellow
    } else {
      if (i >= pcolors.size()) {
        glColor4ub(178, 34, 34, (int)floor(alpha_particles * 255));  // brick color
      } else {
        glColor4ub(pcolors[i].r, pcolors[i].g, pcolors[i].b, (int)floor(alpha_particles * 255));
      }
    }

    vec3r pos = box.Particles[i].pos;

    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);
    quat2GLMatrix<GLfloat>(box.Particles[i].Q, Rot_Matrix);
    glMultMatrixf(Rot_Matrix);
    drawShape(box.Particles[i].shape, box.Particles[i].homothety);
    glPopMatrix();
  }

  if (show_driven == 1) {
    for (size_t i = 0; i < box.nDriven; ++i) {
      if (selectedParticle >= 0 && i == (size_t)selectedParticle) {
        glColor4ub(234, 255, 0, (int)floor(alpha_fixparticles * 255));
      } else
        glColor4ub(128, 128, 128, (int)floor(alpha_fixparticles * 255));

      vec3r pos = box.Particles[i].pos;

      glPushMatrix();
      glTranslatef(pos.x, pos.y, pos.z);
      quat2GLMatrix<GLfloat>(box.Particles[i].Q, Rot_Matrix);
      glMultMatrixf(Rot_Matrix);
      drawShape(box.Particles[i].shape, box.Particles[i].homothety);
      glPopMatrix();
    }
  }
}

void drawTrajectories() {
  if (mouse_mode != NOTHING && releases.size() > 200) return;

  glDisable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glLineWidth(2.0f);

  for (size_t r = 0; r < releases.size(); ++r) {
    for (size_t f = 0; f < releases[r].freeFlights.size(); ++f) {

      glBegin(GL_LINE_STRIP);
      vec3r p;
      for (double s = 0.0; s <= 1.0; s += 0.1) {
        double t = s * releases[r].freeFlights[f].duration;
        p = releases[r].freeFlights[f].pos + releases[r].freeFlights[f].vel * t + 0.5 * box.gravity * t * t;
        glVertex3f(p.x, p.y, p.z);
      }
      glEnd();
    }
  }
}

void drawGlobalFrame() {
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);

  vec3r diag = 0.6 * (box.aabb.max - box.aabb.min);
  glShape::frame(vec3r(0, 0, 0), diag.x, diag.y, diag.z);
}

void drawOBBs() {
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glLineWidth(1.0f);

  OBB obbi, obbj;
  for (size_t i = 0; i < box.Particles.size(); i++) {

    obbi = box.Particles[i].shape->obb;
    obbi.rotate(box.Particles[i].Q);
    obbi.extent *= box.Particles[i].homothety;
    obbi.center *= box.Particles[i].homothety;
    obbi.center += box.Particles[i].pos;
    if (enlarged_obb) obbi.enlarge(0.5 * box.DVerlet);

    glColor4ub(255, 0, 0, 255);  // RED
    glShape::obb(obbi);
  }
}

void drawProbe() {
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glLineWidth(2.0f);

  glColor4ub(0, 255, 0, 60);

  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(probe.min.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.max.y, probe.min.z);
  glVertex3f(probe.min.x, probe.max.y, probe.min.z);
  glEnd();

  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(probe.min.x, probe.min.y, probe.max.z);
  glVertex3f(probe.max.x, probe.min.y, probe.max.z);
  glVertex3f(probe.max.x, probe.max.y, probe.max.z);
  glVertex3f(probe.min.x, probe.max.y, probe.max.z);
  glEnd();

  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(probe.min.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.min.y, probe.max.z);
  glVertex3f(probe.min.x, probe.min.y, probe.max.z);
  glEnd();

  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(probe.min.x, probe.max.y, probe.min.z);
  glVertex3f(probe.max.x, probe.max.y, probe.min.z);
  glVertex3f(probe.max.x, probe.max.y, probe.max.z);
  glVertex3f(probe.min.x, probe.max.y, probe.max.z);
  glEnd();

  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(probe.min.x, probe.min.y, probe.min.z);
  glVertex3f(probe.min.x, probe.min.y, probe.max.z);
  glVertex3f(probe.min.x, probe.max.y, probe.max.z);
  glVertex3f(probe.min.x, probe.max.y, probe.min.z);
  glEnd();

  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(probe.max.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.min.y, probe.max.z);
  glVertex3f(probe.max.x, probe.max.y, probe.max.z);
  glVertex3f(probe.max.x, probe.max.y, probe.min.z);
  glEnd();

  glColor4ub(0, 255, 0, 255);
  glBegin(GL_LINE_LOOP);
  glVertex3f(probe.min.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.max.y, probe.min.z);
  glVertex3f(probe.min.x, probe.max.y, probe.min.z);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f(probe.min.x, probe.min.y, probe.max.z);
  glVertex3f(probe.max.x, probe.min.y, probe.max.z);
  glVertex3f(probe.max.x, probe.max.y, probe.max.z);
  glVertex3f(probe.min.x, probe.max.y, probe.max.z);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(probe.min.x, probe.min.y, probe.min.z);
  glVertex3f(probe.min.x, probe.min.y, probe.max.z);

  glVertex3f(probe.max.x, probe.min.y, probe.min.z);
  glVertex3f(probe.max.x, probe.min.y, probe.max.z);

  glVertex3f(probe.max.x, probe.max.y, probe.min.z);
  glVertex3f(probe.max.x, probe.max.y, probe.max.z);

  glVertex3f(probe.min.x, probe.max.y, probe.min.z);
  glVertex3f(probe.min.x, probe.max.y, probe.max.z);
  glEnd();

  ///// TEST
  /*
  glColor4ub(255, 0, 0, 255);

  OBB zone;
  zone.center = 0.5 * (probe.min + probe.max);
  zone.extent.set(0.5 * (probe.max.x - probe.min.x), 0.5 * (probe.max.y - probe.min.y), 0.5 * (probe.max.z -
  probe.min.z)); std::vector<size_t> pid; for (size_t i = 0; i < box.Particles.size(); ++i) {
    box.Particles[i].updateObb();
    if (zone.intersect(box.Particles[i].obb)) {
      pid.push_back(i);
    }
  }

  vec3r pt3;
  std::vector<double> vv(3);
  Mth::sobolSequence(-3, vv);  // Initialize the Sobol sequence
  //size_t count = 0;
  glBegin(GL_POINTS);
  for (size_t imc = 0; imc < 1000; ++imc) {
    Mth::sobolSequence(3, vv);
    pt3.set(probe.min.x + vv[0] * (probe.max.x - probe.min.x), probe.min.y + vv[1] * (probe.max.y - probe.min.y),
            probe.min.z + vv[2] * (probe.max.z - probe.min.z));

    //bool inSolid = false;
    for (size_t ii = 0; ii < pid.size(); ii++) {
      size_t i = pid[ii];
      vec3r ptTest = pt3 - box.Particles[i].pos;
      quat Qinv = box.Particles[i].Q.get_conjugated();
      ptTest = Qinv * ptTest;
      ptTest /= box.Particles[i].homothety;

      if (box.Particles[i].shape->inside(ptTest)) {
        //inSolid = true;
        glVertex3f(pt3.x, pt3.y, pt3.z);
        break;
      }
    }
    //if (inSolid) count++;
  }

  glEnd();
  */
}

void drawInteractionTypes() {
  if (mouse_mode != NOTHING && box.Particles.size() > 200) return;

  glDisable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  size_t i, j, isub, jsub;
  vec3r orig, dest;
  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {

      i = it->i;
      j = it->j;
      isub = it->isub;
      jsub = it->jsub;

      if (it->type == vvType) {  // ROUGE
        glColor3ub(255, 0, 0);
        orig = box.Particles[i].shape->vertex[isub];
        dest = box.Particles[j].shape->vertex[jsub];
      } else if (it->type == veType) {  // VERT
        glColor3ub(0, 255, 0);
        orig = box.Particles[i].shape->vertex[isub];
        dest = 0.5 * (box.Particles[j].shape->vertex[box.Particles[j].shape->edge[jsub].first] +
                      box.Particles[j].shape->vertex[box.Particles[j].shape->edge[jsub].second]);
      } else if (it->type == eeType) {  // BLEU
        glColor3ub(0, 0, 255);
        orig = 0.5 * (box.Particles[i].shape->vertex[box.Particles[i].shape->edge[isub].first] +
                      box.Particles[i].shape->vertex[box.Particles[i].shape->edge[isub].second]);
        dest = 0.5 * (box.Particles[j].shape->vertex[box.Particles[j].shape->edge[jsub].first] +
                      box.Particles[j].shape->vertex[box.Particles[j].shape->edge[jsub].second]);
      } else if (it->type == vfType) {  // ORANGE
        glColor3ub(255, 131, 0);
        orig = box.Particles[i].shape->vertex[isub];
        dest.reset();
        size_t nf = box.Particles[j].shape->face[jsub].size();
        for (size_t f = 0; f < nf; ++f) {
          dest += box.Particles[j].shape->vertex[box.Particles[j].shape->face[jsub][f]];
        }
        dest /= nf;
      } else
        continue;

      orig = box.Particles[i].Glob(orig);
      dest = box.Particles[j].Glob(dest);

      glLineWidth(2.0f);
      glBegin(GL_LINES);
      glVertex3f(orig.x, orig.y, orig.z);
      glVertex3f(dest.x, dest.y, dest.z);
      glEnd();
    }
  }
}

void drawInteractionFrames() {
  if (mouse_mode != NOTHING && box.Particles.size() > 200) return;

  glDisable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  size_t i, j;
  vec3r dnorm;
  glLineWidth(2.0f);
  glPointSize(5);
  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      i = it->i;
      j = it->j;

      if (it->dn > 0.0 && it->stick == nullptr) continue;

      if (it->type == vvType) {  // ROUGE
        glColor3ub(255, 0, 0);
      } else if (it->type == veType) {  // VERT
        glColor3ub(0, 255, 0);
      } else if (it->type == eeType) {  // BLEU
        glColor3ub(0, 0, 255);
      } else if (it->type == vfType) {  // ORANGE
        glColor3ub(255, 131, 0);
      } else
        continue;

      dnorm = it->n * 2.0 * (box.Particles[i].MinskowskiRadius() + box.Particles[j].MinskowskiRadius());

      glBegin(GL_LINES);
      glVertex3f(it->pos.x - dnorm.x, it->pos.y - dnorm.y, it->pos.z - dnorm.z);
      glVertex3f(it->pos.x + dnorm.x, it->pos.y + dnorm.y, it->pos.z + dnorm.z);
      glEnd();

      glBegin(GL_POINTS);
      glVertex3f(it->pos.x, it->pos.y, it->pos.z);
      glEnd();
    }
  }
}

void drawForces() {
  if (mouse_mode != NOTHING && box.Particles.size() > 200) return;

  glDisable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);

  vec3r force;

  double fnMax = 0.0;
  if (forceFactor == 0.0) {
    for (size_t k = 0; k < box.Interactions.size(); ++k) {
      std::set<Interaction>::iterator it = box.Interactions[k].begin();
      for (; it != box.Interactions[k].end(); ++it) {
        if (fabs(it->fn) > fnMax) fnMax = fabs(it->fn);
      }
    }
    if (fnMax > 0.0) {
      forceFactor = 0.005 / fnMax;  // 0.005 is the maximum length
    } else
      return;
  }

  glLineWidth(2.0f);
  glPointSize(5);
  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      // it is ok but stick is not yet saved in the conf files
      // if (it->dn > 0.0 && it->stick == nullptr) continue;

      if (it->fn > 0.0) {  // ROUGE
        glColor3ub(255, 0, 0);
      } else
        glColor3ub(0, 255, 0);  // VERT

      force = it->fn * it->n + it->ft;
      force *= forceFactor;

      glBegin(GL_LINES);
      glVertex3f(it->pos.x - force.x, it->pos.y - force.y, it->pos.z - force.z);
      glVertex3f(it->pos.x + force.x, it->pos.y + force.y, it->pos.z + force.z);
      glEnd();

      if (it->stick != nullptr) glColor3ub(0, 0, 0);  // NOIR

      glBegin(GL_POINTS);
      glVertex3f(it->pos.x, it->pos.y, it->pos.z);
      glEnd();
    }
  }
}

void drawVelocities() {
  /*
  // Scaling
  double velMax = 0.0;
  double velSqr;
  vec3r Vel;
  for (uint i = 0 ; i < box.Particles.size() ; ++i) {
          Vel = box.Cell.vh * box.Particles[i].pos + box.Cell.h *
  box.Particles[i].vel; // affine + fluctuation velSqr = norm2(Vel); if (velSqr
  > velMax) velMax = velSqr;
  }
  if (velMax == 0.0) return;
  velMax = sqrt(velMax);
  double scal = 5 * radiusMax / velMax;

  arrowSize = radiusMin / 2.0;

  glEnable (GL_LIGHTING);
  glColor4f (0.2f, 0.2f, 0.2f, 1.0f);
  GLColorRGBA color;
  for (uint i = 0 ; i < box.Particles.size() ; ++i) {
          vec3r pos = box.Cell.h * box.Particles[i].pos;
          if (!inSlice(pos)) continue;

          Vel = scal * (box.Cell.vh * box.Particles[i].pos + box.Cell.h *
  box.Particles[i].vel); if (1) { color = colorParticleVelocityMagnitude(i);
                  glColor4f (color.r, color.g, color.b, color.a);
          }
          drawArrow(pos, Vel);
  }
  */
}

void editSelection() {
  using namespace std;
  if (selectedParticle < 0) return;

  // https://en.wikipedia.org/wiki/Box-drawing_character
  cout << "\u256d\u2500" << endl;
  cout << "\u2502 Edit:" << endl;
  cout << "\u2502   1. position" << endl;
  cout << "\u2502   2. velocity" << endl;
  cout << "\u2502   3. angular position" << endl;
  cout << "\u2502   4. angular velocity" << endl;
  cout << "\u2570\u2500" << endl;
  int ans = 0;
  cout << "? \u25b7 ";
  cin >> ans;

  switch (ans) {
    case 1: {
      cout << "current pos = " << box.Particles[selectedParticle].pos << endl;
      cout << "    new pos = ";
      cin >> box.Particles[selectedParticle].pos;
    } break;
    case 2: {
      cout << "current vel = " << box.Particles[selectedParticle].vel << endl;
      cout << "    new vel = ";
      cin >> box.Particles[selectedParticle].vel;
    } break;
    case 3: {
      double angle = box.Particles[selectedParticle].Q.get_angle();
      vec3r axis = box.Particles[selectedParticle].Q.get_axis();
      cout << "current rotation quaternion = " << box.Particles[selectedParticle].Q << endl;
      cout << "current rotation angle = " << angle << endl;
      cout << "current rotation axis  = " << axis << endl;
      cout << "    new rotation angle = ";
      cin >> angle;
      cout << "    new rotation axis  = ";
      cin >> axis;
      box.Particles[selectedParticle].Q.set_axis_angle(axis, angle);
    } break;
    case 4: {
      cout << "current vrot = " << box.Particles[selectedParticle].vrot << endl;
      cout << "    new vrot = ";
      cin >> box.Particles[selectedParticle].vrot;
    } break;
    default:
      break;
  }
}

// This function is for making a screenshot (or a series of screenshots)
// The output format is PNG if linked with libpng, else it is TGA.
int screenshot(const char* filename) {
#ifdef PNG_H
  int i;

  int screenStats[4];
  glGetIntegerv(GL_VIEWPORT, screenStats);
  int x = screenStats[0];
  int y = screenStats[1];
  int width = screenStats[2];
  int height = screenStats[3];

  FILE* fp;
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep* rowp;
  GLubyte* glimage = NULL;

  fp = fopen(filename, "wb");
  if (!fp) return 1;

  glimage = (GLubyte*)malloc(width * height * sizeof(GLubyte) * 3);

  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)glimage);

  rowp = (png_bytep*)malloc(sizeof(png_bytep*) * height);
  if (!rowp) return 1;

  for (i = 0; i < height; i++) {
    rowp[i] = (png_bytep)&glimage[3 * ((height - i - 1) * width)];
  }

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) return 1;

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) return 1;

  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  png_write_image(png_ptr, rowp);
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);

  free(rowp);
  fflush(stdout);
  fclose(fp);
  return 0;

#else

  // http://forum.devmaster.net/t/rendering-a-single-frame-to-a-file-with-opengl/12469/2

  // we will store the image data here
  unsigned char* pixels;
  // the thingy we use to write files
  FILE* shot;
  // we get the width/height of the screen into this array
  int screenStats[4];

  // get the width/height of the window
  glGetIntegerv(GL_VIEWPORT, screenStats);

  // generate an array large enough to hold the pixel data
  // (width*height*bytesPerPixel)
  pixels = new unsigned char[screenStats[2] * screenStats[3] * 3];
  // read in the pixel data, TGA's pixels are BGR aligned
  glReadPixels(0, 0, screenStats[2], screenStats[3], GL_BGR, GL_UNSIGNED_BYTE, pixels);

  // open the file for writing. If unsucessful, return 1
  if ((shot = fopen(filename, "wb")) == NULL) return 1;

  // this is the tga header it must be in the beginning of
  // every (uncompressed) .tga
  unsigned char TGAheader[12] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  // the header that is used to get the dimensions of the .tga
  // header[1]*256+header[0] - width
  // header[3]*256+header[2] - height
  // header[4] - bits per pixel
  // header[5] - ?
  unsigned char header[6] = {((unsigned char)(screenStats[2] % 256)),
                             ((unsigned char)(screenStats[2] / 256)),
                             ((unsigned char)(screenStats[3] % 256)),
                             ((unsigned char)(screenStats[3] / 256)),
                             24,
                             0};

  // write out the TGA header
  fwrite(TGAheader, sizeof(unsigned char), 12, shot);
  // write out the header
  fwrite(header, sizeof(unsigned char), 6, shot);
  // write the pixels
  fwrite(pixels, sizeof(unsigned char), screenStats[2] * screenStats[3] * 3, shot);

  // close the file
  fclose(shot);
  // free the memory
  delete[] pixels;

  // return success
  return 0;
#endif
}

bool tryToReadConf(int num) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << "Read " << file_name << std::endl;
    box.clearMemory();
    box.loadConf(file_name);
    resetColors(colorMode, rescaleColorRange);
    totalNumberOfVertices = 0;
    for (size_t i = 0; i < box.Particles.size(); ++i) {
      totalNumberOfVertices += box.Particles[i].shape->vertex.size();
    }
    textZone.addLine("conf-file: %s (time = %f)", file_name, box.t);
    confNum = box.iconf;
    box.computeAABB();
    adjustClippingPlans();
  } else {
    std::cout << file_name << " does not exist" << std::endl;
    return false;
  }
  lastConfNumOK = num;
  return true;
}

void readTraj(const char* name) {
  if (!fileTool::fileExists(name)) {
    return;  // return silently (most of time it is not required)
  }
  std::cout << "Read trajectory file: " << name << std::endl;

  std::ifstream file(name);
  size_t nbReleases = 0;
  file >> nbReleases;
  std::cout << "Number of releases: " << nbReleases << std::endl;
  releases.clear();
  releases.resize(nbReleases);
  for (size_t r = 0; r < nbReleases; r++) {
    size_t nbFlights = 0;
    file >> nbFlights;
    for (size_t f = 0; f < nbFlights; f++) {
      freeFlight F;
      file >> F.ti >> F.duration >> F.pos >> F.vel >> F.Q >> F.vrot;
      releases[r].freeFlights.push_back(F);
    }
    releases[r].pos = releases[r].freeFlights[0].pos;
    releases[r].vel = releases[r].freeFlights[0].vel;
    releases[r].Q = releases[r].freeFlights[0].Q;
    releases[r].vrot = releases[r].freeFlights[0].vrot;
  }
  show_traj = 1;
}

void error(int error, const char* description) {
  std::cout << "GLFW error code: " << error << ", description: " << description << std::endl;
}

// =====================================================================
// IMGUI INPUT TOOLS FOR ROCKABLE
// =====================================================================

static void HelpPopup(const char* desc) {
  ImGui::TextDisabled("(?)");
  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
    ImGui::TextUnformatted(desc);
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
  }
}

void Input_vec3r(const char* label, vec3r& myVec) {
  vec3<float> tmp(myVec.x, myVec.y, myVec.z);
  if (ImGui::InputFloat3(label, &(tmp.x))) {
    myVec.set(tmp.x, tmp.y, tmp.z);
  }
}

void Input_quat(const char* label, quat& myQuat) {
  float tmp[4];
  tmp[0] = myQuat.v.x;
  tmp[1] = myQuat.v.y;
  tmp[2] = myQuat.v.z;
  tmp[3] = myQuat.s;
  if (ImGui::InputFloat4(label, tmp)) {
    myQuat.v.set(tmp[0], tmp[1], tmp[2]);
    myQuat.s = tmp[3];
    myQuat.normalize();
  }
}

void Input_Shape(Particle& P) {
  size_t item_current_idx = box.shapeId[P.shape->name];
  if (ImGui::BeginCombo("shape", P.shape->name.c_str(), 0)) {
    for (size_t n = 0; n < box.Shapes.size(); n++) {
      const bool is_selected = (item_current_idx == n);
      if (ImGui::Selectable(box.Shapes[n].name.c_str(), is_selected)) {
        item_current_idx = n;
        P.shape = &(box.Shapes[n]);  // Plug to the selected shape
      }
      if (is_selected) {
        ImGui::SetItemDefaultFocus();
      }
    }
    ImGui::EndCombo();
  }
}

void reset_ngroup() {
  size_t new_group_maxId = 0;
  for (size_t i = 0; i < box.Particles.size(); i++) {
    if (box.Particles[i].group > (int)new_group_maxId) new_group_maxId = box.Particles[i].group;
  }
  size_t new_ngroup = new_group_maxId + 1;
  if (new_ngroup != box.dataTable.get_ngroup()) {
    box.dataTable.set_ngroup(new_ngroup);
    box.properties.set_ngroup(new_ngroup);
  }
}

void Input_Properties(const char* propertyName) {
  for (size_t grp = 0; grp < box.properties.ngroup; grp++) {
    size_t idDensity = box.properties.get_id(propertyName);
    char label[256];
    snprintf(label, 256, "density (%ld)", grp);
    ImGui::InputDouble(label, &(box.properties.prop[idDensity][grp]));
  }
}

void Input_LawData(const char* parName) {
  std::string parNameStr(parName);
  size_t id = box.dataTable.data_id[parNameStr];
  for (size_t g1 = 0; g1 < box.dataTable.ngroup; ++g1) {
    for (size_t g2 = g1; g2 < box.dataTable.ngroup; ++g2) {
      if (1 /*box.dataTable.isDefined(id, g1, g2)*/) {
        char label[256];
        snprintf(label, 256, "%s (%ld - %ld)", parNameStr.c_str(), g1, g2);
        ImGui::InputDouble(label, &(box.dataTable.tables[id][g1][g2]));
      }
    }
  }
}

void ForceLaw_Combo() {
  std::vector<std::string> lawNames = {"Default", "Avalanche", "StickedLinks"};
  std::string lawName = box.optionNames["forceLaw"];
  size_t index = 0;
  for (size_t i = 0; i < lawNames.size(); i++) {
    if (lawName == lawNames[i]) {
      index = i;
      break;
    }
  }

  if (ImGui::BeginCombo("forceLaw", lawNames[index].c_str(), 0)) {
    for (size_t n = 0; n < lawNames.size(); n++) {
      const bool is_selected = (index == n);
      if (ImGui::Selectable(lawNames[n].c_str(), is_selected)) {
        index = n;
        // box.setForceLaw(lawNames[n]);
        std::stringstream ss;
        ss << "forceLaw " << lawNames[n];
        box.parser.parseString(ss.str().c_str());
      }
      if (is_selected) {
        ImGui::SetItemDefaultFocus();
      }
    }
    ImGui::EndCombo();
  }
}

void Integrator_Combo() {
  std::vector<std::string> Names = {"velocityVerlet", "Euler", "Beeman", "RungeKutta4"};
  std::string Name = box.optionNames["Integrator"];
  size_t index = 0;
  for (size_t i = 0; i < Names.size(); i++) {
    if (Name == Names[i]) {
      index = i;
      break;
    }
  }

  if (ImGui::BeginCombo("Integrator", Names[index].c_str(), 0)) {
    for (size_t n = 0; n < Names.size(); n++) {
      const bool is_selected = (index == n);
      if (ImGui::Selectable(Names[n].c_str(), is_selected)) {
        index = n;
        box.setIntegrator(Names[n]);
      }
      if (is_selected) {
        ImGui::SetItemDefaultFocus();
      }
    }
    ImGui::EndCombo();
  }
}

void UpdateNL_Combo() {
  std::vector<std::string> Names = {"bruteForce", "linkCells"};
  std::string Name = box.optionNames["UpdateNL"];
  size_t index = 0;
  for (size_t i = 0; i < Names.size(); i++) {
    if (Name == Names[i]) {
      index = i;
      break;
    }
  }

  if (ImGui::BeginCombo("UpdateNL", Names[index].c_str(), 0)) {
    for (size_t n = 0; n < Names.size(); n++) {
      const bool is_selected = (index == n);
      if (ImGui::Selectable(Names[n].c_str(), is_selected)) {
        index = n;
        box.setUpdateNL(Names[n]);
      }
      if (is_selected) {
        ImGui::SetItemDefaultFocus();
      }
    }
    ImGui::EndCombo();
  }
}

void AddOrRemoveInteractions_Combo() {
  std::vector<std::string> Names = {"bruteForce", "OBBtree"};
  std::string Name = box.optionNames["AddOrRemoveInteractions"];
  size_t index = 0;
  for (size_t i = 0; i < Names.size(); i++) {
    if (Name == Names[i]) {
      index = i;
      break;
    }
  }

  if (ImGui::BeginCombo("AddOrRemoveInteractions", Names[index].c_str(), 0)) {
    for (size_t n = 0; n < Names.size(); n++) {
      const bool is_selected = (index == n);
      if (ImGui::Selectable(Names[n].c_str(), is_selected)) {
        index = n;
        box.setAddOrRemoveInteractions(Names[n]);
      }
      if (is_selected) {
        ImGui::SetItemDefaultFocus();
      }
    }
    ImGui::EndCombo();
  }
}

// =====================================================================
// IMGUI WINDOWS
// =====================================================================

void window_conf_data() {
  ImGui::Begin("Configuration data...", &show_window_conf_data);
  ImGui::PushItemWidth(100);

  {
    float g[3];
    g[0] = box.gravity.x;
    g[1] = box.gravity.y;
    g[2] = box.gravity.z;
    ImGui::PushItemWidth(300);
    if (ImGui::InputFloat3("gravity", g)) {
      box.gravity.set(g[0], g[1], g[2]);
    }
    ImGui::PopItemWidth();
  }

  if (ImGui::CollapsingHeader("Options")) {
    ImGui::PushItemWidth(150);
    Integrator_Combo();
    ImGui::SameLine();
    HelpPopup("The integrator to be used");
    AddOrRemoveInteractions_Combo();
    ForceLaw_Combo();
    UpdateNL_Combo();
    ImGui::PopItemWidth();
  }

  if (ImGui::CollapsingHeader("Time flow")) {
    ImGui::InputDouble("current time", &(box.t), 0.0f, 0.0f, "%.8lf");
    ImGui::InputDouble("end time", &(box.tmax), 0.0f, 0.0f, "%.8lf");
    ImGui::InputDouble("dt", &(box.dt), 0.0f, 0.0f, "%.10lf");
    ImGui::InputDouble("Save period", &(box.interConf), 0.0f, 0.0f, "%.8lf");
  }

  if (ImGui::CollapsingHeader("Neighbor list")) {
    ImGui::InputDouble("Verlet update period", &(box.interVerlet));
    ImGui::InputDouble("R-shape alert distance", &(box.DVerlet));
    ImGui::SameLine();
    HelpPopup(
        "The size of each OBB is incremented by this value, "
        "and when two oversized-OBB intersect, they are add in the neighbor-list");
    ImGui::InputDouble("sub-element alert distance", &(box.dVerlet));
    ImGui::SameLine();
    HelpPopup(
        "When two R-shapes are neighbors, this distance is used for checking the proximity "
        "in-between the sub-elements (vertices, edges, and thick faces)");

    ImGui::CheckboxFlags("Dynamic update", &(box.dynamicUpdateNL), 1);
    ImGui::InputDouble("Distance for updating", &(box.dispUpdateNL));
    ImGui::InputDouble("Angle for updating", &(box.angleUpdateNL));
  }

  if (ImGui::CollapsingHeader("Numerical damping (0 = disabled)")) {
    ImGui::InputDouble("Cundall numerical Damping Coeff.", &(box.numericalDampingCoeff));
    ImGui::InputDouble("Velocity Barrier", &(box.velocityBarrier));
    ImGui::InputDouble("Angular Velocity Barrier", &(box.angularVelocityBarrier));
    ImGui::InputDouble("Velocity Barrier Exponent", &(box.velocityBarrierExponent));
    ImGui::InputDouble("Angular Velocity Barrier Exponent", &(box.angularVelocityBarrierExponent));
  }

  if (ImGui::CollapsingHeader("Properties")) {
    if (ImGui::Button("Reset No of groups")) {
      reset_ngroup();
    }
    ImGui::PushItemWidth(150);
    Input_Properties("density");
    ImGui::PopItemWidth();
  }

  if (ImGui::CollapsingHeader("Data tables")) {
    if (ImGui::Button("Reset No of groups")) {
      reset_ngroup();
    }
    ImGui::PushItemWidth(150);

    Input_LawData("knContact");
    Input_LawData("en2Contact");
    Input_LawData("ktContact");
    Input_LawData("muContact");
    Input_LawData("krContact");
    Input_LawData("murContact");

    if (box.optionNames["forceLaw"] == "StickedLinks") {
      Input_LawData("knInnerBond");
      Input_LawData("ktInnerBond");
      Input_LawData("fn0InnerBond");
      Input_LawData("ft0InnerBond");
      Input_LawData("powInnerBond");
      Input_LawData("en2InnerBond");

      Input_LawData("knOuterBond");
      Input_LawData("ktOuterBond");
      Input_LawData("krOuterBond");
      Input_LawData("fn0OuterBond");
      Input_LawData("ft0OuterBond");
      Input_LawData("mom0OuterBond");
      Input_LawData("powOuterBond");
      Input_LawData("en2OuterBond");
    }

    ImGui::PopItemWidth();
  }

  if (ImGui::Button("Close")) {
    show_window_conf_data = false;
  }
  ImGui::End();
}

void window_info() {
  ImGui::Begin("Information...", &show_window_info);
  ImGui::PushItemWidth(250);

  if (selectedParticle >= 0) {
    Particle& P = box.Particles[selectedParticle];
    ImGui::Text("Particle No %d", selectedParticle);
    Input_Shape(P);
    ImGui::InputDouble("homothety", &(P.homothety));
    Input_vec3r("inertia", P.inertia);
    ImGui::InputDouble("mass", &(P.mass));
    ImGui::InputInt("group", &(P.group));
    ImGui::InputInt("cluster", &(P.cluster));
    Input_vec3r("pos", P.pos);
    Input_vec3r("vel", P.vel);
    Input_vec3r("acc", P.acc);
    Input_quat("Q", P.Q);
    Input_vec3r("vrot", P.vrot);
    Input_vec3r("arot", P.arot);
  } else {
    ImGui::Text("cwd: %s", cwd);
    ImGui::Text("Total number of particles: %ld", box.Particles.size());
    ImGui::Text("Number of driven particles: %ld", box.nDriven);

    ImGui::Text("The current scene fits within:");
    ImGui::Text("    Xmin = %.08f   Xmax = %.08f", box.aabb.min.x, box.aabb.max.x);
    ImGui::Text("    Ymin = %.08f   Ymax = %.08f", box.aabb.min.y, box.aabb.max.y);
    ImGui::Text("    Zmin = %.08f   Zmax = %.08f", box.aabb.min.z, box.aabb.max.z);
  }

  if (ImGui::Button("Close")) show_window_info = false;

  ImGui::End();
}

void main_imgui_menu() {

  ImGui::Begin("Main menu");
  ImGui::PushItemWidth(100);

  if (ImGui::InputInt("iconf", &confNum)) {
    if (confNum < 0) confNum = 0;
    if (tryToReadConf(confNum) == false) {
      confNum = lastConfNumOK;
    }
  }

  if (ImGui::Button("Edit")) {
    show_window_conf_data = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("Save")) {
    box.saveConf("exportedConf");
  }
  ImGui::SameLine();
  if (ImGui::Button("Info")) {
    show_window_info = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("Quit")) {
    request_quit = true;
  }

  if (ImGui::CollapsingHeader("Show")) {
    ImGui::CheckboxFlags("Background", &show_background, 1);
    ImGui::CheckboxFlags("Scene frame", &show_globalFrame, 1);
    ImGui::CheckboxFlags("All particles", &show_particles, 1);
    ImGui::CheckboxFlags("Also driven particles", &show_driven, 1);
    ImGui::CheckboxFlags("Only faces (faster display)", &shapeWithoutThickness, 1);
    ImGui::CheckboxFlags("OBB", &show_obb, 1);
    ImGui::CheckboxFlags("OBB + Verlet distance", &enlarged_obb, 1);
  }

  if (ImGui::CollapsingHeader("Transparency")) {
    ImGui::SliderFloat("Particles Trans.", &alpha_particles, 0.0f, 1.0f);
    ImGui::SliderFloat("Driven particles Trans.", &alpha_fixparticles, 0.0f, 1.0f);
  }

  if (ImGui::CollapsingHeader("Camera")) {
    if (ImGui::Button("Fit")) {
      fitView();
      adjustClippingPlans();
    }
    ImGui::SameLine();
    if (ImGui::Button("Vertical gravity")) {
      vec3r d = center - eye;
      if (box.gravity.isnull())
        up.set(-d.x * d.y, d.x * d.x + d.z * d.z, -d.z * d.y);
      else {
        up = -cross(cross(d, box.gravity), d);
      }
      up.normalize();
    }

    if (ImGui::Button("+X +Y")) {
      double d = norm(center - eye);
      eye = center + d * vec3r::unit_z();
      up.set(0.0, 1.0, 0.0);
    }
    ImGui::SameLine();
    if (ImGui::Button("-X +Y")) {
      double d = norm(center - eye);
      eye = center - d * vec3r::unit_z();
      up.set(0.0, 1.0, 0.0);
    }

    if (ImGui::Button("+X -Z")) {
      double d = norm(center - eye);
      eye = center + d * vec3r::unit_y();
      up.set(0.0, 0.0, -1.0);
    }
    ImGui::SameLine();
    if (ImGui::Button("-X +Z")) {
      double d = norm(center - eye);
      eye = center - d * vec3r::unit_y();
      up.set(0.0, 0.0, 1.0);
    }

    if (ImGui::Button("+Z +Y")) {
      double d = norm(center - eye);
      eye = center - d * vec3r::unit_x();
      up.set(0.0, 1.0, 0.0);
    }
    ImGui::SameLine();
    if (ImGui::Button("-Z +Y")) {
      double d = norm(center - eye);
      eye = center + d * vec3r::unit_x();
      up.set(0.0, 1.0, 0.0);
    }
  }

  if (ImGui::CollapsingHeader("Interactions")) {
    ImGui::CheckboxFlags("Interaction frames", &show_interFrames, 1);
    ImGui::CheckboxFlags("Interaction types", &show_interTypes, 1);
    ImGui::CheckboxFlags("Forces", &show_forces, 1);
    ImGui::InputDouble("Force vector-length factor", &forceFactor);
  }

  if (ImGui::CollapsingHeader("Colors")) {
    ImGui::PushItemWidth(150);
    ImGui::ColorEdit3("Background bottom", (float*)&clear_bottom_color);
    ImGui::ColorEdit3("Background top", (float*)&clear_top_color);
    ImGui::PopItemWidth();
  }

  ImGui::End();
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {
	RockableProfiler::ProfilerManager prof;
	
  box.setInteractive(true);

  std::string confFileName;
  std::string trajFileName;

  try {
    TCLAP::CmdLine cmd("Visualization of Rockable simulations", ' ', "0.3");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the conf-file", false, "conf0", "conf-file");
    TCLAP::ValueArg<std::string> trajFileNameArg("t", "traj", "Name of a trajectory file", false, "traj.txt", "string");

    cmd.add(nameArg);
    cmd.add(trajFileNameArg);

    cmd.parse(argc, argv);

    confFileName = nameArg.getValue();
    trajFileName = trajFileNameArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

  box.clearMemory();
  box.loadConf(confFileName.c_str());
  textZone.addLine("conf-file: %s (time = %f)", confFileName.c_str(), box.t);
  box.computeAABB();

  getcwd(cwd, sizeof(cwd));

  if (fileTool::fileExists(trajFileName.c_str())) {
    readTraj(trajFileName.c_str());
  }

  if (fileTool::fileExists("probe.txt")) {
    std::ifstream f("probe.txt");
    f >> probe.min >> probe.max >> probe_MCnsteps;
    std::cout << "Probe: min = " << probe.min << "\n       max = " << probe.max << "\n";
  }

  totalNumberOfVertices = 0;
  for (size_t i = 0; i < box.Particles.size(); ++i) {
    totalNumberOfVertices += box.Particles[i].shape->vertex.size();
  }

  box.System.read();
  confNum = box.iconf;

  if (box.Particles.empty()) {
    std::cerr << "No particles! Goodbye." << std::endl;
    return 1;
  }

  // ==== Setup GLFW
  glfwSetErrorCallback(error);
  if (!glfwInit()) return -1;

  glfwWindowHint(GLFW_SAMPLES, 4);

  GLFWwindow* window = glfwCreateWindow(width, height, "Rockable visualiser (GLFW + Dear ImGui)", NULL, NULL);
  if (!window) {
    glfwTerminate();
    return -1;
  }

  glfwSetKeyCallback(window, keyboard);
  glfwSetMouseButtonCallback(window, mouse);
  glfwSetCursorPosCallback(window, motion);
  glfwSetFramebufferSizeCallback(window, reshape);

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  // ==== Setup Dear ImGui
  // context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  (void)io;
  // style
  ImGui::StyleColorsDark();
  // Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL2_Init();

  // ==== Init the visualizer
  center.set(0.0, 0.0, 0.0);  // where we look at
  eye.set(0.0, 0.0, 1.0);     // from where we look
  up.set(0.0, 1.0, 0.0);      // direction (normalized)

  mouse_mode = NOTHING;
  view_angle = 45.0;
  znear = 0.01;
  zfar = 10.0;

  glText::init();

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

  // print things
  std::cout << "GL VERSION: " << glGetString(GL_VERSION) << '\n';

  // ==== Enter event processing cycle
  adjustClippingPlans();
  fitView();

  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if (ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
      ImGui_window_focused = true;
    else
      ImGui_window_focused = false;

    main_imgui_menu();
    if (show_window_conf_data) window_conf_data();
    if (show_window_info) window_info();

    if (show_demo_window) ImGui::ShowDemoWindow(&show_demo_window);

    // Rendering
    ImGui::Render();

    if (request_quit == true) break;

    display();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

    glfwMakeContextCurrent(window);
    glfwSwapBuffers(window);
  }

  // Dear ImGui cleanup
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  // GLFW cleanup
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
