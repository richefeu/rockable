// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#include "see.hpp"

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
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[+]    load next configuration file");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[-]    load previous configuration file");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[=]    fit the view");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[a][A] decrease/increase alpha (transparence) of driven bodies");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[b]    switch ON/OFF the background color");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[c]    run 5000 steps of computation (used for debugging)");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[e][E] decrease/increase alpha (transparence) of free bodies");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[g]    open another file");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[k]    print this help");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[l]    switch ON/OFF the links (normal vector at contact)");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[m]    switch ON/OFF the links colored by type");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[n]    switch ON/OFF the body displays");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[o]    switch ON/OFF the OBB displays");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_,
                "[O]    switch ON/OFF the enlargement of OBBs by the Verlet distance");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[w]    set the view so that gravity appears vertical");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[x]    print the space limits of the current scene");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[p]    edit selected body");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[q]    quit");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[y]    make the display faster (and less nice)");
#ifdef PNG_SCREENSHOTS
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[z]    make a screenshot (oneshot.png)");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[Z]    make a series of screenshots (shotX.png)");
#else
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[z]    make a screenshot (oneshot.tga)");
  glText::print(GLUT_BITMAP_8_BY_13, 15, _nextLine_, "[Z]    make a series of screenshots (shotX.tga)");
#endif

#undef _nextLine_
  switch2D::back();
}

void keyboardSpecial(int Key, int /*x*/, int /*y*/) {
  switch (Key) {
    case GLUT_KEY_UP:
      textZone.increase_nbLine();
      break;

    case GLUT_KEY_DOWN:
      textZone.decrease_nbLine();
      break;
  };

  glutPostRedisplay();
}

void keyboard(unsigned char Key, int x, int y) {
  switch (Key) {
    case 27:  // ESCAPE
      selectedParticle = -1;
      break;

    case ' ':
      textZone.reset();
      break;

    case '*':
      selection(x, y);
      break;

    case 'a':
      if (alpha_fixparticles > 0.1) alpha_fixparticles -= 0.05;
      break;
    case 'A':
      if (alpha_fixparticles <= 0.95) alpha_fixparticles += 0.05;
      break;

    case 'b':
      show_background = 1 - show_background;
      break;

    case 'c': {  // compute a few steps (for debugging)
      box.UpdateNL();
      for (int i = 0; i < 5000; i++) box.velocityVerletStep();
      textZone.addLine("5000 time steps have been done.");
    } break;

    case 'e':
      if (alpha_particles > 0.1) alpha_particles -= 0.05;
      break;
    case 'E':
      if (alpha_particles <= 0.95) alpha_particles += 0.05;
      break;

    case 'f':
      show_forces = 1 - show_forces;
      break;

    case 'h': {
      textZone.addLine("");
      textZone.addLine("MOUSE:");
      textZone.addLine("        Left button + move = ROTATE");
      textZone.addLine("SHIFT + Left button + move = PAN");
      textZone.addLine("CTRL  + Left button + move = ZOOM (or middle button)");
      textZone.addLine("ALT   + Left button        = SELECT (or SHIFT + middle button)");
      textZone.addLine("The KEYBOARD shortcuts -> key [k]");
    } break;

    case 'i': {
      if (selectedParticle == -1) {  // no particle selected
        textZone.addLine("");
        textZone.addLine("Nb Particles: %d,  t: %g  dt: %g,  tmax: %g", box.Particles.size(), box.t, box.dt, box.tmax);
        textZone.addLine("interVerlet: %g,  DVerlet: %g,  dVerlet: %g", box.interVerlet, box.DVerlet, box.dVerlet);

        double dtc;
        box.getCriticalTimeStep(dtc);
        textZone.addLine("dtc/dt: %g", dtc / box.dt);
      } else {
        Particle* P = &(box.Particles[selectedParticle]);
        textZone.addLine("");
        textZone.addLine("Particle#%d,  group#%d,  cluster#%d,  time: %g ", selectedParticle, P->group, P->cluster,
                         box.t);
        textZone.addLine("pos: %g %g %g", P->pos.x, P->pos.y, P->pos.z);
        textZone.addLine("vel: %g %g %g,  acc: %g %g %g", P->vel.x, P->vel.y, P->vel.z, P->acc.x, P->acc.y, P->acc.z);
        textZone.addLine("Q: %g %g %g %g", P->Q.s, P->Q.v.x, P->Q.v.y, P->Q.v.z);
        textZone.addLine("vrot: %g %g %g,  arot: %g %g %g", P->vrot.x, P->vrot.y, P->vrot.z, P->arot.x, P->arot.y,
                         P->arot.z);
        textZone.addLine("shape: %s,  homothety: %g", P->shape->name.c_str(), P->homothety);
        textZone.addLine("inertia: %g %g %g,  mass: %g", P->inertia.x, P->inertia.y, P->inertia.z, P->mass);
      }
    } break;

    case 'k':
      show_keybinds = 1 - show_keybinds;
      break;

    case 'g': {
      std::cout << ">>>>>>>>>> CONF NUMBER: ";
      int Num;
      std::cin >> Num;
      tryToReadConf(Num);
    } break;

    case 'l':
      show_interFrames = 1 - show_interFrames;
      break;

    case 'm':
      show_interTypes = 1 - show_interTypes;
      break;

    case 'n':
      show_particles = 1 - show_particles;
      break;

    case 'o':
      show_obb = 1 - show_obb;
      break;
    case 'O':
      enlarged_obb = 1 - enlarged_obb;
      break;

    case 'p': {
      editSelection();
    } break;

    case 'q':
      exit(0);
      break;

    case 's':
      show_slice = 1 - show_slice;
      break;

    case 't':
      if (forceTubeFactor > 0.2) forceTubeFactor -= 0.05;
      break;
    case 'T':
      if (forceTubeFactor <= 0.9) forceTubeFactor += 0.05;
      break;

    case 'v':
      show_velocities = 1 - show_velocities;
      break;

    case 'w': {
      vec3r d = center - eye;
      up.set(-d.x * d.y, d.x * d.x + d.z * d.z, -d.z * d.y);
      up.normalize();
    } break;

    case 'x': {
      msg::bestPrecision(std::cout);
      std::cout << "The current scene fits within these limits:" << std::endl;
      std::cout << "Xmin = " << box.aabb.min.x << std::endl;
      std::cout << "Xmax = " << box.aabb.max.x << std::endl;
      std::cout << "Ymin = " << box.aabb.min.y << std::endl;
      std::cout << "Ymax = " << box.aabb.max.y << std::endl;
      std::cout << "Zmin = " << box.aabb.min.z << std::endl;
      std::cout << "Zmax = " << box.aabb.max.z << std::endl;
      msg::normalPrecision(std::cout);

      textZone.addLine("");
      textZone.addLine(
          "The current scene fits within these limits (see console for high "
          "precision):");
      textZone.addLine(" Xmin = %f, Xmax = %f", box.aabb.min.x, box.aabb.max.x);
      textZone.addLine(" Ymin = %f, Ymax = %f", box.aabb.min.y, box.aabb.max.y);
      textZone.addLine(" Zmin = %f, Zmax = %f", box.aabb.min.z, box.aabb.max.z);
    } break;

    case 'y':
      complexMode = 1 - complexMode;
      break;

    case 'z': {
#ifdef PNG_H
      screenshot("oneshot.png");
#else
      screenshot("oneshot.tga");
#endif
    } break;
    case 'Z': {
      // be carreful there's no way to stop this loop
      // if the process if too long
      while (tryToReadConf(confNum + 1)) {
        char name[256];
#ifdef PNG_H
        sprintf(name, "shot%d.png", confNum);
#else
        sprintf(name, "shot%d.tga", confNum);
#endif
        display();
        screenshot(name);
      }
    } break;

    case '-':
      if (confNum > 0) tryToReadConf(confNum - 1);
      break;

    case '+':
      tryToReadConf(confNum + 1);
      break;

    case '=': {
      fitView();
      adjustClippingPlans();
    } break;
  };

  glutPostRedisplay();
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
    for (size_t i = 0; i < box.Particles.size(); ++i) {
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

void mouse(int button, int state, int x, int y) {
  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
          mouse_mode = PAN;
        else if (glutGetModifiers() == GLUT_ACTIVE_ALT)
          selection(x, y);
        else if (glutGetModifiers() == GLUT_ACTIVE_CTRL)
          mouse_mode = ZOOM;
        else
          mouse_mode = ROTATION;
        break;
      case GLUT_MIDDLE_BUTTON:
        mouse_mode = ZOOM;
        break;
    }
  }
}

void motion(int x, int y) {
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

  display();
}

void display() {
  glutTools::clearBackground(show_background);
  adjustClippingPlans();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  // if (show_velocities) drawVelocities();
  if (show_traj) drawTrajectories();
  if (show_forces) drawForces();
  if (show_interFrames) drawInteractionFrames();
  if (show_interTypes) drawInteractionTypes();
  if (show_obb) drawOBBs();
  if (show_particles) drawParticles();

  if (show_keybinds) showKeybinds();
  textZone.draw();
  glFlush();
  glutSwapBuffers();
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

void reshape(int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, width, height);

  adjustClippingPlans();
  glutPostRedisplay();
}

// Draw the shape of the sphero-polyhedron in its own framework
void drawShape(Shape* s, double homothety) {
  double R = homothety * s->radius;
  int nbLevelSphere = 2;
  if (complexityNumber > 10000) {
    nbLevelSphere = 1;
  }

  if (complexMode == 0) {
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
      glutShape::drawTube(orig, arrow, 2.0 * R);
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
  if (mouse_mode != NOTHING && complexityNumber > 5000) {
    drawOBBs();
    return;
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  for (size_t i = box.nDriven; i < box.Particles.size(); ++i) {
    if (selectedParticle >= 0 && i == (size_t)selectedParticle) {
      glColor4ub(234, 255, 0, (int)floor(alpha_particles * 255));
    } else
      glColor4ub(178, 34, 34, (int)floor(alpha_particles * 255));

    vec3r pos = box.Particles[i].pos;

    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);
    quat2GLMatrix<GLfloat>(box.Particles[i].Q, Rot_Matrix);
    glMultMatrixf(Rot_Matrix);
    drawShape(box.Particles[i].shape, box.Particles[i].homothety);
    glPopMatrix();
  }

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

void drawOBBs() {
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glLineWidth(1.0f);
  /*
  for (size_t i = 0 ; i < box.Particles.size() ; ++i) {
          glColor4ub(255, 0, 0, 255); // RED

          double h = box.Particles[i].homothety;
          vec3r pos = box.Particles[i].pos;
          OBB obb = box.Particles[i].shape->obb;
          pos += obb.center * h;
          if (enlarged_obb) obb.enlarge(0.5 * box.DVerlet);
          vec3r ext = obb.extent;
          vec3r e0 = h * obb.e[0] * ext.x;
          vec3r e1 = h * obb.e[1] * ext.y;
          vec3r e2 = h * obb.e[2] * ext.z;

          //if (!inSlice(pos)) continue;
          glPushMatrix ();
          glTranslatef (pos.x, pos.y, pos.z);
          quat2GLMatrix<GLfloat> (box.Particles[i].Q, Rot_Matrix);
          glMultMatrixf (Rot_Matrix);

          vec3r corner;
          glBegin (GL_LINE_LOOP);
          corner = e0 + e1 + e2; glVertex3f (corner.x, corner.y, corner.z);
          corner -= 2.0 * e1; glVertex3f (corner.x, corner.y, corner.z);
          corner -= 2.0 * e2; glVertex3f (corner.x, corner.y, corner.z);
          corner += 2.0 * e1; glVertex3f (corner.x, corner.y, corner.z);
          corner += 2.0 * e2; glVertex3f (corner.x, corner.y, corner.z);
          glEnd();
          glBegin (GL_LINE_LOOP);
          corner -= 2.0 * e0; glVertex3f (corner.x, corner.y, corner.z);
          corner -= 2.0 * e1; glVertex3f (corner.x, corner.y, corner.z);
          corner -= 2.0 * e2; glVertex3f (corner.x, corner.y, corner.z);
          corner += 2.0 * e1; glVertex3f (corner.x, corner.y, corner.z);
          corner += 2.0 * e2; glVertex3f (corner.x, corner.y, corner.z);
          glEnd();

          glBegin (GL_LINES);
          glVertex3f (corner.x, corner.y, corner.z);
          corner += 2.0 * e0; glVertex3f (corner.x, corner.y, corner.z);

          corner -= 2.0 * e1; glVertex3f (corner.x, corner.y, corner.z);
          corner -= 2.0 * e0; glVertex3f (corner.x, corner.y, corner.z);

          corner -= 2.0 * e2; glVertex3f (corner.x, corner.y, corner.z);
          corner += 2.0 * e0; glVertex3f (corner.x, corner.y, corner.z);

          corner += 2.0 * e1; glVertex3f (corner.x, corner.y, corner.z);
          corner -= 2.0 * e0; glVertex3f (corner.x, corner.y, corner.z);
          glEnd();

          glPopMatrix ();
  }
*/
  /// ****

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
    vec3r corner;
    glBegin(GL_LINE_LOOP);
    corner = obbi.center + obbi.extent[0] * obbi.e[0] + obbi.extent[1] * obbi.e[1] + obbi.extent[2] * obbi.e[2];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[1] * obbi.e[1];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[2] * obbi.e[2];
    glVertex3f(corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[1] * obbi.e[1];
    glVertex3f(corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[2] * obbi.e[2];
    glVertex3f(corner.x, corner.y, corner.z);
    glEnd();
    glBegin(GL_LINE_LOOP);
    corner -= 2.0 * obbi.extent[0] * obbi.e[0];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[1] * obbi.e[1];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[2] * obbi.e[2];
    glVertex3f(corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[1] * obbi.e[1];
    glVertex3f(corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[2] * obbi.e[2];
    glVertex3f(corner.x, corner.y, corner.z);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[0] * obbi.e[0];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[1] * obbi.e[1];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[0] * obbi.e[0];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[2] * obbi.e[2];
    glVertex3f(corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[0] * obbi.e[0];
    glVertex3f(corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[1] * obbi.e[1];
    glVertex3f(corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[0] * obbi.e[0];
    glVertex3f(corner.x, corner.y, corner.z);
    glEnd();

    /*
    for (size_t j = i + 1 ; j < box.Particles.size() ; j++) {

            obbj = box.Particles[j].shape->obb;
            obbj.rotate(box.Particles[j].Q);
            obbj.extent *= box.Particles[j].homothety;
            obbj.center *= box.Particles[j].homothety;
            obbj.center += box.Particles[j].pos;
            obbj.enlarge(0.5 * box.DVerlet);

            // check intersection
            if (obbi.intersect(obbj)) {
                    glBegin (GL_LINES);
                    glColor4ub(0, 0, 255, 255); // BLUE
                    glVertex3f (obbi.center.x, obbi.center.y, obbi.center.z);
                    glVertex3f (obbj.center.x, obbj.center.y, obbj.center.z);
                    //glColor4ub(0, 255, 0, 255); // GREEN
                    //glVertex3f (box.Particles[i].pos.x,
    box.Particles[i].pos.y, box.Particles[i].pos.z);
                    //glVertex3f (box.Particles[j].pos.x,
    box.Particles[j].pos.y, box.Particles[j].pos.z);

                    glEnd();
            }
    }
    */
  }
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
  // size_t i, j;
  vec3r force;
  double forceFactor = 0.0;  // tmp #########

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

      // i = it->i;
      // j = it->j;

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
  sprintf(file_name, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << "Read " << file_name << std::endl;
    box.clearMemory();
    box.loadConf(file_name);
    complexityNumber = 0;
    for (size_t i = 0; i < box.Particles.size(); ++i) {
      complexityNumber += box.Particles[i].shape->vertex.size();
    }
    textZone.addLine("conf-file: %s (time = %f, complexityNumber = %ld)", file_name, box.t, complexityNumber);
    confNum = box.iconf;
    box.computeAABB();
    adjustClippingPlans();
  } else {
    std::cout << file_name << " does not exist" << std::endl;
    return false;
  }
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

void menu(int num) {
  switch (num) {

    case 0:
      exit(0);
      break;

    case 100:
      // colorParticle = colorParticleNone;
      break;
    case 101: {
      // colorParticle = colorParticleVelocityMagnitude;
      // ParticleColorTable.setMinMax(0.0, box.VelMax);
      // ParticleColorTable.Rebuild();
    } break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  // int popupmenu, submenu1, submenu2, submenu3;

  /*
  submenu1 = glutCreateMenu(menu); // Particle Colors
  glutAddMenuEntry("None", 100);
  glutAddMenuEntry("Velocity Magnitude", 101);
  glutAddMenuEntry("Sum of Normal Contact Forces", 102);

  submenu2 = glutCreateMenu(menu); // Force Colors
  glutAddMenuEntry("None", 200);
  glutAddMenuEntry("Magnitude", 201);

  submenu3 = glutCreateMenu(menu); // Velocity Colors
  glutAddMenuEntry("None", 300);
  glutAddMenuEntry("Magnitude", 301);
*/
  // popupmenu =
  glutCreateMenu(menu);  // Main menu
  // glutAddSubMenu("Particle Colors", submenu1);
  // glutAddSubMenu("Force Colors", submenu2);
  // glutAddSubMenu("Velocity Colors", submenu2);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {
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

  if (fileTool::fileExists(trajFileName.c_str())) {
    readTraj(trajFileName.c_str());
  }

  complexityNumber = 0;
  for (size_t i = 0; i < box.Particles.size(); ++i) {
    complexityNumber += box.Particles[i].shape->vertex.size();
  }

  box.System.read();
  confNum = box.iconf;

  if (box.Particles.empty()) {
    std::cerr << "No particles! Goodbye." << std::endl;
    return 1;
  }

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("Rockable VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(keyboardSpecial);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

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
  adjustClippingPlans();
  fitView();
  glutMainLoop();
  return 0;
}