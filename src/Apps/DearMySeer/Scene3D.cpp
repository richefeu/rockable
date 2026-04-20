#include "Scene3D.hpp"
#include "glTools.hpp"

Scene3D::Scene3D()
    : box(),
      rotationAngle(0.0f),
      useWireframe(false),
      windowWidth(1280),
      windowHeight(720),
      lastMouseX(0),
      lastMouseY(0),
      mouseDragging(false),
      showAxis(true),                // Afficher le repère 3D par défaut
      axisSize(2.0f),                // Taille du repère
      xAxisColor(1.0f, 0.0f, 0.0f),  // Rouge pour X
      yAxisColor(0.0f, 1.0f, 0.0f),  // Vert pour Y
      zAxisColor(0.0f, 0.0f, 1.0f),  // Bleu pour Z
      arrowFaces(16),                // 16 faces pour les flèches (lisse)
      useLighting(true),             // Activer l'éclairage par défaut
      arrowLengthRatio(0.3f),        // Rapport de longueur de flèche (30% de la taille de l'axe)
      arrowRadiusRatio(0.5f)         // Rapport de rayon de flèche (50% de la longueur)
{
  // Les paramètres de caméra sont déjà initialisés avec les valeurs par défaut
  // définies dans l'en-tête (cameraPosition, cameraTarget, etc.)

  // Calculer la distance de la caméra à la cible
  cameraDistance = glm::distance(cameraPosition, cameraTarget);

  // Initialiser les matrices avec GLM
  updateProjection();
  updateCameraView();
}

bool Scene3D::init() {
  // Initialisation des paramètres OpenGL fixes

  // Désactive le culling (suppression des faces arrière), ce qui signifie que toutes les faces des polygones seront
  // rendues, même celles qui sont orientées vers l'arrière.
  glDisable(GL_CULL_FACE);

  // Définit le mode de rendu des polygones : ici, les polygones sont remplis (GL_FILL) pour les faces avant et arrière.
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Active l'éclairage, ce qui permet d'utiliser des sources de lumière pour calculer les couleurs des objets.
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

  // Active la gestion des matériaux colorés, ce qui permet de définir la couleur des objets directement avec glColor.
  glEnable(GL_COLOR_MATERIAL);

  // Active le test de profondeur (Z-buffer), ce qui permet de gérer correctement l'occlusion des objets en fonction de
  // leur distance par rapport à la caméra.
  glEnable(GL_DEPTH_TEST);

  // Définit la valeur de profondeur par défaut utilisée lors du nettoyage du buffer de profondeur
  // (glClear(GL_DEPTH_BUFFER_BIT)).
  glClearDepth(1.0f);

  // Définit le modèle d'ombrage utilisé pour les polygones : ici, l'ombrage est lissé (GL_SMOOTH), ce qui permet des
  // transitions douces entre les couleurs des sommets.
  glShadeModel(GL_SMOOTH);

  // Active l'antialiasing pour les points, ce qui permet de lisser les bords des points affichés.
  glEnable(GL_POINT_SMOOTH);

  // Définit la qualité de l'antialiasing pour les points : ici, la meilleure qualité est demandée (GL_NICEST).
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  // Active le mélange de couleurs (transparence),
  // ce qui permet de gérer les effets de transparence et de superposition des objets.
  glEnable(GL_BLEND);

  // Définit l'équation de mélange pour les couleurs : ici, l'addition classique est utilisée.
  glBlendEquation(GL_FUNC_ADD);

  // Définit la fonction de mélange pour les couleurs : ici, la transparence est calculée en utilisant la valeur alpha
  // de la source (GL_SRC_ALPHA) et l'inverse de la valeur alpha de la source (GL_ONE_MINUS_SRC_ALPHA).
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glText::init();

  return true;
}

void Scene3D::updateProjection() {
  wh_ratio = static_cast<float>(windowWidth) / static_cast<float>(windowHeight);

  float zf = glm::length(cameraPosition - cameraTarget);
  glm::vec3 aabbmax(box.aabb.max.x, box.aabb.max.y, box.aabb.max.z);
  glm::vec3 aabbmin(box.aabb.min.x, box.aabb.min.y, box.aabb.min.z);
  glm::vec3 mx = aabbmax - aabbmin;
  float diag_length = glm::length(mx);
  float min_length = std::min({mx.x, mx.y, mx.z});
  // float max_length = std::max({mx.x, mx.y, mx.z});
  cameraDistanceMin = min_length * 0.0001f;
  cameraDistanceMax = diag_length * 2.5f;

  znear = zf - diag_length;
  double close_dst = 0.1 * zf;
  if (znear < close_dst) znear = close_dst;
  zfar = zf + diag_length;

  projectionMatrix = glm::perspective(glm::radians(cameraFov), wh_ratio, znear, zfar);
}

void Scene3D::updateCameraView() {
  viewMatrix = glm::lookAt(cameraPosition,  // Position de la caméra
                           cameraTarget,    // Point cible
                           cameraUp         // Vecteur haut
  );
}

void Scene3D::fitView() {
  glm::vec3 direction = glm::normalize(cameraPosition - cameraTarget);
  glm::vec3 aabbmax(box.aabb.max.x, box.aabb.max.y, box.aabb.max.z);
  glm::vec3 aabbmin(box.aabb.min.x, box.aabb.min.y, box.aabb.min.z);
  float maxDimension = std::max({aabbmax.x - aabbmin.x, aabbmax.y - aabbmin.y, aabbmax.z - aabbmin.z});
  cameraTarget = 0.5f * (aabbmax + aabbmin);
  cameraDistance = 0.5f * maxDimension / (float)(tan(glm::radians(0.5f * cameraFov)));
  cameraPosition = cameraTarget + cameraDistance * direction;
}

void Scene3D::selection(int x, int y) {

  selectedParticle = -1;

  // --- 1. Coordonnées normalisées (NDC) ---
  float ndcX = (2.0f * x) / windowWidth - 1.0f;
  float ndcY = 1.0f - (2.0f * y) / windowHeight;

  // Clip space
  glm::vec4 rayClip(ndcX, ndcY, -1.0f, 1.0f);

  // --- 2. Eye space ---
  glm::mat4 invProj = glm::inverse(projectionMatrix);
  glm::vec4 rayEye = invProj * rayClip;
  rayEye = glm::vec4(rayEye.x, rayEye.y, -1.0f, 0.0f);

  // --- 3. World space ---
  glm::mat4 invView = glm::inverse(viewMatrix);
  glm::vec3 rayDir = glm::normalize(glm::vec3(invView * rayEye));

  glm::vec3 rayOrigin = cameraPosition;

  // --- 4. Intersection avec particules ---
  float closestT = std::numeric_limits<float>::max();

  if (show_particles) {

    size_t i0 = 0;
    if (show_driven == 0) i0 = box.nDriven;

    for (size_t i = i0; i < box.Particles.size(); ++i) {

      vec3r pos = box.Particles[i].pos;

      // Approximation : sphère FAKE !!!!!!!!!!!! il faut viser le centre de la particule
      // Il faudrait faire l'intersection avec le OBB
      // mais pour le moment ça dépanne
      float radius = box.Particles[i].shape->radius * box.Particles[i].homothety;

      glm::vec3 center(pos.x, pos.y, pos.z);

      // --- Ray / Sphere intersection ---
      glm::vec3 oc = rayOrigin - center;

      float a = glm::dot(rayDir, rayDir);
      float b = 2.0f * glm::dot(oc, rayDir);
      float c = glm::dot(oc, oc) - radius * radius;

      float discriminant = b * b - 4 * a * c;

      if (discriminant < 0) continue;

      float t = (-b - sqrt(discriminant)) / (2.0f * a);

      if (t > 0.0f && t < closestT) {
        closestT = t;
        selectedParticle = (int)i;
      }
    }
  }

  std::cout << "PARTICLE SELECTED = " << selectedParticle << std::endl;
}

void Scene3D::handleMouseEvent(const SDL_Event& event, int, int) {

  // Vérifier si ImGui est en train de capturer les événements souris
  // Si oui, ne pas traiter les événements pour la caméra
  ImGuiIO& io = ImGui::GetIO();
  if (io.WantCaptureMouse) {
    return;
  }

  // Vérifier l'état des touches modificatrices
  SDL_Keymod mods = SDL_GetModState();
  bool shiftPressed = (mods & KMOD_SHIFT) != 0;
  bool ctrlPressed = (mods & KMOD_CTRL) != 0;

  switch (event.type) {
    case SDL_MOUSEBUTTONDOWN:
      if (event.button.button == SDL_BUTTON_LEFT && !mods) {
        // Bouton gauche seul : préparer pour la rotation
        mouseDragging = true;
        lastMouseX = event.button.x;
        lastMouseY = event.button.y;
      } else if (event.button.button == SDL_BUTTON_RIGHT) {
        // Bouton droit : préparer pour le zoom
        mouseDragging = true;
        lastMouseX = event.button.x;
        lastMouseY = event.button.y;
      } else if (event.button.button == SDL_BUTTON_LEFT && shiftPressed) {
        // Shift + bouton gauche : préparer pour le pan
        mouseDragging = true;
        lastMouseX = event.button.x;
        lastMouseY = event.button.y;
      } else if (event.button.button == SDL_BUTTON_LEFT && ctrlPressed) {
        // Ctrl + bouton gauche : selection
        selection(event.button.x, event.button.y);
      }
      break;

    case SDL_MOUSEBUTTONUP:
      if (event.button.button == SDL_BUTTON_LEFT || event.button.button == SDL_BUTTON_RIGHT) {
        mouseDragging = false;
      }
      break;

    case SDL_MOUSEMOTION:
      if (mouseDragging) {
        // Calculer le déplacement de la souris
        int dx = event.motion.x - lastMouseX;
        int dy = event.motion.y - lastMouseY;

        // Vérifier à nouveau l'état des modificateurs pour le mouvement
        mods = SDL_GetModState();
        shiftPressed = (mods & KMOD_SHIFT) != 0;

        if (!shiftPressed && (event.motion.state & SDL_BUTTON_LMASK)) {
          // Rotation avec bouton gauche (sans Shift)
          float sensitivity = 0.5f;

          // Rotation horizontale (autour de l'axe Y)
          glm::mat4 rotationY = glm::rotate(glm::mat4(1.0f), glm::radians(-dx * sensitivity), cameraUp);

          // Rotation verticale (autour de l'axe X)
          glm::vec3 right = glm::normalize(glm::cross(cameraTarget - cameraPosition, cameraUp));
          glm::mat4 rotationX = glm::rotate(glm::mat4(1.0f), glm::radians(-dy * sensitivity), right);

          // Appliquer les rotations à la direction de la caméra
          glm::vec3 direction = cameraPosition - cameraTarget;
          direction = glm::vec3(rotationX * rotationY * glm::vec4(direction, 0.0f));

          // Mettre à jour la position de la caméra
          cameraPosition = cameraTarget + direction;

          // Mettre à jour la matrice de vue
          updateCameraView();
        } else if ((event.motion.state & SDL_BUTTON_RMASK)) {
          // Zoom avec bouton droit
          float zoomSensitivity = 0.01f;
          float zoomFactor = 1.0f - dy * zoomSensitivity;
          cameraDistance *= zoomFactor;

          // Limiter la distance de la caméra
          cameraDistance = glm::clamp(cameraDistance, cameraDistanceMin, cameraDistanceMax);

          // Mettre à jour la position de la caméra
          glm::vec3 direction = glm::normalize(cameraPosition - cameraTarget);
          cameraPosition = cameraTarget + direction * cameraDistance;

          // Mettre à jour la matrice de vue
          updateCameraView();
        } else if (shiftPressed && (event.motion.state & SDL_BUTTON_LMASK)) {
          // Pan (translation) avec Shift + bouton gauche
          float panSensitivity = 0.005f;

          // Calculer les vecteurs de déplacement
          glm::vec3 direction = glm::normalize(cameraPosition - cameraTarget);
          glm::vec3 right = glm::normalize(glm::cross(direction, cameraUp));
          glm::vec3 up = glm::normalize(cameraUp);

          // Déplacer la caméra et la cible dans le plan perpendiculaire à la
          // direction de vue Correction : inverser le signe pour dy pour avoir un
          // comportement intuitif
          glm::vec3 pan = right * (dx * panSensitivity) + up * (dy * panSensitivity);
          cameraPosition += pan;
          cameraTarget += pan;

          // Mettre à jour la matrice de vue
          updateCameraView();
        }

        // Mettre à jour les coordonnées précédentes de la souris
        lastMouseX = event.motion.x;
        lastMouseY = event.motion.y;
      }
      break;

    case SDL_MOUSEWHEEL: {
      // Zoom avec la molette de la souris
      float zoomFactor = 1.0f - event.wheel.y * 0.1f;
      cameraDistance *= zoomFactor;

      // Limiter la distance de la caméra
      cameraDistance = glm::clamp(cameraDistance, cameraDistanceMin, cameraDistanceMax);
      // cameraDistance = glm::clamp(cameraDistance, 1.0f, 20.0f);
      // if (cameraDistance < znear) cameraDistance = znear;
      // if (cameraDistance > zfar) cameraDistance = zfar;

      // Mettre à jour la position de la caméra
      glm::vec3 direction = glm::normalize(cameraPosition - cameraTarget);
      cameraPosition = cameraTarget + direction * cameraDistance;

      // Mettre à jour la matrice de vue
      updateCameraView();
    } break;
  }
}

void Scene3D::applyProjectionAndView() {
  // Appliquer la matrice de projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMultMatrixf(glm::value_ptr(projectionMatrix));

  // Appliquer la matrice de vue
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixf(glm::value_ptr(viewMatrix));
}

void Scene3D::render() {
  // first are displayed after

  if (particleColoring.show_colorbar) {
    particleColoring.draw_colorbar(box, windowWidth, windowHeight);
  }

  if (show_globalFrame) {
    drawGlobalAABBFrame();
  }

  if (show_obb) {
    drawOBBs();
  }

  if (show_particles || show_driven) {
    drawParticles();
  }
}

void Scene3D::drawShape(Shape* s, double homothety) {
  double R = homothety * s->radius;
  int nbLevelSphere = 2;

  if (1) {
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
    if (s->face[f].size() < 3) {  // At least 3 pts!
      continue;
    }
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

void Scene3D::drawGlobalAABB() {
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  OBB obb;
  obb.extent = 0.5 * (box.aabb.max - box.aabb.min);
  obb.center = 0.5 * (box.aabb.min + box.aabb.max);
  glColor3f(0.0f, 0.0f, 0.0f);
  glShape::obb(obb);
}

void Scene3D::drawGlobalAABBFrame() {
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  OBB obb;
  obb.extent = 0.5 * (box.aabb.max - box.aabb.min);
  double extMin = std::min({obb.extent.x, obb.extent.y, obb.extent.z});

  // X
  glColor3f(1.0f, 0.0f, 0.0f);
  vec3r arrowX(2.0 * extMin, 0.0, 0.0);
  glShape::arrow(-0.5 * arrowX, arrowX, -1.0, 0.7);

  // Y
  glColor3f(0.0f, 1.0f, 0.0f);
  vec3r arrowY(0.0, 2.0 * extMin, 0.0);
  glShape::arrow(-0.5 * arrowY, arrowY, -1.0, 0.7);

  // Z
  glColor3f(0.0f, 0.0f, 1.0f);
  vec3r arrowZ(0.0, 0.0, 2.0 * extMin);
  glShape::arrow(-0.5 * arrowZ, arrowZ, -1.0, 0.7);
}

void Scene3D::drawParticlePositions() {
  glColor3f(0.0f, 0.0f, 0.0f);
  glPointSize(2);
  glBegin(GL_POINTS);
  for (size_t i = box.nDriven; i < box.Particles.size(); ++i) {

    vec3r pos = box.Particles[i].pos;
    glVertex3f(pos.x, pos.y, pos.z);
  }
  glEnd();
}

void Scene3D::drawOBBs() {
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

void Scene3D::drawParticles() {
  if (mouseDragging == true) {
    drawParticlePositions();
    drawGlobalAABB();
    return;
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  if (show_particles) {
    
    // Add OBB for the selected particle
    if (selectedParticle >= 0) {
      OBB obb = box.Particles[selectedParticle].shape->obb;
      obb.rotate(box.Particles[selectedParticle].Q);
      obb.extent *= box.Particles[selectedParticle].homothety;
      obb.center *= box.Particles[selectedParticle].homothety;
      obb.center += box.Particles[selectedParticle].pos;
      glColor4ub(234, 255, 0, (int)floor(alphaParticles * 255));  // shiny yellow
      glShape::obb(obb, 4.0f);
      glEnable(GL_LIGHTING);
      glEnable(GL_DEPTH_TEST);
    }
    
    for (size_t i = box.nDriven; i < box.Particles.size(); ++i) {

      if (selectedParticle >= 0 && i == (size_t)selectedParticle) {
        glColor4ub(234, 255, 0, (int)floor(alphaParticles * 255));  // shiny yellow
      } else {
        if (i >= particleColoring.pcolors.size()) {
          glColor4f(colorParticles.r, colorParticles.g, colorParticles.b, alphaParticles);
        } else {
          glColor4ub(particleColoring.pcolors[i].r, particleColoring.pcolors[i].g, particleColoring.pcolors[i].b,
                     (int)floor(alphaParticles * 255));
        }
      }

      vec3r pos = box.Particles[i].pos;
      if (use_slice && !slice.inSlice(pos)) continue;


      glPushMatrix();
      glTranslatef(pos.x, pos.y, pos.z);
      quat2GLMatrix<GLfloat>(box.Particles[i].Q, Rot_Matrix);
      glMultMatrixf(Rot_Matrix);
      drawShape(box.Particles[i].shape, box.Particles[i].homothety);
      glPopMatrix();
    }
  }

  if (show_driven == 1) {
    for (size_t i = 0; i < box.nDriven; ++i) {

      if (selectedParticle >= 0 && i == (size_t)selectedParticle) {
        glColor4ub(234, 255, 0, (int)floor(alphaDrivens * 255));  // shiny yellow
      } else {
        glColor4f(colorDrivens.r, colorDrivens.g, colorDrivens.b, alphaDrivens);
      }

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
