#include "Scene3D.hpp"
#include "glTools.hpp"

Scene3D::Scene3D() {
  //////////
  Shapes.clear();  // je ne sais pas vraiment pourquoi je suis obligé de faire ça
  Shape fake;
  Shapes.push_back(fake);
  ////////////

  cameraDistance = glm::distance(cameraPosition, cameraTarget);
  cameraUp = glm::normalize(cameraUp);


  updateProjection();
  updateCameraView();
}

int Scene3D::readShapeLib(const char* fileName) {
  if (!fileTool::fileExists(fileName)) {
    std::cout << "Shape Library named '" << fileName << "' has not been found." << std::endl;
    return 0;
  }
  Shapes.clear();
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

  fitView();
  // OBB& obb = Shapes[ishape].obb;

  // center.set(obb.center.x, obb.center.y, obb.center.y);  // where we look at
  // eye.set(obb.center.x + obb.extent.x, obb.center.y,
  //         obb.center.y);  // from where we look
  // up.set(0.0, 1.0, 0.0);  // direction (normalized)

  return 1;
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

  OBB& obb = Shapes[ishape].obb;
  vec3r mx = 2.4 * (obb.extent[0] * obb.e[0] + obb.extent[1] * obb.e[1] + obb.extent[2] * obb.e[2]);

  float diag_length = (GLfloat)(2 * norm(mx));
  float min_length = std::min({mx.x, mx.y, mx.z});
  cameraDistanceMin = min_length * 0.0001f;
  cameraDistanceMax = diag_length * 2.5f;

  znear = zf - diag_length;
  double close_dst = 0.1 * zf;
  if (znear < close_dst) znear = close_dst;
  zfar = zf + diag_length;

  projectionMatrix = glm::perspective(glm::radians(cameraFov), wh_ratio, znear, zfar);
}

void Scene3D::updateCameraView() { viewMatrix = glm::lookAt(cameraPosition, cameraTarget, cameraUp); }

void Scene3D::fitView() {
  glm::vec3 direction = glm::normalize(cameraPosition - cameraTarget);
  OBB& obb = Shapes[ishape].obb;
  float maxDimension = 2.4 * std::max({obb.extent.x, obb.extent.y, obb.extent.z});

  cameraTarget.x = obb.center.x;
  cameraTarget.y = obb.center.y;
  cameraTarget.z = obb.center.z;
  cameraDistance = 0.5f * maxDimension / (float)(tan(glm::radians(0.5f * cameraFov)));
  cameraPosition = cameraTarget + cameraDistance * direction;
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
        // selection(event.button.x, event.button.y);
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

  if (show_frame) drawFrame();
  if (show_obb) drawOBB(ishape);

  glEnable(GL_LIGHTING);
  drawShape(ishape);
}

void Scene3D::drawOBB(size_t ishp) {

  if (ishp >= Shapes.size()) {
    return;
  }

  // glDisable(GL_LIGHTING);
  glColor3f(0.2f, 0.2f, 0.8f);
  glShape::obb(Shapes[ishp].obb);
}

void Scene3D::drawFrame() {
  OBB& obb = Shapes[ishape].obb;
  vec3r diag = 2.0 * (obb.extent[0] * obb.e[0] + obb.extent[1] * obb.e[1] + obb.extent[2] * obb.e[2]);
  double len = diag.length() * 0.333;

  glColor3f(1.0f, 0.0f, 0.0f);
  glShape::arrow(vec3r::zero(), len * vec3r::unit_x());
  glColor3f(0.0f, 1.0f, 0.0f);
  glShape::arrow(vec3r::zero(), len * vec3r::unit_y());
  glColor3f(0.0f, 0.0f, 1.0f);
  glShape::arrow(vec3r::zero(), len * vec3r::unit_z());
}

void Scene3D::drawShape(size_t ishp) {
  if (ishp >= Shapes.size()) {
    return;
  }

  double R = Shapes[ishp].radius;
  // glColor4f(0.666f, 0.729f, 0.09f, alpha); // yellow
  glColor4f(colorParticles.x, colorParticles.y, colorParticles.z, alphaParticles);

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

/*
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
*/
