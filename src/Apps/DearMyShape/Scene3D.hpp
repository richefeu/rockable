/// Structure de base pour la gestion d'une scène 3D avec OpenGL

#pragma once

#include "Core/Shape.hpp"

#include "toofus-gate/glm/glm.hpp"
#include "toofus-gate/glm/gtc/matrix_transform.hpp"
#include "toofus-gate/glm/gtc/type_ptr.hpp"

#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <imgui.h>
#include <iostream>
#include <string>
#include <vector>

#include "fileTool.hpp"

// Déclaration anticipée pour éviter d'inclure SDL.h dans l'en-tête
union SDL_Event;

/// Gère une scène 3D avec projection, vue et rendu
struct Scene3D {

  size_t ishape{0};  // id of the current shape beeing shown
  std::vector<Shape> Shapes;

  std::string shapeFileName{"shapes.txt"};

  int windowWidth{800};       ///< Largeur actuelle de la fenêtre
  int windowHeight{800};      ///< Hauteur actuelle de la fenêtre
  float wh_ratio{(float)windowWidth / (float)windowHeight};
  glm::mat4 projectionMatrix;  ///< Matrice de projection
  glm::mat4 viewMatrix;        ///< Matrice de vue (caméra)

  // flags with default values
  int show_background{1};
  int show_obb{1};
  int show_frame{1};

  glm::vec3 colorParticles{0.6f, 0.6f, 0.6f};
  float alphaParticles{1.0f};

  // Paramètres de caméra courants (avec valeurs par défaut pour la réinitialisation)
  glm::vec3 cameraPosition = glm::vec3(1.0f, 1.0f, 1.0f);
  glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);
  //glm::vec3 cameraUp = glm::vec3(-0.57735f, -0.57735f, 0.57735f);
  glm::vec3 cameraUp = glm::vec3(0.0f, 0.0f, 1.0f);
  float cameraDistance{16.0f};
  float cameraDistanceMin{0.01f};
  float cameraDistanceMax{20.0f};
  float cameraFov{45.0f};
  float znear{0.001f};
  float zfar{50.0f};

  int lastMouseX{0};
  int lastMouseY{0};
  bool mouseDragging{false};

  GLfloat Rot_Matrix[16];

  Scene3D();
  bool init();
  void updateProjection();
  void updateCameraView();
  void fitView();
  void applyProjectionAndView();

  void render();
  void handleMouseEvent(const SDL_Event& event, int windowWidth, int windowHeight);

  void drawFrame();
  void drawShape(size_t ishp);
  void drawOBB(size_t ishp);

  int readShapeLib(const char* fileName);
  // void saveShapeLib(const char* fileName);
  // void exportSample();
};
