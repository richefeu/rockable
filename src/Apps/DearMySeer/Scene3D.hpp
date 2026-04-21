/// Structure de base pour la gestion d'une scène 3D avec OpenGL

#pragma once

#include "Core/Rockable.hpp"

#include "toofus-gate/glm/glm.hpp"
#include "toofus-gate/glm/gtc/matrix_transform.hpp"
#include "toofus-gate/glm/gtc/type_ptr.hpp"

#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <imgui.h>
#include <iostream>

#include "ParticleColoring.hpp"
#include "Slicing.hpp"

// Déclaration anticipée pour éviter d'inclure SDL.h dans l'en-tête
union SDL_Event;

/// Gère une scène 3D avec projection, vue et rendu
struct Scene3D {

  Rockable box;
  int confNum = 0;
  
  ParticleColoring particleColoring;
  std::string particleColorMode{"none"};
  
  Slicing slice;
  int use_slice{0};
  
  float rotationAngle;  ///< Angle de rotation actuel (en degrés)
  bool useWireframe;    ///< Mode wireframe activé/désactivé
  int windowWidth;      ///< Largeur actuelle de la fenêtre
  int windowHeight;     ///< Hauteur actuelle de la fenêtre
  float wh_ratio{(float)windowWidth / (float)windowHeight};
  glm::mat4 projectionMatrix;  ///< Matrice de projection
  glm::mat4 viewMatrix;        ///< Matrice de vue (caméra)

  // flags with default values
  int show_background{1};
  int show_globalFrame{1};
  int show_particles{1};
  int show_driven{1};
  int show_velocities{0};
  int show_slice{0};
  int show_forces{0};
  int show_obb{0};
  int enlarged_obb{0};
  int link_obb{1};
  int show_interFrames{0};
  int show_interTypes{0};
  int show_keybinds{0};
  int show_traj{0};
  int show_probe{0};

  glm::vec3 colorParticles{0.631f, 0.259f, 0.961f};
  float alphaParticles{1.0f};

  glm::vec3 colorDrivens{0.808f, 0.843f, 0.631f};
  float alphaDrivens{0.2f};

  // Paramètres de caméra courants (avec valeurs par défaut pour la réinitialisation)
  glm::vec3 cameraPosition = glm::vec3(0.0f, 0.0f, 1.0f);  ///< Position de la caméra
  glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);    ///< Point cible de la caméra
  glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);        ///< Vecteur haut de la caméra
  float cameraDistance{16.0f};                             ///< Distance de la caméra à la cible
  float cameraDistanceMin{0.01f};
  float cameraDistanceMax{20.0f};
  float cameraFov{45.0f};                                  ///< Champ de vision de la caméra (en degrés)
  float znear{0.001f};
  float zfar{50.0f};

  int lastMouseX;      ///< Position X précédente de la souris
  int lastMouseY;      ///< Position Y précédente de la souris
  bool mouseDragging;  ///< Indique si la souris est en train de glisser

  int selectedParticle{-1};

  // Paramètres du repère 3D
  bool showAxis;           ///< Afficher le repère 3D au lieu du cube
  float axisSize;          ///< Taille du repère 3D
  glm::vec3 xAxisColor;    ///< Couleur de l'axe X
  glm::vec3 yAxisColor;    ///< Couleur de l'axe Y
  glm::vec3 zAxisColor;    ///< Couleur de l'axe Z
  int arrowFaces;          ///< Nombre de faces pour les flèches (par défaut 16)
  bool useLighting;        ///< Activer l'éclairage pour les flèches
  float arrowLengthRatio;  ///< Rapport de longueur de la flèche par rapport à la taille de l'axe
  float arrowRadiusRatio;  ///< Rapport de rayon de la flèche par rapport à sa longueur

  GLfloat Rot_Matrix[16];

  Scene3D();
  bool init();
  void updateProjection();
  void updateCameraView();
  void fitView();
  void applyProjectionAndView();

  void render();
  void selection(int x, int y);
  void handleMouseEvent(const SDL_Event& event, int windowWidth, int windowHeight);

  void drawGlobalAABB();
  void drawGlobalAABBFrame();
  void drawShape(Shape* s, double homothety);
  void drawParticlePositions();
  void drawParticles();
  void drawOBBs();
};
