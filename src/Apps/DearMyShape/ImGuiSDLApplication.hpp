// Classe principale pour une application ImGui + SDL + OpenGL
#pragma once

#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <functional>
#include "toofus-gate/glm/glm.hpp"
#include "toofus-gate/glm/gtc/matrix_transform.hpp"
#include <imgui.h>
#include <imgui_impl_opengl3.h>
#include <imgui_impl_sdl2.h>
#include <implot.h>
#include <memory>
#include <string>
#include <vector>

#include "Scene3D.hpp"

#define ADD_DEMO 0

/// Classe principale de l'application gérant SDL, OpenGL et ImGui
class ImGuiSDLApplication {
 public:
  ImGuiSDLApplication(int width, int height, const std::string& title);
  ~ImGuiSDLApplication();

  bool init();

  int run();

  Scene3D scene3d;

  // Empêcher la copie pour éviter les problèmes de gestion des ressources
  ImGuiSDLApplication(const ImGuiSDLApplication&) = delete;
  ImGuiSDLApplication& operator=(const ImGuiSDLApplication&) = delete;

 private:
  int windowWidth;          ///< Largeur de la fenêtre
  int windowHeight;         ///< Hauteur de la fenêtre
  std::string windowTitle;  ///< Titre de la fenêtre

  SDL_Window* window{nullptr};       ///< Fenêtre SDL
  SDL_GLContext glContext{nullptr};  ///< Contexte OpenGL

  // bool show_demo_window{false};  ///< Afficher la fenêtre de démonstration ImGui
  // bool show_config_panel{false};
  bool show_main_panel{true};
  // bool show_window_info{false};
  // bool show_screenshot_panel{false};
  // bool show_slicing_panel{false};

  // background
  ImVec4 clear_color = ImVec4(0.498f, 0.655f, 0.721f, 1.0f);  ///< Couleur de fond
  float clear_bottom_color[3] = {135.0f / 255.0f, 206.0f / 255.0f, 250.0f / 255.0f};
  float clear_top_color[3] = {1.0f, 1.0f, 1.0f};
  void gradBackground();

  // Gestion des raccourcis clavier ============
  struct KeyBinding {
    SDL_Keycode key;
    bool ctrl;
    bool shift;
    bool alt;
    std::function<void()> action;
  };

  std::vector<KeyBinding> keyBindings;  ///< Liste des raccourcis clavier configurés
  bool shouldQuit{false};               ///< Indicateur de fermeture de l'application

  void setupDefaultKeyBindings();
  void addKeyBinding(SDL_Keycode key, bool ctrl, bool shift, bool alt, std::function<void()> action);
  void checkKeyBindings(const SDL_Event& event);

  // Méthodes d'initialisation ============
  bool initializeSDL();
  bool initializeOpenGL();
  bool initializeImGui();
  bool applyStyle();
  void setupScene3D();

  // screenshots
  // bool SaveScreenshotWithoutImGui(const char* filename, int W, int H);

  // ImGui helper functions
  static void HelpPopup(const char* desc);
  void Input_vec3r(const char* label, vec3r& myVec);
  void Input_quat(const char* label, quat& myQuat);
  // void Input_Shape(Particle& P);

  // Panels
  void mainMenu();
  void mainPanel();

  // Boucle principale ============
  void mainLoop();
  void handleEvents();
  void renderFrame();

  // Configurations d'affichage ============
  void applyConfigToScene();
  void saveConfigToSeparateFile();
  void loadConfigFromSeparateFile();

  void cleanup();  // Nettoie toutes les ressources
};
