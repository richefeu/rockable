// Classe principale pour une application ImGui + SDL + OpenGL
#pragma once

#include <GL/glew.h>
#include <SDL2/SDL.h>

#include "toofus-gate/glm/glm.hpp"
#include "toofus-gate/glm/gtc/matrix_transform.hpp"

#include <imgui.h>
#include <imgui_impl_opengl3.h>
#include <imgui_impl_sdl2.h>
#include <implot.h>

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "Scene3D.hpp"

/// Classe principale de l'application gérant SDL, OpenGL et ImGui
class ImGuiSDLApplication {
 public:
  ImGuiSDLApplication(int width, int height, const std::string& title);
  ~ImGuiSDLApplication();

  bool init();

  int run();

  Rockable& box();
  int confNum{0};
  int lastConfNumOK{confNum};

  bool tryToReadConf(int num);

  Scene3D scene3d;

  // Empêcher la copie pour éviter les problèmes de gestion des ressources
  ImGuiSDLApplication(const ImGuiSDLApplication&) = delete;
  ImGuiSDLApplication& operator=(const ImGuiSDLApplication&) = delete;

 private:
  int windowWidth;
  int windowHeight;
  std::string windowTitle;

  SDL_Window* window{nullptr};
  SDL_GLContext glContext{nullptr};

  bool show_config_panel{false};
  bool show_main_panel{true};
  bool show_window_info{false};
  bool show_screenshot_panel{false};
  bool show_slicing_panel{false};

  // screenshots
  bool should_save_screenshot{false};
  bool should_save_screenshot_series{false};
  char screenshot_filename[256];
  char screenshot_basename[64] = {"screenshot"};
  int screenshot_first{0};
  int screenshot_last{0};

  // background
  ImVec4 clear_color = ImVec4(0.498f, 0.655f, 0.721f, 1.0f);
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
    // std::string description;
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
  bool SaveScreenshotWithoutImGui(const char* filename, int W, int H);

  // ImGui helper functions
  static void HelpPopup(const char* desc);
  void Input_vec3r(const char* label, vec3r& myVec);
  void Input_quat(const char* label, quat& myQuat);
  void Input_Shape(Particle& P);
  void reset_ngroup();
  void Input_Properties(const char* propertyName);
  void Input_LawData(const char* parName);
  void Input_ParticleColorMode();

  void ForceLaw_Combo();
  void Integrator_Combo();
  void UpdateNL_Combo();
  void AddOrRemoveInteractions_Combo();

  // Panels
  void mainMenu();
  void mainPanel();
  void infoPanel();
  void configPanel();
  void screenshotPanel();
  void slicingPanel();

  // Boucle principale ============
  void mainLoop();
  void handleEvents();
  void renderFrame();

  // Configurations d'affichage ============
  void applyConfigToScene();
  void saveConfigToSeparateFile();
  void loadConfigFromSeparateFile();

  void cleanup();  // Nettoie toutes les ressources

 private:
  // Du graphisme sur les buttons
  void LoadTextureFromMemory(unsigned char data[], unsigned int len, GLuint& textureID);
  void genButtons();

  GLuint IDButton_XY;
  GLuint IDButton_XZ;
  GLuint IDButton_YX;
  GLuint IDButton_YZ;
  GLuint IDButton_ZX;
  GLuint IDButton_ZY;
};
