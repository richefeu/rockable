#include "ImGuiSDLApplication.hpp"

#include <fstream>
#include <iostream>
#include <string>

ImGuiSDLApplication::ImGuiSDLApplication(int width, int height, const std::string& title)
    : windowWidth(width), windowHeight(height), windowTitle(title) {

  setupDefaultKeyBindings();
}

ImGuiSDLApplication::~ImGuiSDLApplication() { cleanup(); }

// =======================
// Initializations and Interface configuration
// =======================

// Initialise tous les sous-systèmes de l'application
bool ImGuiSDLApplication::init() {
  if (!initializeSDL()) {
    return false;
  }

  if (!initializeOpenGL()) {
    return false;
  }

  if (!initializeImGui()) {
    return false;
  }

  setupScene3D();

  return true;
}

// Initialise SDL et crée la fenêtre principale
bool ImGuiSDLApplication::initializeSDL() {
  // Initialisation de SDL
  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0) {
    std::cerr << "Erreur SDL_Init: " << SDL_GetError() << std::endl;
    return false;
  }

  // Configuration des attributs OpenGL pour le mode compatible (fonctions fixes)
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_COMPATIBILITY);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  // Création de la fenêtre
  window = SDL_CreateWindow(windowTitle.c_str(), SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, windowWidth,
                            windowHeight, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);

  if (!window) {
    std::cerr << "Erreur SDL_CreateWindow: " << SDL_GetError() << std::endl;
    SDL_Quit();
    return false;
  }

  return true;
}

// Initialise OpenGL et GLEW
bool ImGuiSDLApplication::initializeOpenGL() {
  // Création du contexte OpenGL
  glContext = SDL_GL_CreateContext(window);
  if (!glContext) {
    std::cerr << "Erreur SDL_GL_CreateContext: " << SDL_GetError() << std::endl;
    SDL_DestroyWindow(window);
    SDL_Quit();
    return false;
  }

  // Activation de la synchronisation verticale
  if (SDL_GL_SetSwapInterval(1) < 0) {
    std::cerr << "Warning: Unable to set VSync: " << SDL_GetError() << std::endl;
  }

  // Initialisation de GLEW
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {
    std::cerr << "Erreur glewInit" << std::endl;
    SDL_GL_DeleteContext(glContext);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return false;
  }

  return true;
}

// Initialise ImGui et ses bindings SDL/OpenGL
bool ImGuiSDLApplication::initializeImGui() {
  // Initialisation d'ImGui
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  (void)io;

  // Configuration du style
  ImGui::StyleColorsClassic();
  if (!applyStyle()) {
    std::cerr << "Erreur applyStyle" << std::endl;
    return false;
  }

  // Configuration des plateformes et des renderers
  if (!ImGui_ImplSDL2_InitForOpenGL(window, glContext)) {
    std::cerr << "Erreur ImGui_ImplSDL2_InitForOpenGL" << std::endl;
    return false;
  }

  if (!ImGui_ImplOpenGL3_Init("#version 120")) {
    std::cerr << "Erreur ImGui_ImplOpenGL3_Init" << std::endl;
    return false;
  }

  ImPlot::CreateContext();

  return true;
}

bool ImGuiSDLApplication::applyStyle() {
  ImGuiStyle& style = ImGui::GetStyle();

  style.Colors[ImGuiCol_Text] = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TextDisabled] = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  style.Colors[ImGuiCol_WindowBg] = ImVec4(0.94f, 0.94f, 0.94f, 0.94f);
  style.Colors[ImGuiCol_PopupBg] = ImVec4(1.00f, 1.00f, 1.00f, 0.94f);
  style.Colors[ImGuiCol_Border] = ImVec4(0.00f, 0.00f, 0.00f, 0.39f);
  style.Colors[ImGuiCol_BorderShadow] = ImVec4(1.00f, 1.00f, 1.00f, 0.10f);
  style.Colors[ImGuiCol_FrameBg] = ImVec4(1.00f, 1.00f, 1.00f, 0.94f);
  style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
  style.Colors[ImGuiCol_FrameBgActive] = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
  style.Colors[ImGuiCol_TitleBg] = ImVec4(0.96f, 0.96f, 0.96f, 1.00f);
  style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(1.00f, 1.00f, 1.00f, 0.51f);
  style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.82f, 0.82f, 0.82f, 1.00f);
  style.Colors[ImGuiCol_MenuBarBg] = ImVec4(0.86f, 0.86f, 0.86f, 1.00f);
  style.Colors[ImGuiCol_ScrollbarBg] = ImVec4(0.98f, 0.98f, 0.98f, 0.53f);
  style.Colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.69f, 0.69f, 0.69f, 1.00f);
  style.Colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.59f, 0.59f, 0.59f, 1.00f);
  style.Colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.49f, 0.49f, 0.49f, 1.00f);
  style.Colors[ImGuiCol_CheckMark] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.24f, 0.52f, 0.88f, 1.00f);
  style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_Button] = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
  style.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_ButtonActive] = ImVec4(0.06f, 0.53f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_Header] = ImVec4(0.26f, 0.59f, 0.98f, 0.31f);
  style.Colors[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.80f);
  style.Colors[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_ResizeGrip] = ImVec4(1.00f, 1.00f, 1.00f, 0.50f);
  style.Colors[ImGuiCol_ResizeGripHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
  style.Colors[ImGuiCol_ResizeGripActive] = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
  style.Colors[ImGuiCol_PlotLines] = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
  style.Colors[ImGuiCol_PlotLinesHovered] = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogram] = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogramHovered] = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TextSelectedBg] = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);

  style.WindowRounding = 4;
  style.ChildRounding = 4;
  style.FrameRounding = 3;
  style.PopupRounding = 4;
  style.ScrollbarRounding = 9;
  style.GrabRounding = 2;
  style.TabRounding = 4;
  style.FrameBorderSize = 1;
  style.ScrollbarSize = 20;
  style.Alpha = 0.9f;

  return true;
}

void ImGuiSDLApplication::setupScene3D() {
  scene3d.windowWidth = windowWidth;
  scene3d.windowHeight = windowHeight;
  scene3d.updateProjection();

  if (!scene3d.init()) {
    std::cerr << "Échec de l'initialisation de la scène 3D" << std::endl;
  }

  // Charger la configuration depuis le fichier INI
  // loadConfigFromIni();
  // loadConfigFromSeparateFile();
  applyConfigToScene();
}

/// @brief Applique la configuration à la scène 3D
void ImGuiSDLApplication::applyConfigToScene() {
  // Les paramètres sont déjà chargés directement dans scene3d par
  // loadConfigFromSeparateFile() Nous devons juste nous assurer que la caméra
  // est correctement mise à jour

  //  Mettre à jour la distance de la caméra (au cas où la position ou la cible
  //  aurait changé)
  scene3d.cameraDistance = glm::distance(scene3d.cameraPosition, scene3d.cameraTarget);

  // Mettre à jour les matrices de projection et de vue
  scene3d.updateProjection();
  scene3d.updateCameraView();
}

/// @brief Sauvegarde la configuration dans un fichier séparé
void ImGuiSDLApplication::saveConfigToSeparateFile() {}

/// @brief Charge la configuration depuis un fichier séparé
void ImGuiSDLApplication::loadConfigFromSeparateFile() {}

// =======================
// Key Binding
// =======================

/// @brief Configure les raccourcis clavier par défaut
void ImGuiSDLApplication::setupDefaultKeyBindings() {
  addKeyBinding(SDLK_q, true, false, false, [this]() { shouldQuit = true; });

  addKeyBinding(SDLK_m, true, false, false, [this]() { show_main_panel = !show_main_panel; });
  addKeyBinding(SDLK_ESCAPE, false, false, false, [this]() { shouldQuit = true; });
}

/// Ajoute un raccourci clavier à la liste
void ImGuiSDLApplication::addKeyBinding(SDL_Keycode key, bool ctrl, bool shift, bool alt,
                                        std::function<void()> action) {
  keyBindings.push_back({key, ctrl, shift, alt, action});
}

/// Vérifie si un raccourci clavier doit être exécuté
void ImGuiSDLApplication::checkKeyBindings(const SDL_Event& event) {
  if (event.type != SDL_KEYDOWN) {
    return;
  }

  SDL_Keymod mods = SDL_GetModState();
  bool ctrlPressed = (mods & KMOD_CTRL) != 0;
  bool shiftPressed = (mods & KMOD_SHIFT) != 0;
  bool altPressed = (mods & KMOD_ALT) != 0;

  for (const auto& binding : keyBindings) {
    if (event.key.keysym.sym == binding.key && ctrlPressed == binding.ctrl && shiftPressed == binding.shift &&
        altPressed == binding.alt) {

      binding.action();
      break;  // Un seul raccourci par événement
    }
  }
}

// =======================
// imgui inputs
// =======================

void ImGuiSDLApplication::HelpPopup(const char* desc) {
  ImGui::TextDisabled("(?)");
  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
    ImGui::TextUnformatted(desc);
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
  }
}

void ImGuiSDLApplication::Input_vec3r(const char* label, vec3r& myVec) {
  vec3<float> tmp(myVec.x, myVec.y, myVec.z);
  if (ImGui::InputFloat3(label, &(tmp.x))) {
    myVec.set(tmp.x, tmp.y, tmp.z);
  }
}

void ImGuiSDLApplication::Input_quat(const char* label, quat& myQuat) {
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

/*
void ImGuiSDLApplication::Input_Shape(Particle& P) {
  size_t item_current_idx = box().shapeId[P.shape->name];
  if (ImGui::BeginCombo("shape", P.shape->name.c_str(), 0)) {
    for (size_t n = 0; n < box().Shapes.size(); n++) {
      const bool is_selected = (item_current_idx == n);
      if (ImGui::Selectable(box().Shapes[n].name.c_str(), is_selected)) {
        item_current_idx = n;
        P.shape = &(box().Shapes[n]);  // Plug to the selected shape
      }
      if (is_selected) {
        ImGui::SetItemDefaultFocus();
      }
    }
    ImGui::EndCombo();
  }
}
*/

// =======================
// Main menu and panels
// =======================

void ImGuiSDLApplication::mainMenu() {
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem("Quit", "esc")) {
        shouldQuit = true;
      }
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Panels")) {
      ImGui::MenuItem("Main menu", "ctrl-m", &show_main_panel);
      ImGui::EndMenu();
    }

    ImGui::EndMainMenuBar();
  }
}

void ImGuiSDLApplication::mainPanel() {

  ImGui::Begin("Main menu", &show_main_panel, ImGuiWindowFlags_AlwaysAutoResize);
  ImGui::PushItemWidth(200);

  ImGui::Text("Number of shapes : %zu", scene3d.Shapes.size());

  int N = static_cast<int>(scene3d.Shapes.size());
  if (N > 0) N -= 1;
  int temp_ishape = static_cast<int>(scene3d.ishape);
  if (ImGui::SliderInt("ishape", &temp_ishape, 0, N)) {
    scene3d.ishape = static_cast<size_t>(temp_ishape);
    scene3d.fitView();
  }

  // ImGui::Separator();

  if (ImGui::CollapsingHeader("Show")) {
    ImGui::CheckboxFlags("Background", &scene3d.show_background, 1);
    ImGui::CheckboxFlags("OBB", &scene3d.show_obb, 1);
    // ...
    ImGui::SliderFloat("Transparency", &scene3d.alphaParticles, 0.0f, 1.0f);
  }

  if (ImGui::CollapsingHeader("Camera")) {
    if (ImGui::Button("Fit")) {
      scene3d.fitView();
    }
    ImGui::SameLine();
    if (ImGui::Button("Align")) {
      glm::vec3 d = scene3d.cameraTarget - scene3d.cameraPosition;

      scene3d.cameraUp.x = -d.x * d.y;
      scene3d.cameraUp.y = d.x * d.x + d.z * d.z;
      scene3d.cameraUp.z = -d.z * d.y;
      // Je ne me souviens plus de la logique ici (à vérifier)

      (void)glm::normalize(scene3d.cameraUp);
    }

    ImGui::SliderFloat("FOV", &scene3d.cameraFov, 30.0f, 120.0f, "%.1f°");
  }

  if (ImGui::CollapsingHeader("Colors")) {
    ImGui::PushItemWidth(200);
    ImGui::ColorEdit3("Background bottom", (float*)&clear_bottom_color);
    ImGui::ColorEdit3("Background top", (float*)&clear_top_color);
    ImGui::PopItemWidth();
  }

  // ImGui::Separator();

  if (ImGui::CollapsingHeader("Oriented Bounding Box")) {

    static std::vector<std::string> modes = {"Covariance", "Min volume", "Axis aligned", "Imposed axis"};
    static size_t item_current_idx = 0;
    ImGui::SetNextItemWidth(200.0f);
    if (ImGui::BeginCombo("fitObbOption", modes[scene3d.Shapes[scene3d.ishape].fibObbOption].c_str(), 0)) {
      for (size_t n = 0; n < modes.size(); n++) {
        const bool is_selected = ((size_t)item_current_idx == n);
        if (ImGui::Selectable(modes[n].c_str(), is_selected)) {
          scene3d.Shapes[scene3d.ishape].fibObbOption = (size_t)n;
        }
        if (is_selected) {
          ImGui::SetItemDefaultFocus();
        }
      }
      ImGui::EndCombo();
    }

    if (ImGui::Button("fitObb")) {
      scene3d.Shapes[scene3d.ishape].fitObb();
    }

    ImGui::Separator();
    ImGui::Text("obb.center: %g %g %g", scene3d.Shapes[scene3d.ishape].obb.center.x,
                scene3d.Shapes[scene3d.ishape].obb.center.y, scene3d.Shapes[scene3d.ishape].obb.center.z);
    ImGui::Text("obb.e0: %g %g %g", scene3d.Shapes[scene3d.ishape].obb.e[0].x,
                scene3d.Shapes[scene3d.ishape].obb.e[0].y, scene3d.Shapes[scene3d.ishape].obb.e[0].z);
    ImGui::Text("obb.e0: %g %g %g", scene3d.Shapes[scene3d.ishape].obb.e[1].x,
                scene3d.Shapes[scene3d.ishape].obb.e[1].y, scene3d.Shapes[scene3d.ishape].obb.e[1].z);
    ImGui::Text("obb.e0: %g %g %g", scene3d.Shapes[scene3d.ishape].obb.e[2].x,
                scene3d.Shapes[scene3d.ishape].obb.e[2].y, scene3d.Shapes[scene3d.ishape].obb.e[2].z);
    ImGui::Text("obb.extent: %g %g %g", scene3d.Shapes[scene3d.ishape].obb.extent.x,
                scene3d.Shapes[scene3d.ishape].obb.extent.y, scene3d.Shapes[scene3d.ishape].obb.extent.z);
  }

  if (ImGui::CollapsingHeader("Mass Properties")) {
    int pre_computation_is_done = (scene3d.Shapes[scene3d.ishape].preCompDone == 'y') ? 1 : 0;
    if (ImGui::CheckboxFlags("preCompDone", &pre_computation_is_done, 1)) {
      if (pre_computation_is_done == 1) {
        scene3d.Shapes[scene3d.ishape].preCompDone = 'y';
      } else {
        scene3d.Shapes[scene3d.ishape].preCompDone = 'n';
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Compute")) {
      scene3d.Shapes[scene3d.ishape].massProperties();
    }

    ImGui::Separator();
    ImGui::Text("Monte-Carlo No steps: %zu", scene3d.Shapes[scene3d.ishape].MCnstep);
    ImGui::Text("Volume: %g", scene3d.Shapes[scene3d.ishape].volume);
    ImGui::Text("I/m: %g %g %g", scene3d.Shapes[scene3d.ishape].inertia_mass.x,
                scene3d.Shapes[scene3d.ishape].inertia_mass.y, scene3d.Shapes[scene3d.ishape].inertia_mass.z);
  }

  ImGui::End();
}

// =======================
// Main functions
// =======================

// Traite les événements SDL
void ImGuiSDLApplication::handleEvents() {
  SDL_Event event;
  while (SDL_PollEvent(&event)) {
    ImGui_ImplSDL2_ProcessEvent(&event);

    checkKeyBindings(event);

    scene3d.handleMouseEvent(event, windowWidth, windowHeight);

    if (event.type == SDL_QUIT) {
      shouldQuit = true;
    } else if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_RESIZED) {
      scene3d.windowWidth = event.window.data1;
      scene3d.windowHeight = event.window.data2;

      scene3d.updateProjection();
    }
  }
}

void ImGuiSDLApplication::gradBackground() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_QUADS);
  glColor3f(clear_bottom_color[0], clear_bottom_color[1], clear_bottom_color[2]);
  glVertex2f(-1.0f, -1.0f);
  glVertex2f(1.0f, -1.0f);
  glColor3f(clear_top_color[0], clear_top_color[1], clear_top_color[2]);
  glVertex2f(1.0f, 1.0f);
  glVertex2f(-1.0f, 1.0f);
  glEnd();

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);

  scene3d.updateProjection();
  scene3d.updateCameraView();
  scene3d.applyProjectionAndView();
}

void ImGuiSDLApplication::renderFrame() {
  // Rendu ImGui en premier (il configure son état ET son viewport)
  ImGui::Render();

  int display_w, display_h;
  SDL_GetWindowSize(window, &display_w, &display_h);
  scene3d.windowWidth = display_w;
  scene3d.windowHeight = display_h;

  // Réinitialiser l'état OpenGL pour notre rendu 3D
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (scene3d.show_background) {
    gradBackground();
  } else {
    scene3d.updateProjection();
    scene3d.updateCameraView();
    scene3d.applyProjectionAndView();
  }
  scene3d.render();

  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

  // Échange des buffers
  SDL_GL_SwapWindow(window);
}

void ImGuiSDLApplication::cleanup() {

  saveConfigToSeparateFile();

  ImPlot::DestroyContext();

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  if (glContext) {
    SDL_GL_DeleteContext(glContext);
  }

  if (window) {
    SDL_DestroyWindow(window);
  }

  SDL_Quit();
}

// Boucle principale de l'application
void ImGuiSDLApplication::mainLoop() {

  while (!shouldQuit) {
    handleEvents();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();

    mainMenu();

    if (show_main_panel) {
      mainPanel();
    }

    renderFrame();
  }
}

// Point d'entrée principal de l'application
int ImGuiSDLApplication::run() {
  if (!init()) {
    return 1;
  }

  mainLoop();

  return 0;
}
