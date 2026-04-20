# Makefile pour l'application 3D avec ImGui

# Couleurs ANSI
RED    = \033[0;31m
GREEN  = \033[0;32m
YELLOW = \033[0;33m
BLUE   = \033[0;34m
NC     = \033[0m

# Détection automatique du système
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
    # Configuration pour macOS (Apple Silicon M1-M4 et Intel)
    CXX = g++-15
    CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -Wno-unknown-pragmas -Wno-missing-field-initializers
    
    INCLUDE_PATHS = \
	-I/opt/homebrew/include \
	-I/opt/homebrew/include/SDL2 \
	-I../mini-imgui-implot \
	-I../../../deps/toofus-src \
	-I `brew --prefix libomp`/include

    LDFLAGS = \
	-L/opt/homebrew/lib \
	-L `brew --prefix libomp`/lib
    
    LIBS = \
	-lSDL2 \
	-lSDL2main \
	-lGLEW \
	-framework OpenGL \
	-framework Cocoa \
	-framework IOKit \
	-framework CoreVideo \
	-lomp
else
    # Configuration pour systèmes Unix-like
    CXX = g++
    CXXFLAGS = -std=c++17 -Wall -Wextra
    
    INCLUDE_PATHS = \
	-I/usr/include/SDL2 \
	-I../mini-imgui-implot \
	-I../../../deps/toofus-src \
    
    LDFLAGS = \
	-L/usr/lib \
	-L/usr/lib/x86_64-linux-gnu
    
    LIBS = \
	-lSDL2 \
	-lSDL2main \
	-lGLEW \
	-lGL \
	-lX11 \
	-lXrandr \
	-lXinerama \
	-lXcursor \
	-lXi \
	-lXext
endif

# Fichiers Rockable (librairie déjà compilée)
ROCKABLE_INCLUDES = ../..
ROCKABLE_LIB = ../../../BUILD/libRockable_core.a

# Fichiers sources
SRCS = \
	shape.cpp \
	ImGuiSDLApplication.cpp \
	Scene3D.cpp \
	../mini-imgui-implot/imgui.cpp \
	../mini-imgui-implot/imgui_draw.cpp \
	../mini-imgui-implot/imgui_tables.cpp \
	../mini-imgui-implot/imgui_widgets.cpp \
	../mini-imgui-implot/imgui_demo.cpp \
	../mini-imgui-implot/imgui_impl_sdl2.cpp \
	../mini-imgui-implot/imgui_impl_opengl3.cpp \
	../mini-imgui-implot/implot.cpp \
	../mini-imgui-implot/implot_items.cpp

# Fichiers objets
OBJS = $(SRCS:.cpp=.o)

# Fichier cible
TARGET = shape

# Règle par défaut
all: $(TARGET)
	
# Règle de compilation pour chaque fichier .cpp
%.o: %.cpp
	@echo "$(BLUE)Compilation of $< ...$(NC)"
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATHS) -I$(ROCKABLE_INCLUDES) -c $< -o $@
	@echo "$(GREEN)OK : $@ compiled.$(NC)"

# Règle de compilation
$(TARGET): $(OBJS)
	@echo "$(YELLOW)Linking of objet files to build Dear My Shape...$(NC)"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(OBJS) $(ROCKABLE_LIB) -o $(TARGET)
	@echo "$(GREEN)OK : Dear My Shape is linked.$(NC)"
	
# Règle de nettoyage
clean:
	@echo "$(RED)Cleaning of objets files and app ...$(NC)"
	rm -f *.o $(TARGET)
	@echo "$(GREEN)OK : cleaning done.$(NC)"

clean+: clean
	@echo "$(RED)Clening of Dear ImGui ...$(NC)"
	rm -f ../mini-imgui-implot/*.o
	@echo "$(GREEN)OK : cleaning done.$(NC)"

.PHONY: all clean clean+
