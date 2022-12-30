#  Copyright or © or Copr. Rockable
#  
#  vincent.richefeu@3sr-grenoble.fr
#  
#  This software is a computer program whose purpose is 
#    (i)  to hold sphero-polyhedral shapes,
#    (ii) to manage breakable interfaces. 
#  It is developed for an ACADEMIC USAGE
#  
#  This software is governed by the CeCILL-B license under French law and
#  abiding by the rules of distribution of free software.  You can  use, 
#  modify and/ or redistribute the software under the terms of the CeCILL-B
#  license as circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info". 
#  
#  As a counterpart to the access to the source code and  rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty  and the software's author,  the holder of the
#  economic rights,  and the successive licensors  have only  limited
#  liability. 
#  
#  In this respect, the user's attention is drawn to the risks associated
#  with loading,  using,  modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean  that it is complicated to manipulate,  and  that  also
#  therefore means  that it is reserved for developers  and  experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or 
#  data to be ensured and,  more generally, to use and operate it in the 
#  same conditions as regards security. 
#  
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL-B license and that you accept its terms.
#######################################################################################################################

# An option can be disable by adding 'n' before
OPTIONS = -DnFT_CORR -DnCOMPONENTWISE_NUM_DAMPING -DENABLE_PROFILING

# If you don't want to use libpng in 'see', comment this two lines
# In this case the screenshots will be TGA (not compressed)
#PNGFLAGS = -DPNG_H `pkg-config --cflags libpng`
#PNGLINK = `pkg-config --libs libpng`

# For the compilation, toofus is required. 
# It can be installed with the following command: 
#    cd ~ && rm -rf toofus && git clone https://github.com/richefeu/toofus.git && cd -

# Paths
TOOFUSPATH = ~/toofus
IMGUI_DIR = ./imgui
SPDLOGPATH = /usr/local/include
CHAISCRIPTPATH = ~/ChaiScript-6.1.0/include
LIB = libRockable

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)

  # The compiler and linker to be used
  CXX = g++-12
  LINK = $(CXX)
	
  # The list of flags passed to the compiler
  CXXFLAGS = -fopenmp -O3 -Wall -Wextra -Wshadow -pedantic -std=c++17
  CXXFLAGS += -I $(TOOFUSPATH) -I $(SPDLOGPATH) -I . $(OPTIONS)
  CXXFLAGS += `pkg-config --cflags tclap`

  # need to install nlohmann/json: 'brew install nlohmann-json'
  JSONFLAGS = `pkg-config --cflags nlohmann_json`
	
  OPENMP_LIB = -l gomp
	
  # Link flags for OpenGL and glut
  GLUTFLAGS = `pkg-config --cflags glut`
  GLUTLINK = `pkg-config --libs glut` -framework OpenGL
	
  FLTKFLAGS = `fltk-config --use-gl --cxxflags`
  FLTKLINK = `fltk-config --use-gl --ldstaticflags`
	
  # need to install glfw3 with brew: 'brew install glfw3' or 'brew install glfw' (it depends...)
  GLFWFLAGS = `pkg-config --cflags glfw3`
  GLFWLINK = -framework OpenGL `pkg-config --libs glfw3`

else

  # The compiler and linker to be used
  CXX = g++
  LINK = $(CXX)
	
  # The list of flags passed to the compiler
  CXXFLAGS = -fopenmp -O3 -std=c++17
  CXXFLAGS += -I $(TOOFUSPATH) $(OPTIONS)
	
  # CXXFLAGS += `pkg-config --cflags tclap`
  # Link flags for OpenGL and glut
  GLUTLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
  
  GLFWFLAGS = `pkg-config --cflags glfw3`
  GLFWLINK = -lGLU -lGL `pkg-config --libs glfw3`
	
endif

#######################################################################################################################

# The list of source files needed by Rockable (the library)
SOURCES = \
$(wildcard Core/*.cpp) \
$(wildcard BodyForces/*.cpp) \
$(wildcard DataExtractors/*.cpp) \
$(wildcard ForceLaws/*.cpp) \
$(wildcard PreproCommands/*.cpp) \
$(wildcard ProcessingTools/*.cpp)

IMGUI_SOURCES = $(IMGUI_DIR)/imgui.cpp \
$(IMGUI_DIR)/imgui_demo.cpp \
$(IMGUI_DIR)/imgui_draw.cpp \
$(IMGUI_DIR)/imgui_tables.cpp \
$(IMGUI_DIR)/imgui_widgets.cpp \
$(IMGUI_DIR)/backends/imgui_impl_glfw.cpp \
$(IMGUI_DIR)/backends/imgui_impl_opengl2.cpp

# The list of source files for post-processing (excepted postpro.cpp)
POSTPROSOURCES = $(wildcard PostProcessors/*.cpp) 

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:.cpp=.o)
POSTPROOBJECTS = $(POSTPROSOURCES:.cpp=.o)
IMGUI_OBJECTS = $(addsuffix .o, $(basename $(IMGUI_SOURCES)))

# Each cpp file listed below corresponds to a header file
HEADERS = $(SOURCES:.cpp=.hpp)
POSTPROHEADERS = $(POSTPROSOURCES:.cpp=.hpp)

#######################################################################################################################

$(IMGUI_DIR)/%.o: $(IMGUI_DIR)/%.cpp
	@echo "\033[0;32m-> COMPILING" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -I $(IMGUI_DIR) -I $(IMGUI_DIR)/backends -Wno-unused-function $(GLFWFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	
# All source files (listed in SOURCES or POSTPROSOURCES) will be compiled into an object file
# with the following command
%.o:%.cpp
	@echo "\033[0;32m-> COMPILING" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
.PHONY: all clean clean+

all: lib rockable see

lib: lib/$(LIB).a
	
lib/$(LIB).a: $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p lib
	ar rs $@ $^

# The application that runs a simulation
rockable: bin/rockable
bin/rockable: Apps/run.cpp $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	$(LINK) $(OPENMP_LIB) -o $@ Apps/run.o lib/$(LIB).a
clean_rockable:
	rm -rf Apps/run.o bin/rockable

# The application that converts conf-files into vtk-files
conftovtk: bin/conftovtk
bin/conftovtk: Apps/conftovtk.cpp $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	$(CXX) $(OPENMP_LIB) -o $@ Apps/conftovtk.o lib/$(LIB).a
	
# Scripting with chaiscript
script: bin/script
bin/script: Apps/script.cpp $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -I $(CHAISCRIPTPATH) -DCHAISCRIPT_NO_THREADS -c $< -o $(patsubst %.cpp, %.o, $<)
	$(LINK) $(OPENMP_LIB) -o $@ Apps/script.o lib/$(LIB).a -ldl
clean_script:
	rm -rf Apps/script.o bin/script

# An application for post-processing the conf-files
postpro: bin/postpro
bin/postpro: Apps/postpro.cpp $(OBJECTS) $(POSTPROOBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	$(LINK) $(OPENMP_LIB) -o $@ Apps/postpro.o lib/$(LIB).a $(POSTPROOBJECTS)
clean_postpro:
	rm -rf Apps/postpro.o bin/postpro

# An application to use Rockable as a "trajectography-analysis" tool
patatrac: bin/patatrac
bin/patatrac: Apps/patatrac.cpp $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	$(LINK) $(CXXFLAGS) -o $@ Apps/patatrac.o lib/$(LIB).a
clean_patatrac:
	rm -rf Apps/patatrac.o bin/patatrac

# The application that visualizes the conf files (it uses GLUT)
see: bin/see
bin/see: Apps/see.cpp $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(GLUTFLAGS) $(PNGFLAGS) $(JSONFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	$(LINK) $(OPENMP_LIB) $(GLUTLINK) $(PNGLINK) -o $@ Apps/see.o lib/$(LIB).a
clean_see:
	rm -rf Apps/see.o bin/see
	
# The application that visualizes the conf files (it uses GLFW)
see2: bin/see2
bin/see2: Apps/see2.cpp $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(GLFWFLAGS) $(PNGFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	$(LINK) $(OPENMP_LIB) $(GLFWLINK) $(PNGLINK) -o $@ Apps/see2.o lib/$(LIB).a
clean_see2:
	rm -rf Apps/see2.o bin/see2
	
# The application that visualizes the conf files (it uses GLFW and Dear ImGui)
see3: bin/see3
bin/see3: Apps/see3.cpp $(OBJECTS) $(IMGUI_OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -I $(IMGUI_DIR) -I $(IMGUI_DIR)/backends $(GLFWFLAGS) $(PNGFLAGS) -c $< -o $(patsubst %.cpp, %.o, $<)
	$(LINK) $(OPENMP_LIB) $(GLFWLINK) $(PNGLINK) -o $@ Apps/see3.o lib/$(LIB).a $(IMGUI_OBJECTS)
clean_see3:
	rm -rf Apps/see3.o bin/see3

clean_objects:
	rm -f $(OBJECTS)
clean_postpro_objects:
	rm -f $(POSTPROOBJECTS)
clean_imgui_objects:
	rm -f $(IMGUI_OBJECTS)

clean: clean_objects clean_postpro_objects clean_imgui_objects
	
clean+: clean
	rm -rf lib
	rm -rf bin
	
#######################################################################################################################
