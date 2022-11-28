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
OPTIONS = -DQUAT_ACC -DnFT_CORR -DnROT_MATRIX -DnBREAK_ONCE -DnCOMPONENTWISE_NUM_DAMPING -DnBIND_METHODS

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

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  # The compiler and linker to be used
  CXX = g++-12
	LINK = $(CXX)
	
  # The list of flags passed to the compiler
  CXXFLAGS := -fopenmp -O3 -Wall -Wextra -Wno-deprecated-declarations -std=c++17
  CXXFLAGS += -I $(TOOFUSPATH) -I $(SPDLOGPATH) $(OPTIONS) $(SVNDEV)
  CXXFLAGS += `pkg-config --cflags tclap`
	
  # Link flags for OpenGL and glut
  GLUTFLAGS = `pkg-config --cflags glut`
  GLUTLINK = `pkg-config --libs glut` -framework OpenGL
	
	FLTKFLAGS = `fltk-config --use-gl --cxxflags`
	FLTKLINK = `fltk-config --use-gl --ldstaticflags`
	
  # need to install glfw3 with brew: 'brew install glfw3'
  GLFWFLAGS = `pkg-config --cflags glfw3`
  GLFWLINK = -framework OpenGL `pkg-config --libs glfw3`

else
  # The compiler and linker to be used
  CXX = g++
  LINK = $(CXX)
	
  # The list of flags passed to the compiler
  CXXFLAGS := -fopenmp -O3 -std=c++17
  CXXFLAGS += -I $(TOOFUSPATH) $(OPTIONS) $(SVNDEV)
	
  # CXXFLAGS += `pkg-config --cflags tclap`
  # Link flags for OpenGL and glut
  GLUTLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
  
  GLFWFLAGS = `pkg-config --cflags glfw3`
  GLFWLINK = -lGLU -lGL `pkg-config --libs glfw3`
endif


# The list of source files needed by Rockable (the library)
SOURCES = Shape.cpp \
Particle.cpp \
Interaction.cpp \
Rockable.cpp \
DrivingSystem.cpp \
BreakableInterface.cpp \
ContactPartnership.cpp \
ForceLaw.cpp \
ForceLaw_Default.cpp \
ForceLaw_Avalanche.cpp \
ForceLaw_StickedLinks.cpp \
PreproCommand.cpp \
PreproCommand_stickVerticesInClusters.cpp \
PreproCommand_stickClusters.cpp \
PreproCommand_randomlyOrientedVelocities.cpp \
PreproCommand_randomlyOrientedVelocitiesClusters.cpp \
DataExtractor.cpp \
DataExtractor_TrackBody.cpp \
DataExtractor_TrackRockfall.cpp \
DataExtractor_ClusterAABB.cpp \
DataExtractor_dnStat.cpp \
DataExtractor_DuoBalance.cpp \
DataExtractor_MeanVelocity.cpp \
BodyForce.cpp \
BodyForce_PreferredDirection.cpp \
BodyForce_AttractingPoint.cpp \
BodyForce_ViscousFluid.cpp \
processingTool_probeSolidFraction.cpp \
processingTool_getInteractionGroups.cpp \
processingTool_getBrokenSubClusters.cpp \
processingTool_getMassRange.cpp \
processingTool_getClusters.cpp

IMGUI_SOURCES = $(IMGUI_DIR)/imgui.cpp \
$(IMGUI_DIR)/imgui_demo.cpp \
$(IMGUI_DIR)/imgui_draw.cpp \
$(IMGUI_DIR)/imgui_tables.cpp \
$(IMGUI_DIR)/imgui_widgets.cpp \
$(IMGUI_DIR)/backends/imgui_impl_glfw.cpp \
$(IMGUI_DIR)/backends/imgui_impl_opengl2.cpp

# The list of source files for post-processing (excepted postpro.cpp)
POSTPROSOURCES = PostProcessor.cpp \
PostProcessor_ClusterGranulo.cpp \
PostProcessor_ParticleStress.cpp

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:%.cpp=%.o)
POSTPROOBJECTS = $(POSTPROSOURCES:%.cpp=%.o)
IMGUI_OBJECTS = $(addsuffix .o, $(basename $(notdir $(IMGUI_SOURCES))))

# Each cpp file listed below corresponds to a header file
HEADERS = $(SOURCES:%.cpp=%.hpp)
POSTPROHEADERS = $(POSTPROSOURCES:%.cpp=%.hpp)

# All source files (listed in SOURCES or POSTPROSOURCES) will be compiled into an object file
# with the following command
%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean format

all: run see # postpro patatrac

clean:
	rm -f *~ *.o run see see2 see3 postpro patatrac edt

format:
	clang-format -i $(HEADERS) $(SOURCES) $(POSTPROHEADERS) $(POSTPROSOURCES) run.cpp see.cpp patatrac.cpp postpro.cpp

# The application that runs a simulation
run: run.cpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)
	
# Scripting with chaiscript (TEST!!!)
script: script.cpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -I $(CHAISCRIPTPATH)  -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)
	
# The application that edits and visualize the conf files (it uses FLTK)
edt: edt.cpp edt.hpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(FLTKFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(FLTKLINK) 

# An application for post-processing the conf-files
postpro: postpro.cpp $(HEADERS) $(SOURCES) $(OBJECTS) $(POSTPROHEADERS) $(POSTPROSOURCES) $(POSTPROOBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(POSTPROOBJECTS)

# An application to use Rockable as a "trajectography-analysis" tool
patatrac: patatrac.cpp patatrac.hpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)

# The application that visualizes the conf files (it uses GLUT)
see: see.cpp see.hpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GLUTFLAGS) $(PNGFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(GLUTLINK) $(PNGLINK)
	
# The application that visualizes the conf files (it uses GLFW)
see2: see2.cpp see2.hpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GLFWFLAGS) $(PNGFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(GLFWLINK) $(PNGLINK)
	
# The application that visualizes the conf files (it uses GLFW and Dear ImGui)
$(IMGUI_OBJECTS): $(IMGUI_SOURCES)
	$(CXX) $(CXXFLAGS) -I $(IMGUI_DIR) -I $(IMGUI_DIR)/backends -Wno-unused-function $(GLFWFLAGS) -c $(IMGUI_SOURCES)
	
see3: see3.cpp see3.hpp $(HEADERS) $(OBJECTS) $(IMGUI_OBJECTS)
	$(CXX) $(CXXFLAGS) -I $(IMGUI_DIR) -I $(IMGUI_DIR)/backends $(GLFWFLAGS) $(PNGFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(IMGUI_OBJECTS) $(GLFWLINK) $(PNGLINK)

#######################################################################################################################
