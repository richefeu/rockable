#######################################################################################################################
# get the svn revision number
SVNDEV := -D'SVN_REV="$(shell svnversion -n .)"'

# An option can be disable by adding 'n' before
OPTIONS = -DQUAT_ACC -DnFT_CORR -DnROT_MATRIX -DnBREAK_ONCE

# Paths
COMMONPATH = ../common

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  # The compiler to be used
  CXX = g++-8
  # The list of flags passed to the compiler
  CXXFLAGS = -fopenmp -O3 -Wall -Wextra -Wno-deprecated-declarations -std=c++11 -I $(COMMONPATH) $(OPTIONS) $(SVNDEV)
  # Link flags for OpenGL and glut
  GLLINK = -framework OpenGL -framework GLUT
else
  # The compiler to be used
  CXX = g++
  # The list of flags passed to the compiler
  CXXFLAGS = -fopenmp -O3 -std=c++0x -I $(COMMONPATH) $(OPTIONS) $(SVNDEV)
  # Link flags for OpenGL and glut
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

# If you don't want to use libpng in 'see', comment this two lines
# In this case the screenshots will be TGA (not compressed)
OPTIONS += -DPNG_SCREENSHOTS
GLLINK += -lpng

# The list of source files needed by Rockable
SOURCES = Shape.cpp \
Particle.cpp \
Interaction.cpp \
Rockable.cpp \
DrivingSystem.cpp \
BreakableInterface.cpp \
ContactPartnership.cpp \
DataExtractor.cpp \
DataExtractor_TrackBody.cpp \
DataExtractor_ClusterAABB.cpp \
DataExtractor_dnStat.cpp \
DataExtractor_DuoBalance.cpp \
BodyForce.cpp \
BodyForce_PreferredDirection.cpp \
BodyForce_AttractingPoint.cpp

# The list of source files for post-processing
POSTPROSOURCES = PostProcessor.cpp \
PostProcessor_ClusterGranulo.cpp \
PostProcessor_ParticleStress.cpp

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:%.cpp=%.o)
POSTPROOBJECTS = $(POSTPROSOURCES:%.cpp=%.o)

# Each cpp file listed below corresponds to a header file
HEADERS = $(SOURCES:%.cpp=%.hpp)
POSTPROHEADERS = $(POSTPROSOURCES:%.cpp=%.hpp)

# All source files (listed in SOURCES or POSTPROSOURCES) will be compiled into an object file
# with the following command
%.o:%.cpp Makefile
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean format

all: run see postpro patatrac

clean:
	rm -f *~ *.o run see postpro patatrac

format:
	clang-format -i $(HEADERS) $(SOURCES) $(POSTPROHEADERS) $(POSTPROSOURCES) run.cpp see.cpp patatrac.cpp postpro.cpp

# The application that runs a simulation
run: run.cpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)

# The application for post-processing the conf-files
postpro: postpro.cpp $(HEADERS) $(SOURCES) $(OBJECTS) $(POSTPROHEADERS) $(POSTPROSOURCES) $(POSTPROOBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(POSTPROOBJECTS)

# An application to use Rockable as a "trajectography-analysis" tool
patatrac: patatrac.cpp patatrac.hpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)

# The application that visualizes the conf files
see: see.cpp see.hpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(GLLINK)

#######################################################################################################################
