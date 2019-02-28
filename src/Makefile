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
