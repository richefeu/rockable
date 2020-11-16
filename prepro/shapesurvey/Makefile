UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = clang++
  #GLFLAGS = -framework OpenGL -framework GLUT
  GLUTFLAGS = `pkg-config --cflags glut`
  GLLINK = `pkg-config --libs glut` -framework OpenGL
else
  CXX = g++
  #GLFLAGS = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
	GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

# The list of flags passed to the compiler
CXXFLAGS = -O3 -Wall -std=c++11 -I ../../common

# All source files (listed in SOURCES) will be compiled into an object file
# with the following command
%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean

all: shapeSurvey

clean:
	rm -f *.o shapeSurvey

Shape.o: ../Rockable/Shape.cpp ../Rockable/Shape.hpp
	$(CXX) $(CXXFLAGS) -c $<

# The application that visualize the conf files
shapeSurvey: shapeSurvey.cpp shapeSurvey.hpp Shape.o
	$(CXX) $(CXXFLAGS) $(GLUTFLAGS) -c $< -Wno-deprecated
	$(CXX) $(CXXFLAGS) -o $@ $@.o Shape.o $(GLLINK)
