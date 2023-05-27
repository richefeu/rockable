# Paths
TOOFUSPATH = ~/toofus
ROCKABLEPATH = ../../src

UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = g++-13
  GLUTFLAGS = `pkg-config --cflags glut`
  GLLINK = `pkg-config --libs glut` -framework OpenGL
else
  CXX = g++
	GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

# The list of flags passed to the compiler
CXXFLAGS = -O3 -Wall -Wextra -Wshadow -pedantic -std=c++17 -I $(TOOFUSPATH) -I $(ROCKABLEPATH)/Core

# All source files (listed in SOURCES) will be compiled into an object file
# with the following command
%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean

all: shapeSurvey

clean:
	rm -f *.o shapeSurvey

Shape.o: $(ROCKABLEPATH)/Core/Shape.cpp $(ROCKABLEPATH)/Core/Shape.hpp
	$(CXX) $(CXXFLAGS) -c $<

# The application that visualize the conf files
shapeSurvey: shapeSurvey.cpp shapeSurvey.hpp Shape.o
	$(CXX) $(CXXFLAGS) $(GLUTFLAGS) -c $< -Wno-deprecated
	$(CXX) $(CXXFLAGS) -o $@ $@.o Shape.o $(GLLINK)
