UNAME_S := $(shell uname -s)

COMMONPATH = /usr/local/include/toofus
CPPTOOLSPATH = ../cpptools

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = clang++
  CXXFLAGS = -O3 -Wall -std=c++11 -I $(COMMONPATH) -I $(CPPTOOLSPATH)
  GLFLAGS = `pkg-config --cflags --libs glut` -framework OpenGL 
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++0x -I $(COMMONPATH) -I $(CPPTOOLSPATH)
  GLFLAGS = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

all: pack seePack

clean:
	rm -f *.o pack seePack

pack: TubePacker.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $<
	
seePack: seePack.cpp
	$(CXX) $(CXXFLAGS) -Wno-deprecated $(GLFLAGS) -o $@ $<
