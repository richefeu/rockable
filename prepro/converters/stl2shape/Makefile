# Paths
COMMONPATH = ~/toofus

UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = clang++
  CXXFLAGS = -O3 -Wall -std=c++11 -I $(COMMONPATH) `pkg-config --cflags tclap`
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++0x -I $(COMMONPATH)
endif

all: stl2shape

clean:
	rm -f *.o stl2shape

stl2shape: stl2shape.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $<
