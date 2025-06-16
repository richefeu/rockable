TOOFUSPATH = ~/toofus
ROCKABLEPATH = ../../../src

UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = g++-14
  CXXFLAGS = -O3 -Wall -std=c++11 -I $(TOOFUSPATH) -I $(ROCKABLEPATH)
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++0x -I $(TOOFUSPATH) -I $(ROCKABLEPATH)
endif

all: tess2shape

clean:
	rm -f *.o tess2shape

Shape.o: $(ROCKABLEPATH)/Core/Shape.cpp $(ROCKABLEPATH)/Core/Shape.hpp
	$(CXX) $(CXXFLAGS) -c $<

tess2shape: tess2shape.cpp main.cpp Shape.o
	$(CXX) $(CXXFLAGS) -c tess2shape.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
	$(CXX) $(CXXFLAGS) -o $@ main.cpp tess2shape.o Shape.o

