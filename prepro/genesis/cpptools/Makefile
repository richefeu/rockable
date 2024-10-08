# This c++ application is so simple that it can easily be compiled
# thanks to this simple Makefile. When using this Makefile, make sure
# that the toofus folder is in your home directory 
# (update its content if required)

COMMONPATH = ~/toofus
CPPTOOLSPATH = .

UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = clang++
  CXXFLAGS = -O3 -Wall -std=c++11 -I $(COMMONPATH) -I $(CPPTOOLSPATH)
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++0x -I $(COMMONPATH) -I $(CPPTOOLSPATH)
endif

all: generator

clean:
	rm -f *.o generator

generator: generator.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
	
