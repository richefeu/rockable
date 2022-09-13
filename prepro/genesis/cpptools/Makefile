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
	
