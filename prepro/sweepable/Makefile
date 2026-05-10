CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

TARGET = sweepable
SOURCES = sweepable.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -I. $^ -o $@

clean:
	rm -f $(TARGET)

.PHONY: all clean
