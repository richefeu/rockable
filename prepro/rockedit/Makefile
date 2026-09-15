# rockedit builds with nothing but a C++17 compiler and the standard library,
# so that it can be compiled straight away on a cluster:  make && ./rockedit
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Isrc -I../common

TARGET = rockedit
HDRS   = src/terminal.hpp src/text_buffer.hpp ../common/rockable_lang.hpp

all: $(TARGET)

$(TARGET): src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

clean:
	rm -f $(TARGET)

.PHONY: all clean
