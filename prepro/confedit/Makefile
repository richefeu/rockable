#CXX      = $(shell fltk-config --cxx)
CXX      = c++-15
LINK     = $(CXX)
CXXFLAGS = $(shell fltk-config --cxxflags ) -Wall -std=c++17 -I ~/toofus -I ../genesis/cpptools
LDFLAGS  = $(shell fltk-config --ldflags )
LDSTATIC = $(shell fltk-config --ldstaticflags )
 
TARGET = confedit
OBJS = 
 
.SUFFIXES: .o .cpp
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
 
all: $(TARGET)

$(TARGET): $(TARGET).cpp $(TARGET).hpp
	$(CXX) $(CXXFLAGS) -c $(TARGET).cpp
	$(LINK) $(CXXFLAGS) -o $(TARGET) $(TARGET).o $(OBJS) $(LDSTATIC)
 
clean:
	rm -f $(TARGET) $(TARGET).o $(OBJS)
