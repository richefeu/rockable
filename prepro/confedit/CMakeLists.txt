#file(GLOB confedit_SRCS "*.cpp")
#file(GLOB confedit_HDRS "*.hpp")

#FIND_PACKAGE(FLTK REQUIRED)

#include_directories(${cpptools_HDRS})
#add_executable(confedit confedit.cpp)
#target_link_libraries(confedit PUBLIC Rockable_core ${FLTK_LIBRARIES})
#install(TARGETS confedit)

cmake_minimum_required(VERSION 3.10)
project(confedit)

# Find FLTK package
find_package(FLTK REQUIRED)

if(FLTK_FOUND AND NOT TARGET fltk::fltk)
  add_library(fltk::fltk INTERFACE IMPORTED)
  set_target_properties(
    fltk::fltk
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${FLTK_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${FLTK_LIBRARIES}"
  )
endif()


#if(FLTK_FOUND)
#  message(STATUS "FLTK_INCLUDE_DIR: ${FLTK_INCLUDE_DIR}")
#  message(STATUS "FLTK_LIBRARIES: ${FLTK_LIBRARIES}")
#endif()

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(
    ${FLTK_INCLUDE_DIR}
    ~/toofus
    ../genesis/cpptools
)

# Add executable
add_executable(confedit confedit.cpp confedit.hpp)

# Link FLTK libraries
target_link_libraries(confedit PRIVATE fltk::fltk)

# Set compiler to c++-15 (if needed, adjust for your system)
set(CMAKE_CXX_COMPILER c++-15)



