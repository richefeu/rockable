cmake_minimum_required(VERSION 3.26)

if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")
  # Try to find the latest GNU g++ installed via Homebrew
  execute_process(
    COMMAND
      bash -c
      "ls /usr/local/bin/g++-* /opt/homebrew/bin/g++-* 2>/dev/null | grep -Eo 'g\\+\\+-[0-9]+' | sort -V | tail -1"
    OUTPUT_VARIABLE LATEST_GNU_GPP
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(LATEST_GNU_GPP)
    set(CMAKE_OSX_COMPILER_ID "GNU")
    set(CMAKE_CXX_COMPILER "${LATEST_GNU_GPP}")
    message(STATUS "Using GNU g++: ${LATEST_GNU_GPP}")
  else()
    message(
      WARNING
        "No GNU g++ found in /usr/local/bin/ or /opt/homebrew/bin/. Falling back to default compiler."
    )
  endif()
endif()

project(Rockable LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# All compiled applications will be in the BUILD folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ##############################################################################
# OPTIONS

option(ROCKABLE_USE_FT_CORR "Add objectivity correction to tangent forces" OFF)
option(ROCKABLE_USE_TESTING "Add tests" OFF)

option(ROCKABLE_ENABLE_PROFILING "Enable the time profiling" OFF)
option(ROCKABLE_ENABLE_BOUNDARY "Enable the special boundaries like Ball or Cylinder" OFF)
option(ROCKABLE_ENABLE_SOFT_PARTICLES "Enable strainning of particles" OFF)
option(ROCKABLE_ENABLE_PERIODIC "Enable full periodic boundary conditions" OFF)

option(ROCKABLE_COMPILE_POSTPRO "Compile the application to run post-processing commands" ON)
option(ROCKABLE_COMPILE_SEE "Compile the glut application to visualize the conf-files" ON)
option(ROCKABLE_COMPILE_SEER "Compile Dear My Seer (and Dear My Shape)" OFF)
option(ROCKABLE_COMPILE_PREPRO "Compile the preprocessing tools" ON)
option(ROCKABLE_COMPILE_CONF2VTK "Convert conf files to VTK" OFF)

# ##############################################################################
# compilation flags

if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Build type is set to Release")
endif()

if(ROCKABLE_USE_FT_CORR)
  add_definitions(-DFT_CORR)
endif()

if(ROCKABLE_ENABLE_PROFILING)
  add_compile_options(-DROCKABLE_ENABLE_PROFILING)
endif()

if(ROCKABLE_ENABLE_BOUNDARY)
  add_compile_options(-DROCKABLE_ENABLE_BOUNDARY)
endif()

if(ROCKABLE_ENABLE_PERIODIC)
  add_compile_options(-DROCKABLE_ENABLE_PERIODIC)
endif()

# get the git tag of Rockable
execute_process(
  COMMAND git describe --abbrev=0 --tags --always
  OUTPUT_VARIABLE ROCKABLE_GIT_TAG
  OUTPUT_STRIP_TRAILING_WHITESPACE)
add_definitions(-DROCKABLE_GIT_TAG="${ROCKABLE_GIT_TAG}")
message(STATUS "ROCKABLE_GIT_TAG = " ${ROCKABLE_GIT_TAG})

# get the date of last git commit
execute_process(
  COMMAND git log -1 --format=%cd
  OUTPUT_VARIABLE ROCKABLE_LAST_COMMIT_DATE
  OUTPUT_STRIP_TRAILING_WHITESPACE)
add_definitions(-DROCKABLE_LAST_COMMIT_DATE="${ROCKABLE_LAST_COMMIT_DATE}")
message(STATUS "ROCKABLE_LAST_COMMIT_DATE = " ${ROCKABLE_LAST_COMMIT_DATE})

# the default mode is release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

include(FetchContent)
set(FETCHCONTENT_BASE_DIR ${CMAKE_SOURCE_DIR}/deps)

# ##############################################################################
# Fetching ToOfUs
FetchContent_Declare(
  toofus
  GIT_REPOSITORY https://github.com/richefeu/toofus.git
  GIT_TAG main)
FetchContent_GetProperties(toofus)
if(NOT toofus_POPULATED)
  message(STATUS "Fetching toofus")
  FetchContent_MakeAvailable(toofus)
endif()

# ##############################################################################
# MATools
if(ROCKABLE_ENABLE_PROFILING)
  FetchContent_Declare(
    matools
    GIT_REPOSITORY https://github.com/rprat-pro/MATools
    GIT_TAG v1.1)

  FetchContent_GetProperties(matools)

  if(NOT matools_POPULATED)
    message(STATUS "Fetching MATools")
    FetchContent_MakeAvailable(matools)
    set(MATOOLS_MPI
        OFF
        CACHE BOOL "Use some expat option")
    add_subdirectory(${FETCHCONTENT_BASE_DIR}/matools-src)
  endif()
endif()

# ##############################################################################
# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP CXX has been found")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# ##############################################################################
# OpenGL and glut
find_package(OpenGL REQUIRED)
if(OPENGL_FOUND)
  message(STATUS "OPENGL_INCLUDE_DIR = " ${OPENGL_INCLUDE_DIR})
  message(STATUS "OPENGL_LIBRARIES = " ${OPENGL_LIBRARIES})
  message(STATUS "OPENGL_glu_LIBRARY = ${OPENGL_glu_LIBRARY}")
endif()

find_package(GLUT REQUIRED)
if(GLUT_FOUND)
  message(STATUS "GLUT_INCLUDE_DIR = " ${GLUT_INCLUDE_DIR})
  message(STATUS "GLUT_LIBRARIES = " ${GLUT_LIBRARIES})
endif()

# ##############################################################################
# libRockable_core
set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

include_directories(${CMAKE_SOURCE_DIR} ${toofus_SOURCE_DIR} ${fmt_SOURCE_DIR}/include)

file(
  GLOB
  core_object_files
  ${CMAKE_SOURCE_DIR}/BodyForces/*.cpp
  ${CMAKE_SOURCE_DIR}/Core/*.cpp
  ${CMAKE_SOURCE_DIR}/DataExtractors/*.cpp
  ${CMAKE_SOURCE_DIR}/ForceLaws/*.cpp
  ${CMAKE_SOURCE_DIR}/PostProcessors/*.cpp
  ${CMAKE_SOURCE_DIR}/PreproCommands/*.cpp
  ${CMAKE_SOURCE_DIR}/ProcessingTools/*.cpp)

if(ROCKABLE_ENABLE_BOUNDARY)
  file(GLOB boundaries_files ${CMAKE_SOURCE_DIR} /Boundaries/*.cpp)
  list(APPEND core_object_files ${boundaries_files})
endif() # end of ROCKABLE_ENABLE_BOUNDARY

add_library(Rockable_core STATIC ${core_object_files})
if(ROCKABLE_ENABLE_PROFILING)
  target_link_libraries(Rockable_core MATools)
endif()

# ##############################################################################
# rockable executable (to run simulations from the console)
add_executable(rockable ${CMAKE_SOURCE_DIR}/Apps/run.cpp)
target_link_libraries(rockable PUBLIC Rockable_core)
install(TARGETS rockable)

# ##############################################################################
# postpro
if(ROCKABLE_COMPILE_POSTPRO)
  add_executable(postpro ${CMAKE_SOURCE_DIR}/Apps/postpro.cpp)
  target_link_libraries(postpro PUBLIC Rockable_core)
  install(TARGETS postpro)
endif()

# ##############################################################################
# visualiser see (with glut or freeglut)
if(ROCKABLE_COMPILE_SEE)

  add_executable(see ${CMAKE_SOURCE_DIR}/Apps/see.cpp)

  target_include_directories(
    see
    PRIVATE ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} /usr/local/include)

  target_link_libraries(see PRIVATE OpenGL::GL OpenGL::GLU ${GLUT_LIBRARIES}
                                    Rockable_core)

  install(TARGETS see)

endif() # end of ROCKABLE_COMPILE_SEE

# ##############################################################################
# Dear My Seer
if(ROCKABLE_COMPILE_SEER)
  
  find_package(GLEW REQUIRED)
  message(STATUS "GLEW_INCLUDE_DIRS = " ${GLEW_INCLUDE_DIRS})
  message(STATUS "GLEW_LIBRARIES = " ${GLEW_LIBRARIES})
  
  find_package(SDL2 REQUIRED)
  message(STATUS "SDL2_INCLUDE_DIRS = " ${SDL2_INCLUDE_DIRS})
  message(STATUS "SDL2_LIBRARIES = " ${SDL2_LIBRARIES})
  
  #find_package(OpenGL REQUIRED)
  #message(STATUS "OPENGL_INCLUDE_DIR = " ${OPENGL_INCLUDE_DIR})
  #message(STATUS "OPENGL_LIBRARIES = " ${OPENGL_LIBRARIES})
 
  set(SEER_FILES
      ${CMAKE_SOURCE_DIR}/Apps/DearMySeer/seer.cpp
      ${CMAKE_SOURCE_DIR}/Apps/DearMySeer/ImGuiSDLApplication.cpp
      ${CMAKE_SOURCE_DIR}/Apps/DearMySeer/Scene3D.cpp
      ${CMAKE_SOURCE_DIR}/Apps/DearMySeer/ParticleColoring.cpp
      ${CMAKE_SOURCE_DIR}/Apps/DearMySeer/Slicing.cpp      
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/imgui.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/imgui_demo.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/imgui_draw.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/imgui_tables.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/imgui_widgets.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/imgui_impl_sdl2.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/imgui_impl_opengl3.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/implot.cpp
      ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot/implot_items.cpp)

  add_executable(seer ${SEER_FILES})

  target_include_directories(seer
      PUBLIC ${OPENGL_INCLUDE_DIR} ${SDL2_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS}
             ${CMAKE_SOURCE_DIR}/Apps/DearMySeer ${CMAKE_SOURCE_DIR}/Apps/mini-imgui-implot)

  target_link_libraries(seer PUBLIC OpenGL::GL SDL2::SDL2 GLEW::GLEW Rockable_core)

  install(TARGETS seer)
  
endif() # end of ROCKABLE_COMPILE_SEER

# ##############################################################################
# conf2vtk
if(ROCKABLE_COMPILE_CONF2VTK)
  add_executable(conftovtk ${CMAKE_SOURCE_DIR}/Apps/conftovtk.cpp)
  target_link_libraries(conftovtk Rockable_core)
  install(TARGETS conftovtk)
endif()

# ##############################################################################
# Enable non regression tests
if(ROCKABLE_USE_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

# ##############################################################################
# Pre-processing tools
if(ROCKABLE_COMPILE_PREPRO)
  add_subdirectory(prepro)
endif()
