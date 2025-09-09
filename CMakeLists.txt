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

option(ROCKABLE_USE_FT_CORR "add objectivity correction to tangent forces" OFF)
option(ROCKABLE_USE_TESTING "add tests" OFF)

option(ROCKABLE_ENABLE_PROFILING "enable the time profiling" OFF)
option(ROCKABLE_ENABLE_BOUNDARY
       "enable the special boundaries like Ball or Cylinder" OFF)
option(ROCKABLE_ENABLE_SOFT_PARTICLES "enable strainning of particles" OFF)
option(ROCKABLE_ENABLE_PERIODIC "enable full periodic boundary conditions" OFF)

option(ROCKABLE_COMPILE_POSTPRO
       "Compile the application to run post-processing commands" ON)
option(ROCKABLE_COMPILE_SEE
       "Compile the glut application to visualize the conf-files" ON)
option(ROCKABLE_COMPILE_SEE2
       "Compile the glfw application to visualize the conf-files" OFF)
option(ROCKABLE_COMPILE_SEE3
       "Compile the application to edit graphically the input files" OFF)
option(ROCKABLE_COMPILE_SEE4
       "Compile the application to edit graphically the input files" OFF)
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
# Fetching fmt
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
  GIT_PROGRESS TRUE
  USES_TERMINAL_DOWNLOAD TRUE)
FetchContent_GetProperties(fmt)
if(NOT fmt_POPULATED)
  message(STATUS "Fetching fmt")
  FetchContent_MakeAvailable(fmt)
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
# tclap
FetchContent_Declare(
  tclap
  GIT_REPOSITORY https://github.com/xguerin/tclap.git
  GIT_TAG v1.2.5)

FetchContent_GetProperties(tclap)

if(NOT tclap_POPULATED)
  message(STATUS "Fetching tclap")
  FetchContent_MakeAvailable(tclap)
endif()

# ##############################################################################
# nlohmann-json
FetchContent_Declare(
  nlohmann-json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.12.0)

FetchContent_GetProperties(nlohmann-json)

if(NOT nlohmann-json_POPULATED)
  message(STATUS "Fetching nlohmann-json")
  FetchContent_MakeAvailable(nlohmann-json)
endif()

# ##############################################################################
# imgui

if(ROCKABLE_COMPILE_SEE3 OR ROCKABLE_COMPILE_SEE4)

  FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG master)

  FetchContent_GetProperties(imgui)

  if(NOT imgui_POPULATED)
    message(STATUS "Fetching imgui")
    FetchContent_MakeAvailable(imgui)
  endif()

endif() # end of ROCKABLE_COMPILE_SEE3 or 4

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
if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")

  # here, we do not uses find_package(glut) so that the version installed by
  # apple is NOT used
  if(NOT PkgConfig_FOUND)
    find_package(PkgConfig REQUIRED)
  endif()

  pkg_check_modules(GL REQUIRED gl)
  pkg_check_modules(GLU REQUIRED glu)
  pkg_check_modules(GLUT REQUIRED glut>=3.4.0)

  if(GL_FOUND)
    message(STATUS "GL as been found with PkgConfig")
    message(STATUS "GL_LIBRARY_DIRS = " ${GL_LIBRARY_DIRS})
    message(STATUS "GL_INCLUDE_DIRS = " ${GL_INCLUDE_DIRS})
    message(STATUS "GL_LIBRARIES = " ${GL_LIBRARIES})
  endif()

  if(GLU_FOUND)
    message(STATUS "GLU as been found with PkgConfig")
    message(STATUS "GLU_LIBRARY_DIRS = " ${GLU_LIBRARY_DIRS})
    message(STATUS "GLU_INCLUDE_DIRS = " ${GLU_INCLUDE_DIRS})
    message(STATUS "GLU_LIBRARIES = " ${GLU_LIBRARIES})
  endif()

  if(GLUT_FOUND)
    message(STATUS "GLUT as been found with PkgConfig")
    message(STATUS "GLUT_LIBRARY_DIRS = " ${GLUT_LIBRARY_DIRS})
    message(STATUS "GLUT_INCLUDE_DIRS = " ${GLUT_INCLUDE_DIRS})
    message(STATUS "GLUT_LIBRARIES = " ${GLUT_LIBRARIES})
  endif()

else()

  find_package(OpenGL REQUIRED)
  if(OPENGL_FOUND)
    message(STATUS "OPENGL_INCLUDE_DIR = " ${OPENGL_INCLUDE_DIR})
    message(STATUS "OPENGL_LIBRARIES = " ${OPENGL_LIBRARIES})
  endif()

  find_package(GLUT REQUIRED)
  if(GLUT_FOUND)
    message(STATUS "GLUT_INCLUDE_DIRS = " ${GLUT_INCLUDE_DIRS})
    message(STATUS "GLUT_LIBRARIES = " ${GLUT_LIBRARIES})
  endif()

endif()

# ##############################################################################
# glfw
if(ROCKABLE_COMPILE_SEE2 OR ROCKABLE_COMPILE_SEE3)

  if(NOT PkgConfig_FOUND)
    find_package(PkgConfig REQUIRED)
  endif()

  pkg_check_modules(GLFW REQUIRED glfw3)

  if(GLFW_FOUND)
    message(STATUS "GLFW as been found with PkgConfig")
    message(STATUS "GLFW_INCLUDE_DIRS = " ${GLFW_INCLUDE_DIRS})
    message(STATUS "GLFW_LIBRARY_DIRS = " ${GLFW_LIBRARY_DIRS})
    message(STATUS "GLFW_LIBRARIES = " ${GLFW_LIBRARIES})
  endif()

endif()

# ##############################################################################
# sdl2
if(ROCKABLE_COMPILE_SEE4)

  if(NOT PkgConfig_FOUND)
    find_package(PkgConfig REQUIRED)
  endif()

  # Check for SDL2 instead of GLFW
  pkg_check_modules(SDL2 REQUIRED sdl2)

  if(SDL2_FOUND)
    message(STATUS "SDL2 has been found with PkgConfig")
    message(STATUS "SDL2_LIBRARY_DIRS = " ${SDL2_LIBRARY_DIRS})
    message(STATUS "SDL2_INCLUDE_DIRS = " ${SDL2_INCLUDE_DIRS})
    message(STATUS "SDL2_LIBRARIES = " ${SDL2_LIBRARIES})
  endif()

endif()

# ##############################################################################
# libRockable_core
set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

include_directories(${CMAKE_SOURCE_DIR} ${toofus_SOURCE_DIR}
                    ${tclap_SOURCE_DIR}/include ${fmt_SOURCE_DIR}/include)

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
# visualiser see (it uses freeglut)
if(ROCKABLE_COMPILE_SEE)

  if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux") #=========

    find_package(PNG REQUIRED)
    if(PNG_FOUND)
      message(STATUS "libpng found: ${PNG_INCLUDE_DIRS}")
      message(STATUS "PNG_INCLUDE_DIRS = " ${PNG_INCLUDE_DIRS})
      message(STATUS "PNG_LIBRARY_DIRS = " ${PNG_LIBRARY_DIRS})
      message(STATUS "PNG_LIBRARY = " ${PNG_LIBRARY})
    else()
      message(STATUS "libpng not found")
    endif()

    add_executable(see ${CMAKE_SOURCE_DIR}/Apps/see.cpp)

    if(PNG_FOUND)
      target_compile_options(see PRIVATE -DWITH_LIBPNG)
    endif()

    target_include_directories(
      see
      PRIVATE ${OPENGL_INCLUDE_DIR} ${GLU_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS}
              ${nlohmann-json_SOURCE_DIR}/include ${PNG_INCLUDE_DIRS})

    target_link_directories(see PRIVATE ${GLUT_LIBRARY_DIRS} ${PNG_LIBRARY_DIRS})

    target_link_libraries(see PRIVATE OpenGL::GL OpenGL::GLU ${GLUT_LIBRARIES}
                                     Rockable_core PNG::PNG)

  else() # Darwin =========

    if(NOT PkgConfig_FOUND)
      find_package(PkgConfig REQUIRED)
    endif()

    # Check for SDL2 instead of GLFW
    pkg_check_modules(PNG REQUIRED libpng16)

    if(PNG_FOUND)
      message(STATUS "PNG has been found with PkgConfig")
      message(STATUS "PNG_LIBRARY_DIRS = " ${PNG_LIBRARY_DIRS})
      message(STATUS "PNG_INCLUDE_DIRS = " ${PNG_INCLUDE_DIRS})
      message(STATUS "PNG_LIBRARIES = " ${PNG_LIBRARIES})
    endif()

    add_executable(see ${CMAKE_SOURCE_DIR}/Apps/see.cpp)

    if(PNG_FOUND)
      target_compile_options(see PRIVATE -DWITH_LIBPNG)
    endif()

    target_include_directories(
      see PRIVATE ${GL_INCLUDE_DIRS} ${GLU_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS}
                  ${nlohmann-json_SOURCE_DIR}/include ${PNG_INCLUDE_DIRS})

    target_link_directories(see PRIVATE ${GL_LIBRARY_DIRS} ${GLU_LIBRARY_DIRS}
                            ${GLUT_LIBRARY_DIRS} ${PNG_LIBRARY_DIRS})

    target_link_libraries(
      see PUBLIC ${GL_LIBRARIES} ${GLU_LIBRARIES} ${GLUT_LIBRARIES}
                 ${PNG_LIBRARIES} Rockable_core)

  endif() #=========

  install(TARGETS see)

endif()

# ##############################################################################
# interactive preprocessing (TODO see3 -> will be renamed something like 'craft' or 'seedit' or 'reveal'
# or 'look') Before, the use of glfw will be replaced by glut (to avoid issues
# with macosx)
if(ROCKABLE_COMPILE_SEE3)

  set(SEE3_FILES
      ${CMAKE_SOURCE_DIR}/Apps/see3.cpp
      ${imgui_SOURCE_DIR}/imgui.cpp
      ${imgui_SOURCE_DIR}/imgui_draw.cpp
      ${imgui_SOURCE_DIR}/imgui_tables.cpp
      ${imgui_SOURCE_DIR}/imgui_widgets.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)

  add_executable(see3 ${SEE3_FILES})

  if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux") # ===========

    target_include_directories(
      see3
      PRIVATE ${OPENGL_INCLUDE_DIR} ${imgui_SOURCE_DIR}
             ${imgui_SOURCE_DIR}/backends ${nlohmann-json_SOURCE_DIR}/include)

    target_link_libraries(see3 PRIVATE glfw OpenGL::GL OpenGL::GLU Rockable_core)

  else() # For APPLE  ===========

    target_include_directories(
      see3
      PRIVATE ${GL_INCLUDE_DIRS} ${GLU_INCLUDE_DIRS} ${GLFW_INCLUDE_DIRS}
             ${GLFW_INCLUDE_DIRS} ${imgui_SOURCE_DIR}
             ${imgui_SOURCE_DIR}/backends)

    target_link_directories(see3 PRIVATE ${GL_LIBRARY_DIRS} ${GLU_LIBRARY_DIRS}
                            ${GLFW_LIBRARY_DIRS})

    target_link_libraries(see3 PRIVATE ${GL_LIBRARIES} ${GLU_LIBRARIES}
                                      ${GLFW_LIBRARIES} Rockable_core)

  endif() # ===========

  install(TARGETS see3)

endif() # end of ROCKABLE_COMPILE_SEE3

# ##############################################################################
# interactive preprocessing (TODO see4 -> will be renamed something like 'craft'
# or 'look') Before, the use of glfw will be replaced by glut (to avoid issues
# with macosx)
if(ROCKABLE_COMPILE_SEE4)

  set(SEE4_FILES
      ${CMAKE_SOURCE_DIR}/Apps/see4.cpp
      ${imgui_SOURCE_DIR}/imgui.cpp
      ${imgui_SOURCE_DIR}/imgui_demo.cpp
      ${imgui_SOURCE_DIR}/imgui_draw.cpp
      ${imgui_SOURCE_DIR}/imgui_tables.cpp
      ${imgui_SOURCE_DIR}/imgui_widgets.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_sdl2.cpp
      ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)

  add_executable(see4 ${SEE4_FILES})

  if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")

    target_include_directories(
      see4
      PUBLIC ${OPENGL_INCLUDE_DIR} ${SDL2_INCLUDE_DIRS} ${imgui_SOURCE_DIR}
             ${imgui_SOURCE_DIR}/backends ${nlohmann-json_SOURCE_DIR}/include)

    target_link_libraries(see4 PUBLIC SDL2 OpenGL::GL OpenGL::GLU Rockable_core)

  else()

    target_include_directories(
      see4
      PUBLIC ${GL_INCLUDE_DIRS} ${GLU_INCLUDE_DIRS} ${SDL2_INCLUDE_DIRS}
             ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends
             ${nlohmann-json_SOURCE_DIR}/include)

    target_link_directories(see4 PUBLIC ${GL_LIBRARY_DIRS} ${GLU_LIBRARY_DIRS}
                            ${SDL2_LIBRARY_DIRS})

    target_link_libraries(see4 PUBLIC ${GL_LIBRARIES} ${GLU_LIBRARIES}
                                      ${SDL2_LIBRARIES} Rockable_core)
  endif()

  install(TARGETS see4)

endif() # end of ROCKABLE_COMPILE_SEE4

# ##############################################################################
#
if(ROCKABLE_COMPILE_CONF2VTK)
  add_executable(conftovtk ${CMAKE_SOURCE_DIR}/Apps/conftovtk.cpp)
  target_link_libraries(conftovtk Rockable_core)
  # install(TARGETS conftovtk DESTINATION ./)
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
