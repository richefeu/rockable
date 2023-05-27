cmake_minimum_required(VERSION 3.26)

message(STATUS "CMAKE_HOST_SYSTEM_NAME = " ${CMAKE_HOST_SYSTEM_NAME})
 
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_CXX_COMPILER "g++-13") # do not use AppleClang
  message(STATUS "CMAKE_CXX_COMPILER = " ${CMAKE_CXX_COMPILER})
endif()

project(Rockable LANGUAGES CXX VERSION 0.9)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ROCKABLE_USE_FT_CORR "add objectivity correction to tangent forces" OFF)
option(ROCKABLE_ENABLE_PROFILING "enable the time profiling" OFF)

set(CMAKE_CXX_FLAGS "-O3 -Wno-deprecated-declarations -Wall -Wextra -lpthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

if(ROCKABLE_USE_FT_CORR)
	add_compile_options(-DFT_CORR)
endif()

if(ROCKABLE_ENABLE_PROFILING)
	add_compile_options(-DENABLE_PROFILING)
endif()

######################################################################
##### Fetching TOOFUS
include(FetchContent)
FetchContent_Declare (
	toofus
  GIT_REPOSITORY https://github.com/richefeu/toofus.git
  GIT_TAG        main
)

FetchContent_GetProperties(toofus)
if(NOT toofus_POPULATED)
  message(STATUS "Fetching toofus")
  FetchContent_Populate(toofus)
	#add_subdirectory(${toofus_SOURCE_DIR})
endif()

######################################################################
##### Fetching spdlog
FetchContent_Declare (
	spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.11.0
)

FetchContent_GetProperties(spdlog)
if(NOT spdlog_POPULATED)
  message(STATUS "Fetching spdlog")
  FetchContent_Populate(spdlog)
	#add_subdirectory(${spdlog_SOURCE_DIR})
endif()

######################################################################
##### tclap
FetchContent_Declare (
	tclap
  GIT_REPOSITORY https://github.com/xguerin/tclap.git
  GIT_TAG        v1.2.5
)

FetchContent_GetProperties(tclap)
if(NOT tclap_POPULATED)
  message(STATUS "Fetching tclap")
  FetchContent_Populate(tclap)
	#add_subdirectory(${tclap_SOURCE_DIR})
endif()

######################################################################
##### nlohmann-json
FetchContent_Declare (
	nlohmann-json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.2
)

FetchContent_GetProperties(nlohmann-json)
if(NOT nlohmann-json_POPULATED)
  message(STATUS "Fetching nlohmann-json")
  FetchContent_Populate(nlohmann-json)
	#add_subdirectory(${nlohmann-json_SOURCE_DIR})
endif()

######################################################################
##### imgui
FetchContent_Declare (
	imgui
  GIT_REPOSITORY https://github.com/ocornut/imgui.git
  GIT_TAG        v1.89.5
)

FetchContent_GetProperties(imgui)
if(NOT imgui_POPULATED)
  message(STATUS "Fetching imgui")
  FetchContent_Populate(imgui)
	#add_subdirectory(${imgui_SOURCE_DIR})
endif()

######################################################################
##### OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

######################################################################
##### OpenGL
find_package(OpenGL REQUIRED)
if(OPENGL_FOUND)
  message(STATUS "OPENGL_INCLUDE_DIR = " ${OPENGL_INCLUDE_DIR})
  message(STATUS "OPENGL_LIBRARIES = " ${OPENGL_LIBRARIES})
endif()

######################################################################
#### glut
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")
  # here, we do not uses find_package(glut) so that the version installed by apple is NOT used
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(GLUT REQUIRED glut>=3.4.0)
  if(GLUT_FOUND)
    message(STATUS "GLUT as been found with PkgConfig")
		
		message(STATUS "GLUT_LIBRARY_DIRS = " ${GLUT_LIBRARY_DIRS})
	  message(STATUS "GLUT_INCLUDE_DIRS = " ${GLUT_INCLUDE_DIRS})
		message(STATUS "GLUT_LIBRARIES = " ${GLUT_LIBRARIES})
  endif() 

else()

  find_package(GLUT REQUIRED)
  if(GLUT_FOUND)
    message(STATUS "GLUT_INCLUDE_DIRS = " ${GLUT_INCLUDE_DIRS})
	  message(STATUS "GLUT_LIBRARIES = " ${GLUT_LIBRARIES})
  endif()
	
endif()

######################################################################
#### glfw
#find_package(GLFW3 3.3 REQUIRED)
#if(GLFW3_FOUND)
#  message(STATUS "GLFW3_INCLUDE_DIRS = " ${GLFW_INCLUDE_DIR})
#  message(STATUS "GLFW3_LIBRARIES = " ${GLFW3_LIBRARIES})
#endif()

######################################################################
##### libRockable_core
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR} 
	${toofus_SOURCE_DIR}
	${spdlog_SOURCE_DIR}/include
	${tclap_SOURCE_DIR}/include
)

file(GLOB core_object_files ${CMAKE_CURRENT_SOURCE_DIR}
  BodyForces/*.cpp
  Core/*.cpp
  DataExtractors/*.cpp
  ForceLaws/*.cpp
  PostProcessors/*.cpp
  PreproCommands/*.cpp
  ProcessingTools/*.cpp
)

######################################################################
##### rockable executable (to run simulations from the console)
add_library(Rockable_core STATIC ${core_object_files})
add_executable(rockable Apps/run.cpp)
target_link_libraries(rockable PUBLIC Rockable_core)

######################################################################
##### visualiser see
add_executable(see Apps/see.cpp)
target_include_directories(see PUBLIC 
	${OPENGL_INCLUDE_DIR} 
	${GLUT_INCLUDE_DIRS}
	${nlohmann-json_SOURCE_DIR}/include
)
target_link_directories(see PUBLIC ${GLUT_LIBRARY_DIRS})
target_link_libraries (see PUBLIC OpenGL::GL ${GLUT_LIBRARIES} Rockable_core)

######################################################################
##### interactive preprocessing (TODO see3 -> will be renamed)
#set(imgui_object_files
#  ${imgui_SOURCE_DIR}/imgui.cpp
#  ${imgui_SOURCE_DIR}/imgui_demo.cpp
#  ${imgui_SOURCE_DIR}/imgui_draw.cpp
#  ${imgui_SOURCE_DIR}/imgui_tables.cpp
#  ${imgui_SOURCE_DIR}/imgui_widgets.cpp
#  ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
#  ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp
#)
#add_library(imgui STATIC ${imgui_object_files})
#target_include_directories(imgui PUBLIC 
#	${OPENGL_INCLUDE_DIR} 
#	${imgui_SOURCE_DIR}
#)
#add_executable(see3 Apps/see3.cpp)
#target_include_directories(see3 PUBLIC 
#	${OPENGL_INCLUDE_DIR} 
#	${nlohmann-json_SOURCE_DIR}/include
#	${imgui_SOURCE_DIR}
#)
#target_link_libraries (see3 PUBLIC OpenGL::GL glfw imgui Rockable_core)
