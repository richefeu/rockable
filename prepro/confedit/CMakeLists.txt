cmake_minimum_required(VERSION 3.26)
project(confedit)

# Find FLTK package
find_package(FLTK REQUIRED)

if(FLTK_FOUND)
	message(STATUS "FLTK_INCLUDE_DIR: ${FLTK_INCLUDE_DIR}")
	message(STATUS "FLTK_LIBRARIES: ${FLTK_LIBRARIES}")
endif()

# prepro/CMakeLists.txt exports this; fall back to the sibling directory so that
# this CMakeLists can still be configured on its own.
if(NOT rockable_lang_DIR)
  set(rockable_lang_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../common")
endif()

# Executable
add_executable(confedit confedit.cpp)

# Include directories
target_include_directories( confedit PRIVATE
    ${FLTK_INCLUDE_DIR}
    ${rockable_lang_DIR}
    ~/toofus
    ../genesis/cpptools
)

# Link FLTK libraries
target_link_libraries(confedit  ${FLTK_LIBRARIES})

install(TARGETS confedit)



