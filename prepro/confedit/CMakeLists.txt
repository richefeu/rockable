cmake_minimum_required(VERSION 3.26)
project(confedit)

# Find FLTK package
find_package(FLTK REQUIRED)

if(FLTK_FOUND)
	message(STATUS "FLTK_INCLUDE_DIR: ${FLTK_INCLUDE_DIR}")
	message(STATUS "FLTK_LIBRARIES: ${FLTK_LIBRARIES}")
endif()

# Executable
add_executable(confedit confedit.cpp)

# Include directories
target_include_directories( confedit PRIVATE
    ${FLTK_INCLUDE_DIR}
    ~/toofus
    ../genesis/cpptools
)

# Link FLTK libraries
target_link_libraries(confedit  ${FLTK_LIBRARIES})

install(TARGETS confedit)



