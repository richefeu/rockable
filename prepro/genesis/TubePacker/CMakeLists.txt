
file(GLOB sphere_packer_SRCS "*.cpp")
file(GLOB sphere_packer_HDRS "*.hpp")

add_executable(seePack-tube seePack.cpp)
target_link_libraries(seePack-tube ${GLUT_LIBRARY} glut GLU GL)
install(TARGETS seePack-tube DESTINATION ./)

include_directories(${cpptools_HDRS})
add_executable(pack-tube TubePacker.cpp)
install(TARGETS pack-tube DESTINATION ./)