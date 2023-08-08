
file(GLOB sphere_packer_SRCS "*.cpp")
file(GLOB sphere_packer_HDRS "*.hpp")

add_executable(seePack seePack.cpp)
target_link_libraries(seePack ${GLUT_LIBRARY} glut GLU GL)
install(TARGETS seePack DESTINATION ./)

include_directories(${cpptools_HDRS})
add_executable(pack SpherePacker.cpp)
install(TARGETS pack DESTINATION ./)