
file(GLOB sphere_packer_SRCS "*.cpp")
file(GLOB sphere_packer_HDRS "*.hpp")

add_executable(seePack-tube seePack.cpp)
target_include_directories(seePack-tube PUBLIC ${COMMON_HEADERS})
target_link_libraries(seePack-tube ${GLUT_LIBRARY} glut GLU GL)
install(TARGETS seePack-tube DESTINATION bin)

add_executable(pack-tube TubePacker.cpp)
target_include_directories(pack-tube PUBLIC ${cpptools_SOURCE})
install(TARGETS pack-tube DESTINATION bin)