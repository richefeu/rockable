
file(GLOB sphere_packer_SRCS "*.cpp")
file(GLOB sphere_packer_HDRS "*.hpp")

add_executable(seePack seePack.cpp)
target_include_directories(seePack PUBLIC ${COMMON_HEADERS})
target_link_libraries(seePack ${GLUT_LIBRARY} glut GLU GL)
install(TARGETS seePack DESTINATION bin)


add_executable(pack SpherePacker.cpp)
target_include_directories(pack PUBLIC ${cpptools_SOURCE})
install(TARGETS pack DESTINATION bin)