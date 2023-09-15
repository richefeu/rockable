#file(GLOB stl2shape_SRCS "*.cpp")
add_executable(stl2shape stl2shape.cpp)
#target_link_libraries(stl2shape PUBLIC Rockable_core)
install(TARGETS stl2shape DESTINATION ./)
