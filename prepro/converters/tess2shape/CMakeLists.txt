add_executable(tess2shape main.cpp tess2shape.cpp)
target_link_libraries(tess2shape PUBLIC Rockable_core)

install(TARGETS tess2shape)
