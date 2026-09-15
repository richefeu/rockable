add_executable(shape2mesh shape2mesh.cpp)
target_link_libraries(shape2mesh PUBLIC Rockable_core)

install(TARGETS shape2mesh)
