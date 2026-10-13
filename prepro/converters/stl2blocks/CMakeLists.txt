add_executable(stl2blocks stl2blocks.cpp)
target_link_libraries(stl2blocks PUBLIC Rockable_core)

install(TARGETS stl2blocks)
