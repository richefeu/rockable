
file(GLOB confedit_SRCS "*.cpp")
file(GLOB confedit_HDRS "*.hpp")

FIND_PACKAGE(FLTK REQUIRED)

include_directories(${cpptools_HDRS})
add_executable(confedit confedit.cpp)
target_link_libraries(confedit Rockable_core ${FLTK_LIBRARIES})
install(TARGETS confedit DESTINATION bin)
