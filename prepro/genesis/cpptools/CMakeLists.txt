file(GLOB cpptools_SRC "*.cpp")
message(STATUS "cpptools_SRC: ${cpptools_SRC}")
#file(GLOB cpptools_HDR "*.hpp")

add_executable(generator ${cpptools_SRC})
install(TARGETS generator DESTINATION ./)
