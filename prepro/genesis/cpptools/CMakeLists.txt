file(GLOB cpptools_SRC "*.cpp")
file(GLOB cpptools_HDR "*.hpp")

add_executable(generator ${cpptools_SRC} ${cpptools_HDRS})
install(TARGETS generator DESTINATION ./)
