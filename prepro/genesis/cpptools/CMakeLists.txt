
file(GLOB cpptools_SRCS "*.cpp")
file(GLOB cpptools_HDRS "*.hpp")

add_executable(generator ${cpptools_SRCS} ${cpptools_HDRS})
install(TARGETS generator DESTINATION bin)
