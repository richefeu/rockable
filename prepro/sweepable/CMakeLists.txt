add_executable(sweepable sweepable.cpp)
install(TARGETS sweepable)

# Install the script file (copy during configuration)
#file(COPY batchable.sh DESTINATION ${CMAKE_BINARY_DIR})
