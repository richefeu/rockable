FIND_PACKAGE(TIFF REQUIRED)
FIND_PACKAGE(X11 REQUIRED)

add_executable(tif2rockable tif2rockable.cpp)
target_link_libraries(tif2rockable ${TIFF_LIBRARIES} ${X11_LIBRARIES})
install(TARGETS tif2rockable DESTINATION bin)      
