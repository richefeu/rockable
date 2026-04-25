file(GLOB shapeSurvey_SRCS "*.cpp")
add_executable(shapeSurvey shapeSurvey.cpp)

if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")

  target_include_directories(shapeSurvey PRIVATE ${OPENGL_INCLUDE_DIR}
                                                ${GLUT_INCLUDE_DIRS})

  target_link_libraries(shapeSurvey PRIVATE OpenGL::GL OpenGL::GLU
                                           ${GLUT_LIBRARIES} Rockable_core)

else()

  target_include_directories(
    shapeSurvey PRIVATE ${OPENGL_INCLUDE_DIR} ${GLU_INCLUDE_DIRS}
                        ${GLUT_INCLUDE_DIRS} /usr/local/include)

  target_link_libraries(shapeSurvey PRIVATE OpenGL::GL OpenGL::GLU ${GLUT_LIBRARIES} Rockable_core)

endif()

install(TARGETS shapeSurvey)
