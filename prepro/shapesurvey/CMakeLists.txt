file(GLOB shapeSurvey_SRCS "*.cpp")
add_executable(shapeSurvey shapeSurvey.cpp)

if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")

  target_include_directories(shapeSurvey PUBLIC 
  	${OPENGL_INCLUDE_DIR}
  	${GLUT_INCLUDE_DIRS}
  )

  target_link_directories(shapeSurvey PUBLIC ${GLUT_LIBRARY_DIRS})

  target_link_libraries (shapeSurvey PUBLIC OpenGL::GL OpenGL::GLU ${GLUT_LIBRARIES} Rockable_core)
  
else()

  target_include_directories(shapeSurvey PUBLIC 
  	${GL_INCLUDE_DIRS}
    ${GLU_INCLUDE_DIRS}
  	${GLUT_INCLUDE_DIRS}
  )

  target_link_directories(shapeSurvey PUBLIC     
    ${GL_LIBRARY_DIRS}
    ${GLU_LIBRARY_DIRS}
    ${GLUT_LIBRARY_DIRS}
  )

  target_link_libraries (shapeSurvey PUBLIC ${GL_LIBRARIES} ${GLU_LIBRARIES} ${GLUT_LIBRARIES} Rockable_core)
  
endif()

install(TARGETS shapeSurvey)
