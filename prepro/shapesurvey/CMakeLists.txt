file(GLOB shapeSurvey_SRCS "*.cpp")
add_executable(shapeSurvey shapeSurvey.cpp)

target_include_directories(shapeSurvey PUBLIC 
	${OPENGL_INCLUDE_DIR}
	${GLUT_INCLUDE_DIRS}
)

target_link_directories(shapeSurvey PUBLIC ${GLUT_LIBRARY_DIRS})

if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
  target_link_libraries (shapeSurvey PUBLIC OpenGL::GL OpenGL::GLU ${GLUT_LIBRARIES} Rockable_core)
else()
  target_link_libraries (shapeSurvey PUBLIC OpenGL::GL OpenGL::GLU ${GLUT_LIBRARIES} Rockable_core)
endif()

install(TARGETS shapeSurvey DESTINATION ./)    