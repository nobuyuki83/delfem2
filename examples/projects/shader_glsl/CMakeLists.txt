project(test_glsl)

cmake_minimum_required(VERSION 2.8)
set( CMAKE_C_FLAGS "-Wall -g -std=c++11" )

set(DELFEM2_INCLUDE "../../../include")
set(DELFEM2_SRC     "../../../src")
set(DELFEM2_INC     "../../../include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${GLEW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE}
)

add_executable(${PROJECT_NAME}
  main.cpp
  ${DELFEM2_SRC}/funcs_gl.cpp
  ${DELFEM2_SRC}/camera_gl.cpp    ${DELFEM2_INC}/camera_gl.h
)

target_link_libraries(${PROJECT_NAME}  
  ${GLUT_LIBRARY} 
  ${GLEW_LIBRARY}  
  ${OPENGL_LIBRARY})

add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/../../bin/${PROJECT_NAME} 
)
