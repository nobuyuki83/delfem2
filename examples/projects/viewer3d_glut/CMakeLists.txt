project(viewer3d_glut)

cmake_minimum_required(VERSION 2.8)
set( CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )

set(DELFEM2_INCLUDE "../../../include")
set(DELFEM2_SRC     "../../../src")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE}
)

add_executable(${PROJECT_NAME}
  main.cpp
  ${DELFEM2_SRC}/quat.cpp
  ${DELFEM2_SRC}/funcs_gl.cpp
  ${DELFEM2_SRC}/color_gl.cpp
  ${DELFEM2_SRC}/camera_gl.cpp           
  ${DELFEM2_INCLUDE}/delfem2/camera_gl.h
)

target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY})

add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD  
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${PROJECT_NAME}  
  ${PROJECT_SOURCE_DIR}/../../bin/${PROJECT_NAME}
)
