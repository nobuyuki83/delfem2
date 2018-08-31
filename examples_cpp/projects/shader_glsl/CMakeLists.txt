project(test_glsl)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g -std=c++11" )

set(DELFEM2_INCLUDE "../../../src_cpp/include")
set(DELFEM2_SRC     "../../../src_cpp/cpp")
set(DELFEM2_INC     "../../../src_cpp/include/delfem2")

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
  ${DELFEM2_SRC}/funcs_gl.cpp     ${DELFEM2_INC}/funcs_gl.h
  ${DELFEM2_SRC}/camera_gl.cpp    ${DELFEM2_INC}/camera_gl.h
  main.cpp  
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
