project(paramgeo2d)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )

set(DELFEM2_INCLUDE_DIR "../../../include")
set(DELFEM2_SRC         "../../../src")
set(DELFEM2_INC         "../../../include/delfem2")


find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_SRC}/mat3.cpp              ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/vec2.cpp              ${DELFEM2_INC}/vec2.h
  ${DELFEM2_SRC}/vec3.cpp              ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/quat.cpp              ${DELFEM2_INC}/quat.h  

  ${DELFEM2_SRC}/paramgeo_v23.cpp      ${DELFEM2_INC}/paramgeo_v23.h

  ${DELFEM2_SRC}/funcs_gl.cpp
  ${DELFEM2_SRC}/v23_gl.cpp
  ${DELFEM2_SRC}/camera_gl.cpp         ${DELFEM2_INC}/camera_gl.h

  main.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${GLUT_LIBRARY}
  ${OPENGL_LIBRARY}
)

add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/../../bin/${PROJECT_NAME} 
)

