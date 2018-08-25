cmake_minimum_required(VERSION 3.0)

project(pydelfem2)
set(PYBIND11_CPP_STANDARD -std=c++14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O2")

set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

set(DELFEM2_INCLUDE_DIR "../include")
set(DELFEM2_SRC         "../src")
set(DELFEM2_INC         "../include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

pybind11_add_module(dfm2 
  ${DELFEM2_SRC}/vec2.cpp        ${DELFEM2_INC}/vec2.h
  ${DELFEM2_SRC}/vec3.cpp        ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/mat3.cpp        ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/quat.cpp        ${DELFEM2_INC}/quat.h
  ${DELFEM2_SRC}/msh.cpp         ${DELFEM2_INC}/msh.h
  ${DELFEM2_SRC}/mshtopo.cpp     ${DELFEM2_INC}/mshtopo.h  
  ${DELFEM2_SRC}/mshio.cpp       ${DELFEM2_INC}/mshio.h
  ${DELFEM2_SRC}/funcs.cpp       ${DELFEM2_INC}/funcs.h
  ${DELFEM2_SRC}/voxel.cpp       ${DELFEM2_INC}/voxel.h

  ${DELFEM2_SRC}/v23m3q.cpp      ${DELFEM2_INC}/v23m3q.h

  ${DELFEM2_SRC}/v23_gl.cpp      ${DELFEM2_INC}/v23_gl.h
  ${DELFEM2_SRC}/funcs_gl.cpp    ${DELFEM2_INC}/funcs_gl.h
  ${DELFEM2_SRC}/color_gl.cpp    ${DELFEM2_INC}/color_gl.h

  ${DELFEM2_INC}/mshtopoio_gl.h
  py_mshtopoio_gl.cpp
  py_delfem2.cpp
)

target_link_libraries(dfm2
PRIVATE
  ${GLUT_LIBRARY}
  ${OPENGL_LIBRARY}
)

add_custom_command(TARGET dfm2
POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/dfm2.*
  ${PROJECT_SOURCE_DIR}/dfm2/dfm2.so
  COMMENT "Post Event"
)
