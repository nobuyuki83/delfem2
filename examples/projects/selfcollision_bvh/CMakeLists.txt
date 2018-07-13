project(selfcollision_bvh)

cmake_minimum_required(VERSION 2.8)
set( CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -std=c++11" )

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
  main.cpp  
  ${DELFEM2_SRC}/vec3.cpp
  ${DELFEM2_SRC}/msh.cpp
  ${DELFEM2_SRC}/mshtopo.cpp
                                            ${DELFEM2_INC}/bv.h
  ${DELFEM2_SRC}/bvh.cpp                    ${DELFEM2_INC}/bvh.h

  ${DELFEM2_SRC}/collisiontri_v3bvh.cpp     ${DELFEM2_INC}/collisiontri_v3bvh.h

  ${DELFEM2_SRC}/funcs_gl.cpp
  ${DELFEM2_SRC}/camera_gl.cpp              ${DELFEM2_INC}/camera_gl.h
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