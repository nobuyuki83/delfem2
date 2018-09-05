project(fbx_read)

cmake_minimum_required(VERSION 2.8)
set( CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )

set(DELFEM2_INCLUDE_DIR "../../../src_cpp/include")
set(DELFEM2_SRC         "../../../src_cpp/cpp")
set(DELFEM2_INC         "../../../src_cpp/include/delfem2")

IF(APPLE)
  MESSAGE("platform apple os x")
  SET(FBX_ROOT "/Applications/Autodesk/FBX SDK/2019.0")
  SET(FBX_INCLUDE_DIR ${INCLUDE_DIR} ${FBX_ROOT}/include)
  SET(FBX_LINK_LIBS ${FBX_ROOT}/lib/clang/release/libfbxsdk.a)
ENDIF(APPLE)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${FBX_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)


add_executable(${PROJECT_NAME}
  ${DELFEM2_SRC}/funcs.cpp          ${DELFEM2_INC}/funcs.h
  ${DELFEM2_SRC}/vec2.cpp           ${DELFEM2_INC}/vec2.h
  ${DELFEM2_SRC}/vec3.cpp           ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/mat3.cpp           ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/quat.cpp           ${DELFEM2_INC}/quat.h
  ${DELFEM2_SRC}/mshio.cpp          ${DELFEM2_INC}/mshio.h
  ${DELFEM2_SRC}/msh.cpp          ${DELFEM2_INC}/msh.h

  ${DELFEM2_SRC}/v23m3q.cpp         ${DELFEM2_INC}/v23m3q.h
  ${DELFEM2_SRC}/rigmesh.cpp        ${DELFEM2_INC}/rigmesh.h

  ${DELFEM2_SRC}/camera_gl.cpp      ${DELFEM2_INC}/camera_gl.h
  ${DELFEM2_SRC}/funcs_gl.cpp       ${DELFEM2_INC}/funcs_gl.h
  ${DELFEM2_SRC}/color_gl.cpp       ${DELFEM2_INC}/color_gl.h
  ${DELFEM2_SRC}/v23_gl.cpp       ${DELFEM2_INC}/v23_gl.h

   ${DELFEM2_SRC}/../external/io_fbx.cpp  ${DELFEM2_INC}/../../external/io_fbx.h

  main.cpp
)


target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY}
  ${FBX_LINK_LIBS}
)

add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/../../bin/${PROJECT_NAME} 
)

