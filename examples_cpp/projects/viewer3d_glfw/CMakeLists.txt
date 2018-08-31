project(viewer3d_glfw)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")

set(DELFEM2_INCLUDE_DIR "../../../include")
set(DELFEM2_SRC         "../../../src")
set(DELFEM2_INC         "../../../include/delfem2")

set(GLFW_DIR "../../../src_cpp/external/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

find_package(OpenGL REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_SRC}/vec2.cpp            ${DELFEM2_INC}/vec2.h
  ${DELFEM2_SRC}/vec3.cpp            ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/mat3.cpp            ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/quat.cpp            ${DELFEM2_INC}/quat.h

  ${DELFEM2_SRC}/funcs_gl.cpp        ${DELFEM2_INC}/funcs_gl.h
  ${DELFEM2_SRC}/v23_gl.cpp          ${DELFEM2_INC}/v23_gl.h

  main.cpp
)

IF(APPLE)
  find_library( COCOA_LIBRARY Cocoa )
  find_library( OpenGL_LIBRARY OpenGL )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  SET( EXTRA_LIBS ${COCOA_LIBRARY} ${OpenGL_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} )
ELSEIF(UNIX)
  SET( EXTRA_LIBS   -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp )
ENDIF()

target_link_libraries(${PROJECT_NAME} 
  ${OPENGL_LIBRARY}
  ${GLFW_DIR}/src/libglfw3.a
  ${EXTRA_LIBS}
 
)

add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/../../bin/${PROJECT_NAME} 
)
