project(401_read_gltf)
cmake_minimum_required(VERSION 2.8)

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")
set(DELFEM2_EXT         "../../src_cpp/external")
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(GLFW_DIR "../../src_cpp/external/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

find_package(OpenGL REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
 ${DELFEM2_INC}/vec2.h          ${DELFEM2_SRC}/vec2.cpp
  ${DELFEM2_INC}/vec3.h         ${DELFEM2_SRC}/vec3.cpp
  ${DELFEM2_INC}/mat3.h         ${DELFEM2_SRC}/mat3.cpp
  ${DELFEM2_INC}/quat.h         ${DELFEM2_SRC}/quat.cpp
  ${DELFEM2_INC}/mshmisc.h      ${DELFEM2_SRC}/mshmisc.cpp
  ${DELFEM2_INC}/mshio.h        ${DELFEM2_SRC}/mshio.cpp
  ${DELFEM2_INC}/funcs.h        ${DELFEM2_SRC}/funcs.cpp
  ${DELFEM2_INC}/camera.h       ${DELFEM2_SRC}/camera.cpp
  ${DELFEM2_INC}/color.h        ${DELFEM2_SRC}/color.cpp
  
  ${DELFEM2_INC}/v23m3q.h       ${DELFEM2_SRC}/v23m3q.cpp
  ${DELFEM2_INC}/rig_v3q.h      ${DELFEM2_SRC}/rig_v3q.cpp

  ${DELFEM2_INC}/opengl/gl2_funcs.h          ${DELFEM2_INC}/opengl/gl2_funcs.cpp
  ${DELFEM2_INC}/opengl/gl2_color.h          ${DELFEM2_INC}/opengl/gl2_color.cpp
  ${DELFEM2_INC}/opengl/gl2_v23.h            ${DELFEM2_INC}/opengl/gl2_v23.cpp
  ${DELFEM2_INC}/opengl/gl_rig_v23q.h        ${DELFEM2_INC}/opengl/gl_rig_v23q.cpp
  ${DELFEM2_INC}/opengl/glfw_viewer.hpp      ${DELFEM2_INC}/opengl/glfw_viewer.cpp

  ${DELFEM2_EXT}/io_gltf.h      ${DELFEM2_EXT}/io_gltf.cpp
  main.cpp
)

IF(APPLE)
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g -std=c++11")
  find_library( COCOA_LIBRARY Cocoa )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    ${COCOA_LIBRARY}
    ${IOKIT_LIBRARY}
    ${COREVIDEO_LIBRARY}
  )
ELSEIF(UNIX)
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g -std=c++11")
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
  )
ELSEIF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2")  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "64bit")
    target_link_libraries(${PROJECT_NAME} 
      ${GLFW_DIR}/buildVS64/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  else()
    message(STATUS "32bit")
    target_link_libraries(${PROJECT_NAME} 
      ${GLFW_DIR}/buildVS32/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  endif()
ENDIF()




