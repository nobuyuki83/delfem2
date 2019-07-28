project(viewer3d_glfw)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
endif()

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

set(GLFW_DIR "../../src_cpp/external/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

find_package(OpenGL REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h         ${DELFEM2_SRC}/vec2.cpp             
  ${DELFEM2_INC}/vec3.h         ${DELFEM2_SRC}/vec3.cpp             
  ${DELFEM2_INC}/mat3.h         ${DELFEM2_SRC}/mat3.cpp             
  ${DELFEM2_INC}/quat.h         ${DELFEM2_SRC}/quat.cpp             
  
  ${DELFEM2_INC}/gl_funcs.h     ${DELFEM2_SRC}/gl_funcs.cpp             
  ${DELFEM2_INC}/gl_v23q.h      ${DELFEM2_SRC}/gl_v23q.cpp              

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

