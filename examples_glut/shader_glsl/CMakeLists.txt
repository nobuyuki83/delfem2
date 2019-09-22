project(shader_glsl)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g -std=c++11" )
endif()
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")
add_definitions(-DUSE_GLEW=1)

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

set(GLAD_INCLUDE_DIR "../glad2/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${GLEW_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_SRC}/gl_tex.cpp      ${DELFEM2_INC}/gl_tex.h
  ${DELFEM2_SRC}/gl2_funcs.cpp   ${DELFEM2_INC}/gl2_funcs.h
  ${DELFEM2_SRC}/gl24_camera.cpp ${DELFEM2_INC}/gl24_camera.h
  ${DELFEM2_SRC}/gl24_funcs.cpp  ${DELFEM2_INC}/gl24_funcs.h
  main.cpp  
)

target_link_libraries(${PROJECT_NAME}
  ${GLEW_LIBRARY}
  ${GLUT_LIBRARY}
  ${OPENGL_LIBRARY})

