project(shader_glsl)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g -std=c++11" )
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${GLEW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_SRC}/gl_tex.cpp       ${DELFEM2_INC}/gl_tex.h
  ${DELFEM2_SRC}/gl_funcs.cpp     ${DELFEM2_INC}/gl_funcs.h
  ${DELFEM2_SRC}/gl_camera.cpp    ${DELFEM2_INC}/gl_camera.h
  main.cpp  
)

target_link_libraries(${PROJECT_NAME}  
  ${GLUT_LIBRARY} 
  ${GLEW_LIBRARY}  
  ${OPENGL_LIBRARY})
