project(adf3d)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -std=c++11" )
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
  ${DELFEM2_SRC}/vec3.cpp                ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/msh.cpp
  ${DELFEM2_SRC}/mshio.cpp
  
  ${DELFEM2_SRC}/adf.cpp
  ${DELFEM2_SRC}/sdf.cpp
  ${DELFEM2_SRC}/hash3d.cpp

  ${DELFEM2_SRC}/gl_funcs.cpp
  ${DELFEM2_SRC}/gl_camera.cpp           ${DELFEM2_INC}/gl_camera.h
)

target_link_libraries(${PROJECT_NAME}
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY}
)
