project(cloth_selfcollision)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall  -Wno-deprecated-declarations -std=c++11" )

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
  ${DELFEM2_SRC}/vec3.cpp                  ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/mshtopo.cpp               ${DELFEM2_INC}/mshtopo.h
  ${DELFEM2_SRC}/bvh.cpp                   ${DELFEM2_INC}/bvh.h
  ${DELFEM2_SRC}/collisiontri_v3bvh.cpp    ${DELFEM2_INC}/collisiontri_v3bvh.h

  ${DELFEM2_SRC}/fem_ematrix.cpp
  ${DELFEM2_SRC}/fem.cpp
  ${DELFEM2_SRC}/matrix_sparse.cpp
  ${DELFEM2_SRC}/ilu_sparse.cpp

  ${DELFEM2_SRC}/gl_funcs.cpp
  ${DELFEM2_SRC}/gl_camera.cpp             ${DELFEM2_INC}/gl_camera.h

  ${DELFEM2_SRC}/cloth_selfcollision.cpp   ${DELFEM2_INC}/cloth_selfcollision.h
  ${DELFEM2_INC}/cloth_internal.h  
  main.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY}
)
