project(dyntri2d_triangulation)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")

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
  ${DELFEM2_SRC}/vec2.cpp           ${DELFEM2_INC}/vec2.h
  ${DELFEM2_SRC}/vec3.cpp           ${DELFEM2_INC}/vec3.h  
  ${DELFEM2_SRC}/dyntri.cpp         ${DELFEM2_INC}/dyntri.h
  ${DELFEM2_SRC}/dyntri_v2.cpp      ${DELFEM2_INC}/dyntri_v2.h

  ${DELFEM2_SRC}/gl_funcs.cpp      ${DELFEM2_INC}/gl_funcs.h  
  ${DELFEM2_SRC}/gl_cad_dyntri_v23.cpp      
  ${DELFEM2_INC}/gl_cad_dyntri_v23.h  
  main.cpp  
)

target_link_libraries(${PROJECT_NAME}  
 ${GLUT_LIBRARY} 
 ${OPENGL_LIBRARY}
)
