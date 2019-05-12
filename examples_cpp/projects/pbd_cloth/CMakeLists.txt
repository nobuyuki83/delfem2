project(pbd_cloth)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")

set(DELFEM2_INCLUDE_DIR "../../../src_cpp/include")
set(DELFEM2_SRC         "../../../src_cpp/cpp")
set(DELFEM2_INC         "../../../src_cpp/include/delfem2")

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
  ${DELFEM2_SRC}/mat3.cpp           ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/quat.cpp           ${DELFEM2_INC}/quat.h
  ${DELFEM2_SRC}/dyntri.cpp         ${DELFEM2_INC}/dyntri.h
  ${DELFEM2_SRC}/mshtopo.cpp        ${DELFEM2_INC}/mshtopo.h

  ${DELFEM2_SRC}/v23m3q.cpp         ${DELFEM2_INC}/v23m3q.h
  ${DELFEM2_SRC}/dyntri_v2.cpp      ${DELFEM2_INC}/dyntri_v2.h
  ${DELFEM2_INC}/cloth_internal.h

  ${DELFEM2_SRC}/funcs_gl.cpp       ${DELFEM2_INC}/funcs_gl.h
  ${DELFEM2_SRC}/camera_gl.cpp      ${DELFEM2_INC}/camera_gl.h

  main.cpp  
)

target_link_libraries(${PROJECT_NAME}  
 ${GLUT_LIBRARY} 
 ${OPENGL_LIBRARY}
)

