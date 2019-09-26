project(pbdcloth_cad)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations")
endif()

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
  ${DELFEM2_INC}/funcs.h         ${DELFEM2_SRC}/funcs.cpp
  ${DELFEM2_INC}/vec2.h             ${DELFEM2_SRC}/vec2.cpp                 
  ${DELFEM2_INC}/vec3.h             ${DELFEM2_SRC}/vec3.cpp                 
  ${DELFEM2_INC}/mat3.h             ${DELFEM2_SRC}/mat3.cpp                 
  ${DELFEM2_INC}/quat.h             ${DELFEM2_SRC}/quat.cpp
  ${DELFEM2_INC}/dtri.h           ${DELFEM2_SRC}/dtri.cpp                 
  ${DELFEM2_INC}/mshtopo.h          ${DELFEM2_SRC}/mshtopo.cpp
  ${DELFEM2_INC}/msh.h              ${DELFEM2_SRC}/msh.cpp
  ${DELFEM2_INC}/bv.h
  ${DELFEM2_INC}/bvh.h              ${DELFEM2_SRC}/bvh.cpp
  ${DELFEM2_INC}/primitive.h        ${DELFEM2_SRC}/primitive.cpp

  ${DELFEM2_INC}/v23m3q.h           ${DELFEM2_SRC}/v23m3q.cpp 
  ${DELFEM2_INC}/objfunc_v23.h      ${DELFEM2_SRC}/objfunc_v23.cpp
  ${DELFEM2_INC}/objfunc_v23dtri.h  ${DELFEM2_SRC}/objfunc_v23dtri.cpp
  ${DELFEM2_INC}/dtri_v2.h        ${DELFEM2_SRC}/dtri_v2.cpp
  ${DELFEM2_INC}/cad2d.h            ${DELFEM2_SRC}/cad2d.cpp
  ${DELFEM2_INC}/srchuni_v3.h       ${DELFEM2_SRC}/srchuni_v3.cpp
  ${DELFEM2_INC}/srch_v3bvhmshtopo.h

  ${DELFEM2_INC}/gl2_funcs.h         ${DELFEM2_SRC}/gl2_funcs.cpp                 
  ${DELFEM2_INC}/gl24_camera.h        ${DELFEM2_SRC}/gl24_camera.cpp
  ${DELFEM2_INC}/gl_v23.h              ${DELFEM2_SRC}/gl_v23.cpp
  ${DELFEM2_SRC}/gl_cad_dyntri_v23.cpp       
  ${DELFEM2_INC}/gl_cad_dyntri_v23.h  

  main.cpp  
)

target_link_libraries(${PROJECT_NAME}  
 ${GLUT_LIBRARY} 
 ${OPENGL_LIBRARY}
)

