project(pbdcloth_cad_gltf)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
set(CMAKE_CXX_FLAGS "-O2 /W4\
  /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 \
  /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 ")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations  -std=c++11")
endif()

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(DELFEM2_EXT         "../../src_cpp/external")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h                   ${DELFEM2_INC}/vec2.cpp                 
  ${DELFEM2_INC}/vec3.h                   ${DELFEM2_INC}/vec3.cpp                 
  ${DELFEM2_INC}/mat3.h                   ${DELFEM2_INC}/mat3.cpp                 
  ${DELFEM2_INC}/quat.h                   ${DELFEM2_INC}/quat.cpp 
  ${DELFEM2_INC}/camera.h                 ${DELFEM2_INC}/camera.cpp                  
  ${DELFEM2_INC}/dtri.h                   ${DELFEM2_INC}/dtri.cpp
  ${DELFEM2_INC}/mshtopo.h                ${DELFEM2_INC}/mshtopo.cpp
  ${DELFEM2_INC}/mshmisc.h                ${DELFEM2_INC}/mshmisc.cpp
  ${DELFEM2_INC}/bvh.h                    ${DELFEM2_INC}/bvh.cpp
  ${DELFEM2_INC}/funcs.h                  ${DELFEM2_INC}/funcs.cpp
  ${DELFEM2_INC}/bv.h        
      
  ${DELFEM2_INC}/dtri_v2.h                ${DELFEM2_INC}/dtri_v2.cpp
  ${DELFEM2_INC}/objfunc_v23.h            ${DELFEM2_INC}/objfunc_v23.cpp
  ${DELFEM2_INC}/objfunc_v23dtri.h        ${DELFEM2_INC}/objfunc_v23dtri.cpp
  ${DELFEM2_INC}/cad2d.h                  ${DELFEM2_INC}/cad2d.cpp
  ${DELFEM2_INC}/rig_v3q.h                ${DELFEM2_INC}/rig_v3q.cpp
  ${DELFEM2_INC}/v23m3q.h                 ${DELFEM2_INC}/v23m3q.cpp
  ${DELFEM2_INC}/srchuni_v3.h             ${DELFEM2_INC}/srchuni_v3.cpp
  ${DELFEM2_INC}/srch_v3bvhmshtopo.h
  
  ${DELFEM2_INC}/opengl/gl2_funcs.h       ${DELFEM2_INC}/opengl/gl2_funcs.cpp
  ${DELFEM2_INC}/opengl/gl2_v23.h         ${DELFEM2_INC}/opengl/gl2_v23.cpp
  ${DELFEM2_INC}/opengl/gl2_v23dtricad.h  ${DELFEM2_INC}/opengl/gl2_v23dtricad.cpp
  ${DELFEM2_INC}/opengl/gl_rig_v23q.h     ${DELFEM2_INC}/opengl/gl_rig_v23q.cpp
  
  ${DELFEM2_EXT}/io_gltf.h
  ${DELFEM2_EXT}/io_gltf.cpp

  main.cpp  
)

target_link_libraries(${PROJECT_NAME}  
 ${GLUT_LIBRARY} 
 ${OPENGL_LIBRARY}
)

