project(edit_cubes)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
else()
set( CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g -std=c++11" )
endif()

add_definitions(-DUSE_GL)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h                  ${DELFEM2_INC}/vec2.cpp       
  ${DELFEM2_INC}/vec3.h                  ${DELFEM2_INC}/vec3.cpp      
  ${DELFEM2_INC}/mat3.h                  ${DELFEM2_INC}/mat3.cpp      
  ${DELFEM2_INC}/quat.h                  ${DELFEM2_INC}/quat.cpp      
  ${DELFEM2_INC}/camera.h                ${DELFEM2_INC}/camera.cpp
  ${DELFEM2_INC}/color.h                 ${DELFEM2_INC}/color.cpp  
  ${DELFEM2_INC}/voxel.h                 ${DELFEM2_INC}/voxel.cpp      

  ${DELFEM2_INC}/opengl/gl_voxbv.h       ${DELFEM2_INC}/opengl/gl_voxbv.cpp
  ${DELFEM2_INC}/opengl/gl2_v23.h        ${DELFEM2_INC}/opengl/gl2_v23.cpp
  ${DELFEM2_INC}/opengl/gl2_funcs.h      ${DELFEM2_INC}/opengl/gl2_funcs.cpp    
  ${DELFEM2_INC}/opengl/gl2_color.h      ${DELFEM2_INC}/opengl/gl2_color.cpp

  main.cpp
)

target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY}
)
