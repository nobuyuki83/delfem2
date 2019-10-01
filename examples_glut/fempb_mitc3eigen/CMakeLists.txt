project(fempb_mitc3eigen)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )
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
  ${DELFEM2_INC}/vec2.h         ${DELFEM2_SRC}/vec2.cpp
  ${DELFEM2_INC}/msh.h          ${DELFEM2_SRC}/msh.cpp
  ${DELFEM2_INC}/mshtopo.h      ${DELFEM2_SRC}/mshtopo.cpp
  ${DELFEM2_INC}/emat.h         ${DELFEM2_SRC}/emat.cpp
  ${DELFEM2_INC}/dtri.h         ${DELFEM2_SRC}/dtri.cpp
  ${DELFEM2_INC}/mats.h         ${DELFEM2_SRC}/mats.cpp
  ${DELFEM2_INC}/camera.h       ${DELFEM2_SRC}/camera.cpp

  ${DELFEM2_INC}/dtri_v2.h      ${DELFEM2_SRC}/dtri_v2.cpp
  ${DELFEM2_INC}/ilu_mats.h     ${DELFEM2_SRC}/ilu_mats.cpp
  ${DELFEM2_INC}/fem_emats.h    ${DELFEM2_SRC}/fem_emats.cpp

  ${DELFEM2_INC}/gl2_funcs.h    ${DELFEM2_SRC}/gl2_funcs.cpp
  ${DELFEM2_INC}/gl_color.h     ${DELFEM2_SRC}/gl_color.cpp                     
  
  main.cpp
)

target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY})
