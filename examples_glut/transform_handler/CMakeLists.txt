cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

#########################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
set(CMAKE_CXX_FLAGS "/W4  -O2 \
  /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4456 \
  /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
endif()

#########################

project(trans_hndlr)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

IF(NOT USE_HEADERONLY)
  message("compile as STATIC libary: ${PROJECT_NAME} ${USE_HEADERONLY}")
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  set(SRC_DFM2
      ${DELFEM2_INC}/vec2.h                    ${DELFEM2_INC}/vec2.cpp
      ${DELFEM2_INC}/vec3.h                    ${DELFEM2_INC}/vec3.cpp
      ${DELFEM2_INC}/mat3.h                    ${DELFEM2_INC}/mat3.cpp
      ${DELFEM2_INC}/mat4.h                    ${DELFEM2_INC}/mat4.cpp
      ${DELFEM2_INC}/quat.h                      ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/color.h                   ${DELFEM2_INC}/color.cpp
      ${DELFEM2_INC}/camera.h               ${DELFEM2_INC}/camera.cpp

      ${DELFEM2_INC}/v23m34q.h                  ${DELFEM2_INC}/v23m34q.cpp
      )
ELSE()
  message("compile as HEADER-ONLY library: ${PROJECT_NAME} ${USE_HEADERONLY}")
ENDIF()

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/opengl/glold_funcs.h        ${DELFEM2_INC}/opengl/glold_funcs.cpp
  ${DELFEM2_INC}/opengl/glold_color.h        ${DELFEM2_INC}/opengl/glold_color.cpp   
  ${DELFEM2_INC}/opengl/v3q_glold.h          ${DELFEM2_INC}/opengl/v3q_glold.cpp
  ${SRC_DFM2}
  main.cpp
)

target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY})
