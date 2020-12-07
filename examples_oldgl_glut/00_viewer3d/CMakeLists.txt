cmake_minimum_required(VERSION 3.12)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

#########################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
set(CMAKE_CXX_FLAGS "-O2 /W4 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996\
    /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )
endif()

########################

project(00_viewer3d)

set(DELFEM2_INCLUDE_DIR "../../include")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

set(DELFEM2_INC         "../../include/delfem2")
IF(NOT USE_HEADERONLY)
  message("compile as STATIC libary: ${PROJECT_NAME} ${USE_HEADERONLY}")
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  set(SRC_DFM2
      ${DELFEM2_INC}/quat.h                 ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/color.h                ${DELFEM2_INC}/color.cpp
      ${DELFEM2_INC}/mat4.h                 ${DELFEM2_INC}/mat4.cpp
      ${DELFEM2_INC}/quat.h                 ${DELFEM2_INC}/quat.cpp
      
      ${DELFEM2_INC}/cam3_m4q.h             ${DELFEM2_INC}/cam3_m4q.cpp
      
      ${DELFEM2_INC}/opengl/color_glold.h   ${DELFEM2_INC}/opengl/color_glold.cpp
      ${DELFEM2_INC}/opengl/funcs_glold.h   ${DELFEM2_INC}/opengl/funcs_glold.cpp
      )
ELSE()
  message("compile as HEADER-ONLY library: ${PROJECT_NAME} ${USE_HEADERONLY}")
ENDIF()

add_executable(${PROJECT_NAME}
  ${SRC_DFM2}
  main.cpp
)

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY})
