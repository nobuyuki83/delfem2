cmake_minimum_required(VERSION 3.0)

# compile as "cmake -DEMSCRIPTEN=1 ."  to build with Emscripten
option(EMSCRIPTEN "Use Emscripten:" OFF)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

##########################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(NOT EMSCRIPTEN)
  IF(MSVC)
  ELSE()
    set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
  ENDIF()
ENDIF()

##########################################

project(00_openwin)

##########################################

# glfw
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

#########################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  #/Users/nobuyuki/emsdk/fastcomp/emscripten/system/include/compat
)

add_executable(${PROJECT_NAME}
  main.cpp
)

target_link_libraries(${PROJECT_NAME} 
  ${OPENGL_LIBRARY}
  ${GLFW_LIBRARIES}
)
