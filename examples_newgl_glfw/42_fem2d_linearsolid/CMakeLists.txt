cmake_minimum_required(VERSION 3.12)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

##############################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "-O2 -std=c++11 /W4 \
      /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 \
      /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 \
      /wd4838 /wd4305 /wd4505 /wd4456 /wd4055")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

##############################################

project(42_fem2d_linearsolid)

# glad
set(GLAD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

# dfm2
if(USE_HEADERONLY)
  add_definitions(-DDFM2_HEADER_ONLY=ON)
else()
  if(MSVC)
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/dfm2_static_glfwnew.lib)
  else()
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/libdfm2_static_glfwnew.a )
  endif()
endif()
set(DELFEM2_INCLUDE_DIR "../../include")

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/GLFW_Lib)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

###########################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

set(3RD_PARTY  "../../3rd_party")
add_executable(${PROJECT_NAME}
  ${3RD_PARTY}/glad4/include/glad/glad.h   
  ${3RD_PARTY}/glad4/src/glad.c
  ${SRC_DFM2}
  main.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${DFM2_LIBRARY}
  ${OPENGL_LIBRARY}
  ${GLFW_LIBRARIES}
)