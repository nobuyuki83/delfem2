cmake_minimum_required(VERSION 3.12)

option(DFM2_HEADER_ONLY "compile delfem2 as a header-only library?" ON)

####################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS
  "/W4 -O2 -std=c++11 \
  /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 \
  /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 \
  /wd4505 /wd4456 /wd4055")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

####################################

project(41_fem2d_poisson)

# glad
get_filename_component(GLAD_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/../../3rd_party/glad4/include ABSOLUTE)

# dfm2
if(DFM2_HEADER_ONLY)
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
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

#############################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

set(3RD_PARTY           "../../3rd_party")
add_executable(${PROJECT_NAME}
  ${3RD_PARTY}/glad4/include/glad/glad.h   
  ${3RD_PARTY}/glad4/src/glad.c
  main.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${DFM2_LIBRARY}
  ${OPENGL_LIBRARY}
  ${GLFW_LIBRARIES}
)