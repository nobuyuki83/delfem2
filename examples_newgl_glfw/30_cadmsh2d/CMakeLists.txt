cmake_minimum_required(VERSION 3.12)

option(USE_STATIC_LIB "compile delfem2 as a static library?" OFF)

#########################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")
ENDIF()

#########################################

project(30_cadmsh2d)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
IF (USE_STATIC_LIB)
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  if(MSVC)
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/dfm2_static_glfwnew.lib)
  else()
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/libdfm2_static_glfwnew.a )
  endif()
endif()
set(DELFEM2_INCLUDE_DIR "../../include")

# glad
set(GLAD_INCLUDE_DIR "../../3rd_party/glad4/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

################################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

set(3RD_PARTY "../../3rd_party")
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