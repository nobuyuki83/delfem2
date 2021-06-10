cmake_minimum_required(VERSION 3.12)

#########################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )

########################

project(149_PathTracingSimple)

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/GLFW_Lib)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
add_definitions(-DDFM2_HEADER_ONLY=ON)      
set(DELFEM2_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../include")

########################


add_executable(${PROJECT_NAME}
  main.cpp
)

include_directories(
  ${DELFEM2_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
)
