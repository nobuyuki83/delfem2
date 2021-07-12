cmake_minimum_required(VERSION 3.12)

#########################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800  /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE ()
  set(CMAKE_CXX_FLAGS "-Wall  -g")
ENDIF ()

########################

project(633_mpm2snow)

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