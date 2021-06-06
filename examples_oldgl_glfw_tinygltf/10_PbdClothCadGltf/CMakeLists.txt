cmake_minimum_required(VERSION 3.12)

######################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
set(CMAKE_CXX_FLAGS "-O2 /W4\
  /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 \
  /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 ")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations  -std=c++11")
endif()

####################################

project(10_PbdClothCadGltf)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

####################################

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/GLFW_Lib)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
add_definitions(-DDFM2_HEADER_ONLY=ON)            
set(DELFEM2_INCLUDE_DIR "../../include")

#####################################

set(3RD_PARTY           "../../3rd_party")
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}  
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}  
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_INCLUDE_DIR}/delfem2/tinygltf/io_gltf.h
  ${DELFEM2_INCLUDE_DIR}/delfem2/tinygltf/io_gltf.cpp
  main.cpp  
)

target_link_libraries(${PROJECT_NAME}
  ${OPENGL_LIBRARY}
  ${GLFW_LIBRARIES}
)