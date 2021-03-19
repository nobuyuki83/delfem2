cmake_minimum_required(VERSION 3.12)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

###################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

####################################

project(111_ExponentialMapElem)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/GLFW_Lib)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
if(USE_HEADERONLY)
  add_definitions(-DDFM2_HEADER_ONLY=ON)  
else()
  if(MSVC)
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/dfm2_static_glfwold.lib)
  else()
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/libdfm2_static_glfwold.a)
  endif()
endif()
set(DELFEM2_INCLUDE_DIR "../../include")

#########################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${DFM2_LIBRARY}
  ${OPENGL_LIBRARY}
  ${GLFW_LIBRARIES}
)
