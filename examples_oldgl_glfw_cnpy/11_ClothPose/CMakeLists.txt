cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -g \
    -Wno-unused-variable \
    -Wno-deprecated-declarations")
ENDIF()

####################################

project(11_ClothPose)

set(TEST_INPUTS "../../test_inputs")
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/${TEST_INPUTS}")
add_definitions(-DPATH_OUTPUT_DIR="${PROJECT_SOURCE_DIR}")

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/GLFW_Lib)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
add_definitions(-DDFM2_HEADER_ONLY=ON)            
set(DELFEM2_INCLUDE_DIR "../../include")

# zlib
find_package(ZLIB REQUIRED)

######################################

set(3RD_PARTY           "../../3rd_party")
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}
  ${TEST_INPUTS}
	${ZLIB_INCLUDE_DIRS}
)

set(DELFEM2_INC         "../../include/delfem2")
add_executable(${PROJECT_NAME}
  ${3RD_PARTY}/cnpy/cnpy.h                  
  ${3RD_PARTY}/cnpy/cnpy.cpp
  ${TEST_INPUTS}/inputs_imgboneloc.h
  main.cpp
)

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
    ${ZLIB_LIBRARY}
)