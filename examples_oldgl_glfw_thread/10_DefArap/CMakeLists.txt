cmake_minimum_required(VERSION 3.12)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")
ENDIF ()

##########################################

project(00_DefArap)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
add_definitions(-DDFM2_HEADER_ONLY=ON)
set(DELFEM2_INCLUDE_DIR "../../include")

#######################################

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
    )

add_executable(${PROJECT_NAME}
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
    )