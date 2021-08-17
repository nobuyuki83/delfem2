cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if (MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
else ()
  set(CMAKE_CXX_FLAGS "-Wall -g")
endif ()

############################################

project(23_SphereTracing)

add_definitions(-DPATH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# glad
get_filename_component(GLAD_INCLUDE_DIR
		${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad3/include ABSOLUTE)

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# delfem2
set(DELFEM2_INCLUDE_DIR "../../include")

############################################

include_directories(
		${OPENGL_INCLUDE_DIR}
		${GLAD_INCLUDE_DIR}
		${GLFW_INCLUDE_DIR}
		${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
		${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad3/include/glad/glad.h
		${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad3/src/glad.c
		main.cpp
		)

target_link_libraries(${PROJECT_NAME}
		${DFM2_LIBRARY}
		${OPENGL_LIBRARY}
		${GLFW_LIBRARIES}
		)
