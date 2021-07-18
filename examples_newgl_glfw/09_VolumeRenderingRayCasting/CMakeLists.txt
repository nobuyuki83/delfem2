cmake_minimum_required(VERSION 3.12)

#######################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if (MSVC)
	set(CMAKE_CXX_FLAGS "-O2 /W4 \
     /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 \
     /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456")
else ()
	set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
endif ()

############################################

project(09_VolRenRayCast)

add_definitions(-DPATH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# glad
get_filename_component(GLAD_INCLUDE_DIR
		${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/include ABSOLUTE)
SET(GLAD_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/include/glad/glad.h
		${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/src/glad.c )

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
add_definitions(-DDFM2_HEADER_ONLY=ON)
set(DELFEM2_INCLUDE_DIR "../../include")

############################################

include_directories(
		${OPENGL_INCLUDE_DIR}
		${GLAD_INCLUDE_DIR}
		${GLFW_INCLUDE_DIR}
		${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
		${GLAD_SOURCES}
		main.cpp
		)

target_link_libraries(${PROJECT_NAME}
		${DFM2_LIBRARY}
		${OPENGL_LIBRARY}
		${GLFW_LIBRARIES}
		)
