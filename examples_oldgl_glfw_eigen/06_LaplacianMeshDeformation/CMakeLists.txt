cmake_minimum_required(VERSION 3.12)

#############################
# set C++ detail
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
ENDIF()

#############################
# set project name

project(06_LaplacianMeshDeformation)

# define a macro to get absolute path in C++
add_definitions(-DSOURCE_DIR="${PROJECT_SOURCE_DIR}")

#############################
# specifying libraries to use

# use opengl
find_package(OpenGL REQUIRED)

# use glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../external/glfwlib)
find_package(glfw3 REQUIRED)

# use eigen
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libeigen)
FIND_PACKAGE(Eigen3 REQUIRED)

set(DELFEM2_INCLUDE_DIR "../../include")

########################
# include, build, & link

include_directories(${PROJECT_NAME}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
)

target_link_libraries(${PROJECT_NAME}
  OpenGL::GL
  glfw
  Eigen3::Eigen
)