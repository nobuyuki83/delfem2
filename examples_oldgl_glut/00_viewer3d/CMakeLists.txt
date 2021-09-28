cmake_minimum_required(VERSION 3.12)

#########################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
endif()

########################

project(00_viewer3d)

set(DELFEM2_INCLUDE_DIR "../../include")

find_package(OpenGL REQUIRED)

find_package(GLUT REQUIRED)

add_executable(${PROJECT_NAME}
  main.cpp
)

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME}
  GLUT::GLUT
  OpenGL::GL
  OpenGL::GLU
  )
