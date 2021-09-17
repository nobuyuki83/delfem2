cmake_minimum_required(VERSION 3.12)

option(USE_STATIC_LIB "compile delfem2 with static library?" OFF)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
ENDIF()

######################################################

project(000_OpenWin)

# glfw
set(CMAKE_PREFIX_PATH ../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

# opengl
find_package(OpenGL REQUIRED)

######################################################

include_directories(
  ${DELFEM2_INCLUDE_DIR}
  )

add_executable(${PROJECT_NAME}
  main.cpp
  )

target_link_libraries(${PROJECT_NAME}
  OpenGL::GL
  glfw
  )