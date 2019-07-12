project(mathexpeval)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_SRC}/mathexpeval.cpp          ${DELFEM2_INC}/mathexpeval.h
  main.cpp
)

target_link_libraries(${PROJECT_NAME} 
)