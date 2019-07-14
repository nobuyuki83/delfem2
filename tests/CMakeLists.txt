cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
project(testGeoSim)

# Setup testing
enable_testing()

set(DELFEM2_INCLUDE_DIR "../src_cpp/include")
set(DELFEM2_INC         "../src_cpp/include/delfem2")
set(DELFEM2_SRC         "../src_cpp/cpp")

set(GTEST_INCLUDE_DIR  "googletest/googletest/include")

include_directories(
  ${GTEST_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR})

# Add test cpp file
add_executable(runUnitTests
  main.cpp
  ${DELFEM2_SRC}/vec3.cpp    ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/mat3.cpp    ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/quat.cpp    ${DELFEM2_INC}/quat.h
)

find_library(GTEST_LIBRARY
  NAMES libgtest.a
  PATHS ${CMAKE_CURRENT_SOURCE_DIR}/googletest/lib)  

find_library(GTEST_MAIN_LIBRARY
  NAMES libgtest_main.a
  PATHS ${CMAKE_CURRENT_SOURCE_DIR}/googletest/lib)  

# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests 
  ${GTEST_LIBRARY}
  ${GTEST_MAIN_LIBRARY}
  -pthread)

add_test(
  NAME runUnitTests
  COMMAND runUnitTests
)
