cmake_minimum_required(VERSION 3.12)

option(DFM2_HEADER_ONLY "compile delfem2 as a header only library" ON)

################################

set(CMAKE_CXX_STANDARD 17)
if(MSVC)
  # GTest is static compiled by default
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")  
  set(CMAKE_CXX_FLAGS "/W4")
endif()

################################

project(test_eigen)
set(MY_BINARY_NAME runUnitTests)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../test_inputs")

################################

# google test
set(GTEST_ROOT ${PROJECT_SOURCE_DIR}/../3rd_party/libgtest)
find_package( GTest REQUIRED )

# thread
find_package(Threads REQUIRED)

# eigen
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../3rd_party/libeigen)
FIND_PACKAGE(Eigen3 REQUIRED)

# delfem2
set(DELFEM2_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/../include")

#############################

enable_testing()

include_directories(
  ${GTEST_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

file(GLOB
    TEST_SRC
    *.cpp)
add_executable(${MY_BINARY_NAME}
  ${TEST_SRC}
)

target_link_libraries(${MY_BINARY_NAME}
  ${GTEST_LIBRARY}
  ${GTEST_MAIN_LIBRARY}
  Threads::Threads
  Eigen3::Eigen
)

add_test(
  NAME ${MY_BINARY_NAME}
  COMMAND ${MY_BINARY_NAME}
)
