cmake_minimum_required(VERSION 3.12)

option(USE_STATIC_LIB "compile delfem2 as a static library" OFF)

################################

set(CMAKE_CXX_STANDARD 17)
if(MSVC)
  # GTest is static compiled by default
  # GTest requires multi-threading
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

################################

if(USE_STATIC_LIB)
  project(test_static)
ELSE()
  project(test_hdronly)
ENDIF()
set(MY_BINARY_NAME runUnitTests)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../test_inputs")

set(DELFEM2_INCLUDE_DIR "../include")

# google test
set(GTEST_ROOT ${PROJECT_SOURCE_DIR}/../3rd_party/libgtest)
find_package( GTest REQUIRED )

# thread
find_package(Threads REQUIRED)

enable_testing()

include_directories(
  ${GTEST_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ../3rd_party
)

# Add test cpp file
IF(USE_STATIC_LIB)
  message("compile as STATIC libary: ${PROJECT_NAME} USE_STATIC_LIB=${USE_STATIC_LIB}")
  add_definitions(-DDFM2_STATIC_LIBRARY)
  file(GLOB 
    DFM2_SRC
    ${PROJECT_SOURCE_DIR}/../include/delfem2/*.cpp
    ${PROJECT_SOURCE_DIR}/../include/delfem2/*.h)
ELSE()
  message("compile delfem2 as HEADER-ONLY library: ${PROJECT_NAME} ${USE_STATIC_LIB}}")
ENDIF()

file(GLOB
    TEST_SRC
    *.cpp)
add_executable(${MY_BINARY_NAME}
  ${DFM2_SRC}
  ../include/delfem2/tinygltf/io_gltf.h
  ../include/delfem2/tinygltf/io_gltf.cpp
  ${TEST_SRC})

target_link_libraries(${MY_BINARY_NAME}
  ${GTEST_LIBRARY}
  ${GTEST_MAIN_LIBRARY}
  Threads::Threads
)

add_test(
  NAME ${MY_BINARY_NAME}
  COMMAND ${MY_BINARY_NAME}
)
