cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
project(testsDelFEM2)

if(MSVC)
    # GTest is static compiled by default
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
endif()

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../test_inputs")

set(DELFEM2_INCLUDE_DIR "../src_cpp/include")
set(DELFEM2_INC         "../src_cpp/include/delfem2")
set(DELFEM2_SRC         "../src_cpp/cpp")

enable_testing()

if(MSVC)
  if(CMAKE_CL_64)
    set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS64/lib/Release/gtest_main.lib)
    set(GTEST_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS64/lib/Release/gtest.lib)
  else()
    set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS32/lib/Release/gtest_main.lib)
    set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS32/lib/Release/gtest.lib)
  endif()
else()
  set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/lib/libgtest_main.a)
  set(GTEST_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/lib/libgtest.a)
endif()
set(GTEST_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/googletest/googletest/include)

set(MY_BINARY_NAME runUnitTests)
message(STATUS "${GTEST_INCLUDE_DIR}")

include_directories(
  ${GTEST_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

# Add test cpp file
add_executable(${MY_BINARY_NAME}
  ${DELFEM2_INC}/vec3.h        ${DELFEM2_SRC}/vec3.cpp
  ${DELFEM2_INC}/funcs.h       ${DELFEM2_SRC}/funcs.cpp
  ${DELFEM2_INC}/vec2.h        ${DELFEM2_SRC}/vec2.cpp
  ${DELFEM2_INC}/mat3.h        ${DELFEM2_SRC}/mat3.cpp
  ${DELFEM2_INC}/quat.h        ${DELFEM2_SRC}/quat.cpp
  ${DELFEM2_INC}/mathfuncs.h   ${DELFEM2_SRC}/mathfuncs.cpp
  ${DELFEM2_INC}/msh.h         ${DELFEM2_SRC}/msh.cpp
  ${DELFEM2_INC}/mshtopo.h     ${DELFEM2_SRC}/mshtopo.cpp
  ${DELFEM2_INC}/mshio.h       ${DELFEM2_SRC}/mshio.cpp
  ${DELFEM2_INC}/voxel.h       ${DELFEM2_SRC}/voxel.cpp
  ${DELFEM2_INC}/mathexpeval.h ${DELFEM2_SRC}/mathexpeval.cpp
  ${DELFEM2_INC}/lp.h
  main.cpp  
)

target_link_libraries(${MY_BINARY_NAME}
  ${GTEST_LIBRARY}
  ${GTEST_MAIN_LIBRARY}
  -pthread
)

add_test(
  NAME ${MY_BINARY_NAME}
  COMMAND ${MY_BINARY_NAME}
)
