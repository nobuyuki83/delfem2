cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
project(testsDelFEM2)

set(DELFEM2_INCLUDE_DIR "../src_cpp/include")
set(DELFEM2_INC         "../src_cpp/include/delfem2")
set(DELFEM2_SRC         "../src_cpp/cpp")

enable_testing()

set(GTEST_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/lib/libgtest.a)
set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/lib/libgtest_main.a)
set(GTEST_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/googletest/googletest/include)

set(MY_BINARY_NAME runUnitTests)
message(STATUS "${GTEST_INCLUDE_DIR}")

include_directories(
  ${GTEST_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

# Add test cpp file
add_executable(${MY_BINARY_NAME}
  main.cpp
  ${DELFEM2_SRC}/vec3.cpp      ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/mat3.cpp      ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/quat.cpp      ${DELFEM2_INC}/quat.h
  ${DELFEM2_SRC}/mathfuncs.cpp ${DELFEM2_INC}/mathfuncs.h
  ${DELFEM2_SRC}/msh.cpp       ${DELFEM2_INC}/msh.h
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