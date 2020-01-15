cmake_minimum_required(VERSION 2.8)

set(PROJ_NAME test_cuda)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../include")
set(DELFEM2_INC         "../include/delfem2")
set(3RD_PARTY           "../3rd_party")

get_filename_component(GTEST_DIR ${3RD_PARTY}/googletest ABSOLUTE)

enable_testing()

add_executable(${PROJ_NAME}
    main.cpp
)

include_directories(
    ${GTEST_DIR}/googletest/include
    ${DELFEM2_INCLUDE_DIR})

# Add sub project
add_subdirectory(dfm2_cuda_static)
target_include_directories(${PROJ_NAME} PUBLIC
    dfm2_cuda_static)
target_link_libraries(${PROJ_NAME}
    dfm2_cuda_static
    ${GTEST_DIR}/lib/libgtest.a
    ${GTEST_DIR}/lib/libgtest_main.a
    -pthread)

add_test(
    NAME ${PROJ_NAME}
    COMMAND ${PROJ_NAME}
)