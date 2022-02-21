cmake_minimum_required(VERSION 3.10)
project(VDB_SAMPLE)

# build option
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O1 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -DDELEASE")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -Og -DNDEBUG")
message("Build type: " ${CMAKE_BUILD_TYPE})

# build path
message("Current directory: " ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})

# option
option(BUILD_AS_SHARED "whether to build as shared lib" OFF)

list(APPEND SAMPLE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB SAMPLE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# tbb
find_library(TBB_LIBRARY tbb HINTS /usr/local/lib)

# openexr
#find_library(HALF_LIB half HINTS /usr/local/lib)

# openvdb
list(APPEND SAMPLE_INCLUDE_DIR /usr/local/include)
list(APPEND SAMPLE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../3rd_party/openvdblib/include)
find_library(OpenVDB_LIBRARY
    openvdb HINTS ${PROJECT_SOURCE_DIR}/../../3rd_party/openvdblib/lib)

## include headers
include_directories(
    ${SAMPLE_INCLUDE_DIR})
add_executable(sample
    ${SAMPLE_SOURCES})
target_link_libraries(
    sample
    ${TBB_LIBRARY}
    #${HALF_LIB}
    ${OpenVDB_LIBRARY}
)
