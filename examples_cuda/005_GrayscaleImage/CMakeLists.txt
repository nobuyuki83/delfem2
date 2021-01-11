cmake_minimum_required(VERSION 3.12)
project(005_GrayscaleImage)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")
add_definitions(-DPATH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# For CUDA
find_package(CUDA REQUIRED)

message("Cuda version: " ${CUDA_VERSION})

INCLUDE_DIRECTORIES(
        ${CUDA_INCLUDE_DIRS}
        ../../3rd_party)

CUDA_ADD_EXECUTABLE(${PROJECT_NAME}
        main.cu)

