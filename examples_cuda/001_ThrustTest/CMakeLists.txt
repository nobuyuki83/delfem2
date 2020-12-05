cmake_minimum_required(VERSION 3.18)
project(001_thrusttest)

# For CUDA
find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    INCLUDE_DIRECTORIES(
            ${CUDA_INCLUDE_DIRS})
    CUDA_ADD_EXECUTABLE(${PROJECT_NAME}
            main.cu
            )
endif()
