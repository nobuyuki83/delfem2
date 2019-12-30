cmake_minimum_required(VERSION 2.8)
project(000_simple)

# For CUDA
find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_EXECUTABLE(${PROJECT_NAME}
        main.cu
    )
endif()
