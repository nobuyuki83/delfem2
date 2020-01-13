cmake_minimum_required(VERSION 2.8)
project(003_vecdot)

# For CUDA
find_package(CUDA REQUIRED)
message("-- CUDA --")
message("Version: " ${CUDA_VERSION})
message("Library: " ${CUDA_CUDA_LIBRARY})
message("Runtime: " ${CUDA_CUDART_LIBRARY})

if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_EXECUTABLE(${PROJECT_NAME}
        main.cu
    )
endif()
