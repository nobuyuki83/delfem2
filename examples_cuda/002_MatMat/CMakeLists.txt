cmake_minimum_required(VERSION 3.18)
project(002_matmat)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")

# For CUDA
find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_EXECUTABLE(${PROJECT_NAME}
        ${DELFEM2_INC}/cuda/cu_matvec.cu
        ${DELFEM2_INC}/cuda/cu_matvec.h
        main.cpp
    )
    include_directories(${DELFEM2_INCLUDE_DIR})
endif()
