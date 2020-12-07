cmake_minimum_required(VERSION 3.12)
project(003_vecdot)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")

# For CUDA
find_package(CUDA REQUIRED)
message("-- CUDA --")
message("Version: " ${CUDA_VERSION})
message("Library: " ${CUDA_CUDA_LIBRARY})
message("Runtime: " ${CUDA_CUDART_LIBRARY})

if(CUDA_FOUND)
  message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_EXECUTABLE(${PROJECT_NAME}
        ${DELFEM2_INC}/cuda/cu_matvec.cu
        ${DELFEM2_INC}/cuda/cu_matvec.h
        main.cpp)
endif()

include_directories(
    ${DELFEM2_INCLUDE_DIR})