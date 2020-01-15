set(PROJ_NAME dfm2_cuda_static)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")

find_package(CUDA REQUIRED)
message("Cuda version: " ${CUDA_VERSION})

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(${PROJ_NAME} STATIC
        ${DELFEM2_INC}/cuda/cu_matvec.cu
        ${DELFEM2_INC}/cuda/cu_matvec.h
    )
endif()

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(WIN32)
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND copy
        "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\c_core.*.pyd"
        ${NATIVE_PROJECT_SOURCE_DIR}\\..\\..\\PyDelFEM2\\
        COMMENT "Post Event"
        )
else()
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/libdfm2_cuda_static.*
        ${PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
        )
endif()