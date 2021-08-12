cmake_minimum_required(VERSION 3.12)

####################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "-Wall -g")
ENDIF ()

####################################

project(00_EigenModes)

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# eigen
FILE(GLOB EIGEN3_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/eigen*) # find downloaded eigen directories in 3rd_party
SET(EIGEN3_ROOT
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libeigen/include # installed eigen is higher priority
    ${EIGEN3_ROOT}) # specify candidate for eigen location
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindEigen3.cmake) # call script to find eigen3
#find_package(Eigen3 REQUIRED) # this works only for installation from package manager
message(STATUS "eigen3 locatoin: ${EIGEN3_INCLUDE_DIR}")

# dfm2
set(DELFEM2_INCLUDE_DIR "../../include")

#######################################################

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
    )
