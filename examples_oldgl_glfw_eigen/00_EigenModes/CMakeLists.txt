cmake_minimum_required(VERSION 3.12)

####################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
ENDIF ()

####################################

project(00_EigenModes)

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

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
    ${EIGEN3_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    OpenGL::GL
    glfw
    )
