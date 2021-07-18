cmake_minimum_required(VERSION 3.12)

option(DFM2_HEADER_ONLY "compile delfem2 as a header-only library?" ON)

######################################################

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "-Wall -g")
ENDIF ()

###########################################

project(02_FemSolidLinear2)


add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# oepngl
find_package(OpenGL REQUIRED)

# eigen
# find downloaded eigen directories in 3rd_party
FILE(GLOB EIGEN3_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/eigen*) 
# specify candidate for eigen location
SET(EIGEN3_ROOT
  ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/Eigen_Lib/include # installed eigen is higher priority
  ${EIGEN3_ROOT} ) 
# call script to find eigen3  
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindEigen3.cmake) 
#find_package(Eigen3 REQUIRED) # this works only for installation from package manager
message(STATUS "eigen3 locatoin: ${EIGEN3_INCLUDE_DIR}")

# dfm2
IF (DFM2_HEADER_ONLY)
  add_definitions(-DDFM2_HEADER_ONLY=ON)
else ()
  if (MSVC)
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/dfm2_static_glfwold.lib)
  else ()
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/libdfm2_static_glfwold.a)
  endif ()
ENDIF ()
set(DELFEM2_INCLUDE_DIR "../../include")

############################################

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
    )