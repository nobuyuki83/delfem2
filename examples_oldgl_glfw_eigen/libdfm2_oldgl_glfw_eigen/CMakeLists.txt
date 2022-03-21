cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
ENDIF()

####################################

set(PROJ_NAME dfm2_oldgl_glfw_eigen)
project(${PROJ_NAME})

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)
get_target_property(GLFW_INCLUDE_DIR glfw INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "glfw include dir: ${GLFW_INCLUDE_DIR}")

# glad
get_filename_component(GLAD_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad3/include ABSOLUTE)

# opengl
find_package(OpenGL REQUIRED)

# eigen
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libeigen)
FIND_PACKAGE(Eigen3 REQUIRED)

# delfem2
add_definitions(-DDFM2_STATIC_LIBRARY=ON)
set(DELFEM2_INCLUDE_DIR "../../include")

#####################################

include_directories(
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party
)

file(GLOB 
    DFM2_SRC
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/*.cpp
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/*.h)
file(GLOB 
    DFM2_OPENGL_SRC
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/opengl/*.cpp
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/opengl/*.h)        
file(GLOB 
    DFM2_OPENGL_OLD_SRC
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/opengl/old/*.cpp
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/opengl/old/*.h)
file(GLOB 
    DFM2_GLFW_SRC
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/glfw/*.cpp
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/glfw/*.h)    
file(GLOB
    DFM2_EIGEN_SRC
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/eigen/*.cpp
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/eigen/*.h)    
set(DELFEM2_INC         "../../include/delfem2")    
add_library(${PROJ_NAME} STATIC
    ${DFM2_SRC}
    ${DFM2_OPENGL_SRC}
    ${DFM2_OPENGL_OLD_SRC}
    ${DFM2_GLFW_SRC}
    ${DFM2_EIGEN_SRC}
    ${DELFEM2_INC}/stb_opengl/glyph.h
    ${DELFEM2_INC}/stb_opengl/glyph.cpp
)

target_link_libraries (${PROJ_NAME}
    Eigen3::Eigen)


file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(WIN32)
    message("copy lib from: ${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\dfm2_static_glfwold.lib")
    message("copy lib to: ${NATIVE_PROJECT_SOURCE_DIR}")
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND copy
        "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\dfm2_static_glfwold.lib"
        ${NATIVE_PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
        )
else()
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/libdfm2_oldgl_glfw_eigen.*
        ${PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
        )
endif()
