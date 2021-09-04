cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -g")
ENDIF()

####################################

set(PROJ_NAME dfm2_static_glfwnew)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../../include")
set(3RD_PARTY           "../../3rd_party")
add_definitions(-DDFM2_STATIC_LIBRARY=ON)

# glfw
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)


####################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}/glad4/include
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
    DFM2_OPENGL_NEW_SRC
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/opengl/new/*.cpp
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/opengl/new/*.h)            
file(GLOB 
    DFM2_GLFW_SRC
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/glfw/*.cpp
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/glfw/*.h)        
set(DELFEM2_INC         "../../include/delfem2")
add_library(${PROJ_NAME} STATIC
    ${DFM2_SRC}
    ${DFM2_OPENGL_SRC}
    ${DFM2_OPENGL_NEW_SRC}    
    ${DFM2_GLFW_SRC}        
)

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(WIN32)
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND copy
        ${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\dfm2_static_glfwnew.lib
        ${NATIVE_PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
    )
else()
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/libdfm2_static_glfwnew.*
        ${PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
    )
endif()
