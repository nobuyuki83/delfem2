cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")
ENDIF()

####################################

set(PROJ_NAME dfm2_static_glfwold)
project(${PROJ_NAME})

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# glad
get_filename_component(GLAD_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad3/include ABSOLUTE)

# opengl
find_package(OpenGL REQUIRED)

add_definitions(-DDFM2_STATIC_LIBRARY=ON)
set(DELFEM2_INCLUDE_DIR "../../include")

#####################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party
)

file(GLOB 
    DFM2_CPP_FILES
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/*.cpp)
file(GLOB 
    DFM2_HDR_FILES
    ${PROJECT_SOURCE_DIR}/../../include/delfem2/*.h)    
set(DELFEM2_INC         "../../include/delfem2")    
add_library(${PROJ_NAME} STATIC
    ${DFM2_CPP_FILES}
    ${DFM2_HDR_FILES}    
    ${DELFEM2_INC}/opengl/funcs.h           ${DELFEM2_INC}/opengl/funcs.cpp
    ${DELFEM2_INC}/opengl/tex.h             ${DELFEM2_INC}/opengl/tex.cpp
    ${DELFEM2_INC}/opengl/r2t.h             ${DELFEM2_INC}/opengl/r2t.cpp
    
    ${DELFEM2_INC}/opengl/old/v2.h           ${DELFEM2_INC}/opengl/old/v2.cpp
    ${DELFEM2_INC}/opengl/old/v3q.h          ${DELFEM2_INC}/opengl/old/v3q.cpp
    ${DELFEM2_INC}/opengl/old/color.h        ${DELFEM2_INC}/opengl/old/color.cpp
    ${DELFEM2_INC}/opengl/old/funcs.h        ${DELFEM2_INC}/opengl/old/funcs.cpp
    ${DELFEM2_INC}/opengl/old/mshuni.h       ${DELFEM2_INC}/opengl/old/mshuni.cpp
    ${DELFEM2_INC}/opengl/old/cad2dtriv2.h   ${DELFEM2_INC}/opengl/old/cad2dtriv2.cpp
    ${DELFEM2_INC}/opengl/old/caddtri_v3.h   ${DELFEM2_INC}/opengl/old/caddtri_v3.cpp
    ${DELFEM2_INC}/opengl/old/rigv3.h        ${DELFEM2_INC}/opengl/old/rigv3.cpp
    ${DELFEM2_INC}/opengl/old/colorv3.h      ${DELFEM2_INC}/opengl/old/colorv3.cpp
    ${DELFEM2_INC}/opengl/old/r2tglo.h       ${DELFEM2_INC}/opengl/old/r2tglo.cpp
    ${DELFEM2_INC}/opengl/old/gridcube.h     ${DELFEM2_INC}/opengl/old/gridcube.cpp
    ${DELFEM2_INC}/opengl/old/gizmo.h        ${DELFEM2_INC}/opengl/old/gizmo.cpp
    
    ${DELFEM2_INC}/glfw/util.h     
    ${DELFEM2_INC}/glfw/viewer2.h     
    ${DELFEM2_INC}/glfw/viewer2.cpp
    ${DELFEM2_INC}/glfw/viewer3.h     
    ${DELFEM2_INC}/glfw/viewer3.cpp

    ${DELFEM2_INC}/openglstb/glyph.h         ${DELFEM2_INC}/openglstb/glyph.cpp
)

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
        ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/libdfm2_static_glfwold.*
        ${PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
        )
endif()
