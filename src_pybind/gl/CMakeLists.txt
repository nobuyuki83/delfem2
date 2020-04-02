cmake_minimum_required(VERSION 3.0)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

#########################################################

if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4  -O2 \
    /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4701 \
    /wd4267 /wd4305 /wd4389 /wd4456 /wd4702 /wd4996 \
    /wd4505 /wd4245 /wd4457 /wd4800 /wd4838")
else()
  set(PYBIND11_CPP_STANDARD -std=c++11)
  set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra \
    -Wno-deprecated-declarations \
    -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable")
endif()

################################

set(PROJ_NAME c_gl)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")

set(GLAD_INCLUDE_DIR "${3RD_PARTY}/glad3/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

set(PYBIND11_DIR ${3RD_PARTY}/pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

if(NOT USE_HEADERONLY)
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  set(SRC_DFM2
      ${DELFEM2_INC}/quat.h                      ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/vec2.h                      ${DELFEM2_INC}/vec2.cpp
      ${DELFEM2_INC}/vec3.h                      ${DELFEM2_INC}/vec3.cpp
      ${DELFEM2_INC}/mat3.h                      ${DELFEM2_INC}/mat3.cpp
      ${DELFEM2_INC}/mat4.h                      ${DELFEM2_INC}/mat4.cpp
      ${DELFEM2_INC}/camera.h                    ${DELFEM2_INC}/camera.cpp
      ${DELFEM2_INC}/color.h                     ${DELFEM2_INC}/color.cpp
      ${DELFEM2_INC}/v23m34q.h                   ${DELFEM2_INC}/v23m34q.cpp
      )
endif()

set(INPUT_SRC
    ${DELFEM2_INC}/mshmisc.h                   ${DELFEM2_INC}/mshmisc.cpp
    ${DELFEM2_INC}/voxel.h                     ${DELFEM2_INC}/voxel.cpp
    ${DELFEM2_INC}/bv.h

    ${DELFEM2_INC}/opengl/gl_funcs.h           ${DELFEM2_INC}/opengl/gl_funcs.cpp
    ${DELFEM2_INC}/opengl/gl_tex.h             ${DELFEM2_INC}/opengl/gl_tex.cpp
    ${DELFEM2_INC}/opengl/render2tex_gl.h      ${DELFEM2_INC}/opengl/render2tex_gl.cpp

    ${DELFEM2_INC}/opengl/render2tex_glold.h   ${DELFEM2_INC}/opengl/render2tex_glold.cpp        
    ${DELFEM2_INC}/opengl/glold_voxbv.h        ${DELFEM2_INC}/opengl/glold_voxbv.cpp
    ${DELFEM2_INC}/opengl/glold_funcs.h        ${DELFEM2_INC}/opengl/glold_funcs.cpp
    ${DELFEM2_INC}/opengl/glold_color.h        ${DELFEM2_INC}/opengl/glold_color.cpp
    ${DELFEM2_INC}/opengl/v3q_glold.h          ${DELFEM2_INC}/opengl/v3q_glold.cpp
    ${DELFEM2_INC}/opengl/v2_glold.h           ${DELFEM2_INC}/opengl/v2_glold.cpp    
    ${DELFEM2_INC}/opengl/caddtri_v2_glold.h   ${DELFEM2_INC}/opengl/caddtri_v2_glold.cpp
    ${DELFEM2_INC}/opengl/caddtri_v3_glold.h   ${DELFEM2_INC}/opengl/caddtri_v3_glold.cpp
  
    ${3RD_PARTY}/glad3/include/glad/glad.h
    ${3RD_PARTY}/glad3/src/glad.c

    ${SRC_DFM2}

    py_sampler.cpp
    py_texture.cpp
    main.cpp
)


find_package(OpenGL REQUIRED)

set(INPUT_INCLUDE_DIR
  ${OPENGL_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}  
  ${DELFEM2_INCLUDE_DIR}
)
set(INPUT_LIBRARY
  ${OPENGL_LIBRARY}
)

include_directories(
  ${INPUT_INCLUDE_DIR}
  ${3RD_PARTY}/glad3/include
  ${3RD_PARTY}/pybind11/include)
pybind11_add_module(${PROJ_NAME}
  ${INPUT_SRC} )

target_link_libraries(${PROJ_NAME}
  PRIVATE ${INPUT_LIBRARY})

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(MSVC)
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND copy "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\c_gl.*.pyd"   ${NATIVE_PROJECT_SOURCE_DIR}\\..\\..\\PyDelFEM2\\gl\\
  COMMENT "Post Event"
)
else()
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_gl.*
  ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/gl
  COMMENT "Post Event"
)
endif()
