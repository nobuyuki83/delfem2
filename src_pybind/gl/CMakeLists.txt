cmake_minimum_required(VERSION 3.0)

set(PROJ_NAME c_gl)

project(${PROJ_NAME})
set(PYBIND11_CPP_STANDARD -std=c++14)
if(MSVC)
set(CMAKE_CXX_FLAGS "/W4  -O2 \
    /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4701 \
    /wd4267 /wd4305 /wd4389 /wd4456 /wd4702 /wd4996 \
    /wd4505 /wd4245 /wd4457 /wd4800 /wd4838")
else()
set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra \
    -Wno-deprecated-declarations \
    -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable")
endif()

set(GLAD_INCLUDE_DIR "../glad3/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

set(INPUT_SRC
  ${DELFEM2_INC}/vec2.h          ${DELFEM2_SRC}/vec2.cpp
  ${DELFEM2_INC}/vec3.h          ${DELFEM2_SRC}/vec3.cpp
  ${DELFEM2_INC}/mat3.h          ${DELFEM2_SRC}/mat3.cpp
  ${DELFEM2_INC}/quat.h          ${DELFEM2_SRC}/quat.cpp
  ${DELFEM2_INC}/msh.h           ${DELFEM2_SRC}/msh.cpp
  ${DELFEM2_INC}/voxel.h         ${DELFEM2_SRC}/voxel.cpp
  ${DELFEM2_INC}/camera.h        ${DELFEM2_SRC}/camera.cpp  
  ${DELFEM2_INC}/bv.h
  ${DELFEM2_INC}/color.h        ${DELFEM2_SRC}/color.cpp    

  ${DELFEM2_INC}/gl_gpusampler.h ${DELFEM2_SRC}/gl_gpusampler.cpp
  ${DELFEM2_INC}/gl_tex.h        ${DELFEM2_SRC}/gl_tex.cpp
  ${DELFEM2_INC}/gl24_funcs.h    ${DELFEM2_SRC}/gl24_funcs.cpp
  ${DELFEM2_INC}/gl2_funcs.h     ${DELFEM2_SRC}/gl2_funcs.cpp
  ${DELFEM2_INC}/gl2ew_funcs.h   ${DELFEM2_SRC}/gl2ew_funcs.cpp
  ${DELFEM2_INC}/gl2_color.h      ${DELFEM2_SRC}/gl2_color.cpp  
  ${DELFEM2_INC}/gl_v23.h        ${DELFEM2_SRC}/gl_v23.cpp
  ${DELFEM2_INC}/gl_voxbv.h      ${DELFEM2_SRC}/gl_voxbv.cpp

  ${DELFEM2_INC}/gl_cad_dyntri_v23.h
  ${DELFEM2_SRC}/gl_cad_dyntri_v23.cpp
  
  ../glad3/include/glad/glad.h  ../glad3/src/glad.c

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

include_directories( ${INPUT_INCLUDE_DIR} )
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
