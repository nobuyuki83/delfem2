cmake_minimum_required(VERSION 3.0)

set(PROJ_NAME c_gl)

project(${PROJ_NAME})
set(PYBIND11_CPP_STANDARD -std=c++14)
if(WIN32)
set(CMAKE_CXX_FLAGS "-Wall -O2")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -O2")
endif()		

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
  ${DELFEM2_INC}/bv.h

  ${DELFEM2_INC}/gl_v23q.h       ${DELFEM2_SRC}/gl_v23q.cpp
  ${DELFEM2_INC}/gl_funcs.h      ${DELFEM2_SRC}/gl_funcs.cpp
  ${DELFEM2_INC}/gl_color.h      ${DELFEM2_SRC}/gl_color.cpp
  ${DELFEM2_INC}/glew_funcs.h
  ${DELFEM2_INC}/gl_camera.h     ${DELFEM2_SRC}/gl_camera.cpp
  ${DELFEM2_INC}/gl_voxbv.h      ${DELFEM2_SRC}/gl_voxbv.cpp
  ${DELFEM2_INC}/gl_gpusampler.h ${DELFEM2_SRC}/gl_gpusampler.cpp
  ${DELFEM2_INC}/gl_tex.h        ${DELFEM2_SRC}/gl_tex.cpp

  ${DELFEM2_INC}/gl_cad_dyntri_v23.h
  ${DELFEM2_SRC}/gl_cad_dyntri_v23.cpp

  py_sampler.cpp
  py_texture.cpp 
  main.cpp
)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)

set(INPUT_INCLUDE_DIR
  ${OPENGL_INCLUDE_DIR}
  ${GLEW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)
set(INPUT_LIBRARY
  ${OPENGL_LIBRARY}
  ${GLEW_LIBRARY}
)

include_directories( ${INPUT_INCLUDE_DIR} )
pybind11_add_module(${PROJ_NAME}
  ${INPUT_SRC} )

target_link_libraries(${PROJ_NAME}
  PRIVATE ${INPUT_LIBRARY})

if(WIN32)
else()
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_gl.*
  ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/gl
  COMMENT "Post Event"
)
endif()
