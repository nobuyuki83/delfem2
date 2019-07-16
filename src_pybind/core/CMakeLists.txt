cmake_minimum_required(VERSION 3.0)

set(PROJ_NAME c_core)

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
  ${DELFEM2_INC}/vec2.h             ${DELFEM2_SRC}/vec2.cpp
  ${DELFEM2_INC}/vec3.h             ${DELFEM2_SRC}/vec3.cpp
  ${DELFEM2_INC}/mat3.h             ${DELFEM2_SRC}/mat3.cpp
  ${DELFEM2_INC}/quat.h             ${DELFEM2_SRC}/quat.cpp
  ${DELFEM2_INC}/msh.h              ${DELFEM2_SRC}/msh.cpp
  ${DELFEM2_INC}/mshtopo.h          ${DELFEM2_SRC}/mshtopo.cpp
  ${DELFEM2_INC}/mshio.h            ${DELFEM2_SRC}/mshio.cpp
  ${DELFEM2_INC}/funcs.h            ${DELFEM2_SRC}/funcs.cpp
  ${DELFEM2_INC}/voxel.h            ${DELFEM2_SRC}/voxel.cpp
  ${DELFEM2_INC}/bv.h
  ${DELFEM2_INC}/mathexpeval.h      ${DELFEM2_SRC}/mathexpeval.cpp
  ${DELFEM2_INC}/dyntri.h           ${DELFEM2_SRC}/dyntri.cpp

  ${DELFEM2_INC}/matrix_sparse.h    ${DELFEM2_SRC}/matrix_sparse.cpp
  ${DELFEM2_INC}/ilu_sparse.h       ${DELFEM2_SRC}/ilu_sparse.cpp
  ${DELFEM2_INC}/fem.h              ${DELFEM2_SRC}/fem.cpp
  ${DELFEM2_INC}/fem_ematrix.h      ${DELFEM2_SRC}/fem_ematrix.cpp
  ${DELFEM2_INC}/objfunc_v23.h      ${DELFEM2_SRC}/objfunc_v23.cpp

  ${DELFEM2_INC}/v23m3q.h           ${DELFEM2_SRC}/v23m3q.cpp
  ${DELFEM2_INC}/dyntri_v3.h        ${DELFEM2_SRC}/dyntri_v3.cpp
  ${DELFEM2_INC}/dyntri_v2.h        ${DELFEM2_SRC}/dyntri_v2.cpp
  ${DELFEM2_INC}/mshtopoio.h        ${DELFEM2_SRC}/mshtopoio.cpp
  ${DELFEM2_INC}/cad2d.h            ${DELFEM2_SRC}/cad2d.cpp
  ${DELFEM2_INC}/sdf.h              ${DELFEM2_SRC}/sdf.cpp
  ${DELFEM2_INC}/hash3d.h           ${DELFEM2_SRC}/hash3d.cpp

  ${DELFEM2_INC}/isosurface_stuffing.h
  ${DELFEM2_SRC}/isosurface_stuffing.cpp

  py_mshtopoio.cpp
  py_polyline.cpp 
  py_field.cpp
  py_fem.cpp
  main.cpp
)

set(INPUT_INCLUDE_DIR
  ${DELFEM2_INCLUDE_DIR}
)
set(INPUT_LIBRARY
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
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_core.*
  ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/
  COMMENT "Post Event"
)
endif()
