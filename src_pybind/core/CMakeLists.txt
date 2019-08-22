cmake_minimum_required(VERSION 3.0)

set(PROJ_NAME c_core)

project(${PROJ_NAME})
set(PYBIND11_CPP_STANDARD -std=c++11)

if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4701 /wd4267 /wd4305 /wd4389 /wd4456 /wd4702 /wd4996 /wd4505 -O2")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -Wno-unused-function -O2")
endif()		

set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")
set(DELFEM2_EXT         "../../src_cpp/external")

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
  ${DELFEM2_INC}/bvh.h              ${DELFEM2_SRC}/bvh.cpp
  ${DELFEM2_INC}/mathexpeval.h      ${DELFEM2_SRC}/mathexpeval.cpp
  ${DELFEM2_INC}/dyntri.h           ${DELFEM2_SRC}/dyntri.cpp
  ${DELFEM2_INC}/mats.h             ${DELFEM2_SRC}/mats.cpp
  ${DELFEM2_INC}/emat.h             ${DELFEM2_SRC}/emat.cpp
  ${DELFEM2_INC}/primitive.h        ${DELFEM2_SRC}/primitive.cpp
  ${DELFEM2_INC}/iss.h              ${DELFEM2_SRC}/iss.cpp

  ${DELFEM2_INC}/ilu_mats.h         ${DELFEM2_SRC}/ilu_mats.cpp
  ${DELFEM2_INC}/fem_emats.h        ${DELFEM2_SRC}/fem_emats.cpp
  ${DELFEM2_INC}/objfunc_v23.h      ${DELFEM2_SRC}/objfunc_v23.cpp
  ${DELFEM2_INC}/objfunc_v23dtri.h  ${DELFEM2_SRC}/objfunc_v23dtri.cpp
  ${DELFEM2_INC}/srchuni_v3.h       ${DELFEM2_SRC}/srchuni_v3.cpp
  ${DELFEM2_INC}/srch_v3bvhmshtopo.h

  ${DELFEM2_INC}/v23m3q.h           ${DELFEM2_SRC}/v23m3q.cpp
  ${DELFEM2_INC}/dyntri_v3.h        ${DELFEM2_SRC}/dyntri_v3.cpp
  ${DELFEM2_INC}/dyntri_v2.h        ${DELFEM2_SRC}/dyntri_v2.cpp
  ${DELFEM2_INC}/mshtopoio.h        ${DELFEM2_SRC}/mshtopoio.cpp
  ${DELFEM2_INC}/cad2d.h            ${DELFEM2_SRC}/cad2d.cpp

  ${DELFEM2_EXT}/io_gltf.h
  ${DELFEM2_EXT}/io_gltf.cpp

  ../py_funcs.h

  py_mshtopoio.cpp
  py_polyline.cpp 
  py_field.cpp
  py_fem.cpp
  py_sdf.cpp
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

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(WIN32)
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND copy 
    "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\c_core.*.pyd"   
    ${NATIVE_PROJECT_SOURCE_DIR}\\..\\..\\PyDelFEM2\\
  COMMENT "Post Event"
)
else()
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
    ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_core.*
    ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/
  COMMENT "Post Event"
)
endif()
