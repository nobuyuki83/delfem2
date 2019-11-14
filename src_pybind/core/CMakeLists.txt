cmake_minimum_required(VERSION 3.0)

set(PROJ_NAME c_core)
project(${PROJ_NAME})

if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4702 /wd4996 \
    /wd4701 /wd4267 /wd4305 /wd4389 /wd4456 /wd4505 ")
else()
  set(PYBIND11_CPP_STANDARD -std=c++11)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -Wno-unused-function -O2")
endif()		

set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(DELFEM2_EXT         "../../src_cpp/external")

set(INPUT_SRC
  ${DELFEM2_INC}/vec2.h               ${DELFEM2_INC}/vec2.cpp
  ${DELFEM2_INC}/vec3.h               ${DELFEM2_INC}/vec3.cpp
  ${DELFEM2_INC}/mat3.h               ${DELFEM2_INC}/mat3.cpp
  ${DELFEM2_INC}/quat.h               ${DELFEM2_INC}/quat.cpp
  ${DELFEM2_INC}/mshmisc.h            ${DELFEM2_INC}/mshmisc.cpp
  ${DELFEM2_INC}/mshtopo.h            ${DELFEM2_INC}/mshtopo.cpp
  ${DELFEM2_INC}/mshio.h              ${DELFEM2_INC}/mshio.cpp
  ${DELFEM2_INC}/funcs.h              ${DELFEM2_INC}/funcs.cpp
  ${DELFEM2_INC}/voxel.h              ${DELFEM2_INC}/voxel.cpp
  ${DELFEM2_INC}/evalmathexp.h        ${DELFEM2_INC}/evalmathexp.cpp
  ${DELFEM2_INC}/dtri.h               ${DELFEM2_INC}/dtri.cpp
  ${DELFEM2_INC}/mats.h               ${DELFEM2_INC}/mats.cpp
  ${DELFEM2_INC}/emat.h               ${DELFEM2_INC}/emat.cpp
  ${DELFEM2_INC}/primitive.h          ${DELFEM2_INC}/primitive.cpp
  ${DELFEM2_INC}/iss.h                ${DELFEM2_INC}/iss.cpp
  ${DELFEM2_INC}/color.h              ${DELFEM2_INC}/color.cpp
  ${DELFEM2_INC}/slice.h              ${DELFEM2_INC}/slice.cpp
  ${DELFEM2_INC}/bvh.h                ${DELFEM2_INC}/bvh.cpp
  ${DELFEM2_INC}/bv.h
  
  ${DELFEM2_INC}/ilu_mats.h           ${DELFEM2_INC}/ilu_mats.cpp
  ${DELFEM2_INC}/fem_emats.h          ${DELFEM2_INC}/fem_emats.cpp
  ${DELFEM2_INC}/objfunc_v23.h        ${DELFEM2_INC}/objfunc_v23.cpp
  ${DELFEM2_INC}/objfunc_v23dtri.h    ${DELFEM2_INC}/objfunc_v23dtri.cpp
  ${DELFEM2_INC}/rig_v3q.h            ${DELFEM2_INC}/rig_v3q.cpp
  ${DELFEM2_INC}/v23m3q.h             ${DELFEM2_INC}/v23m3q.cpp
  ${DELFEM2_INC}/dtri_v3.h            ${DELFEM2_INC}/dtri_v3.cpp
  ${DELFEM2_INC}/dtri_v2.h            ${DELFEM2_INC}/dtri_v2.cpp
  ${DELFEM2_INC}/mshtopoio.h          ${DELFEM2_INC}/mshtopoio.cpp
  ${DELFEM2_INC}/cad2d.h              ${DELFEM2_INC}/cad2d.cpp
  ${DELFEM2_INC}/eigen_rigidbody.h    ${DELFEM2_INC}/eigen_rigidbody.cpp
  ${DELFEM2_INC}/srchuni_v3.h         ${DELFEM2_INC}/srchuni_v3.cpp
  ${DELFEM2_INC}/srch_v3bvhmshtopo.h

  ${DELFEM2_EXT}/io_gltf.h            ${DELFEM2_EXT}/io_gltf.cpp

  ../py_funcs.h

  py_mshtopoio.cpp
  py_polyline.cpp 
  py_field.cpp
  py_fem.cpp
  py_sdf.cpp
  py_rigidbody.cpp
  
  main.cpp
)

set(INPUT_INCLUDE_DIR
  "../../src_cpp/external/eigen-git-mirror"
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
