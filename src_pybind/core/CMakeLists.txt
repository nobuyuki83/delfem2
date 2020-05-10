cmake_minimum_required(VERSION 3.0)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

#########################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4702 /wd4996 \
    /wd4701 /wd4267 /wd4305 /wd4389 /wd4456 /wd4505 ")
else()
  set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra \
    -Wno-deprecated-declarations \
    -Wno-unused-parameter \
    -Wno-unused-function")
endif()

##########################################################

set(PROJ_NAME c_core)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")

set(PYBIND11_DIR ${3RD_PARTY}/pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

########################

include_directories(
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}/glad3/include
  ${3RD_PARTY}/pybind11/include
  ${3RD_PARTY}
)

if(NOT USE_HEADERONLY)
  set(SRC_DFM2
      ${DELFEM2_INC}/vec2.h               ${DELFEM2_INC}/vec2.cpp
      ${DELFEM2_INC}/vec3.h               ${DELFEM2_INC}/vec3.cpp
      ${DELFEM2_INC}/mat3.h               ${DELFEM2_INC}/mat3.cpp
      ${DELFEM2_INC}/mat4.h               ${DELFEM2_INC}/mat4.cpp
      ${DELFEM2_INC}/quat.h               ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/primitive.h          ${DELFEM2_INC}/primitive.cpp
      ${DELFEM2_INC}/mshtopo.h            ${DELFEM2_INC}/mshtopo.cpp
      ${DELFEM2_INC}/mshmisc.h            ${DELFEM2_INC}/mshmisc.cpp
      ${DELFEM2_INC}/mshio.h              ${DELFEM2_INC}/mshio.cpp
      ${DELFEM2_INC}/funcs.h              ${DELFEM2_INC}/funcs.cpp
      ${DELFEM2_INC}/color.h              ${DELFEM2_INC}/color.cpp
      ${DELFEM2_INC}/dtri.h               ${DELFEM2_INC}/dtri.cpp
      ${DELFEM2_INC}/mats.h               ${DELFEM2_INC}/mats.cpp
      ${DELFEM2_INC}/bvh.h                ${DELFEM2_INC}/bvh.cpp
      ${DELFEM2_INC}/emat.h               ${DELFEM2_INC}/emat.cpp
      ${DELFEM2_INC}/iss.h                ${DELFEM2_INC}/iss.cpp      
      ${DELFEM2_INC}/vecxitrsol.h         ${DELFEM2_INC}/vecxitrsol.cpp
      ${DELFEM2_INC}/gridvoxel.h          ${DELFEM2_INC}/gridvoxel.cpp
      ${DELFEM2_INC}/gridcube.h           ${DELFEM2_INC}/gridcube.cpp            
      ${DELFEM2_INC}/bv.h
      
      ${DELFEM2_INC}/geo3_v23m34q.h       ${DELFEM2_INC}/geo3_v23m34q.cpp
      ${DELFEM2_INC}/rig_geo3.h           ${DELFEM2_INC}/rig_geo3.cpp
      ${DELFEM2_INC}/dtri2_v2dtri.h       ${DELFEM2_INC}/dtri2_v2dtri.cpp
      ${DELFEM2_INC}/dtri3_v3dtri.h       ${DELFEM2_INC}/dtri3_v3dtri.cpp
      ${DELFEM2_INC}/cad2_dtri2.h         ${DELFEM2_INC}/cad2_dtri2.cpp
      ${DELFEM2_INC}/fem_emats.h          ${DELFEM2_INC}/fem_emats.cpp
      ${DELFEM2_INC}/srchuni_v3.h         ${DELFEM2_INC}/srchuni_v3.cpp
      ${DELFEM2_INC}/ilu_mats.h           ${DELFEM2_INC}/ilu_mats.cpp      
      ${DELFEM2_INC}/objf_geo3.h             ${DELFEM2_INC}/objf_geo3.cpp
      ${DELFEM2_INC}/objfdtri_objfdtri23.h   ${DELFEM2_INC}/objfdtri_objfdtri23.cpp
      ${DELFEM2_INC}/srch_v3bvhmshtopo.h
      )
else()
  add_definitions(-DDFM2_HEADER_ONLY=ON)
endif()

set(INPUT_SRC
    ${DELFEM2_INC}/evalmathexp.h        ${DELFEM2_INC}/evalmathexp.cpp
    ${DELFEM2_INC}/slice.h              ${DELFEM2_INC}/slice.cpp
    ${DELFEM2_INC}/mshtopoio.h          ${DELFEM2_INC}/mshtopoio.cpp
    ${3RD_PARTY}/io_gltf.h              ${3RD_PARTY}/io_gltf.cpp
    
    ${SRC_DFM2}
    
    ../py_funcs.h

    py_mshtopoio.cpp
    py_polyline.cpp
    py_field.cpp
    py_fem.cpp
    py_sdf.cpp

    main.cpp
)

###########################################

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
