cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

################################

set(CMAKE_CXX_STANDARD 11)
if(MSVC)
  # GTest is static compiled by default
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")  
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
endif()

################################

if(USE_HEADERONLY)
  project(test_hdronly)
ELSE()
  project(test_static)
ENDIF()
set(MY_BINARY_NAME runUnitTests)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../test_inputs")

set(DELFEM2_INCLUDE_DIR "../include")
set(DELFEM2_INC         "../include/delfem2")
set(3RD_PARTY           "../3rd_party")

get_filename_component(GTEST_DIR ${3RD_PARTY}/googletest ABSOLUTE)

enable_testing()

include_directories(
  ${GTEST_DIR}/googletest/include
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}
)

# Add test cpp file

IF(NOT USE_HEADERONLY)
  message("compile as STATIC libary: ${PROJECT_NAME} USE_HEADERONLY=${USE_HEADERONLY}")
  set(SRC_DFM2
      ${DELFEM2_INC}/vec2.h                 ${DELFEM2_INC}/vec2.cpp
      ${DELFEM2_INC}/vec3.h                 ${DELFEM2_INC}/vec3.cpp
      ${DELFEM2_INC}/mat3.h                 ${DELFEM2_INC}/mat3.cpp
      ${DELFEM2_INC}/mat4.h                 ${DELFEM2_INC}/mat4.cpp
      ${DELFEM2_INC}/quat.h                 ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/camera.h               ${DELFEM2_INC}/camera.cpp
      ${DELFEM2_INC}/primitive.h            ${DELFEM2_INC}/primitive.cpp
      ${DELFEM2_INC}/specialfuncs.h         ${DELFEM2_INC}/specialfuncs.cpp      
      ${DELFEM2_INC}/mshtopo.h              ${DELFEM2_INC}/mshtopo.cpp
      ${DELFEM2_INC}/mshmisc.h              ${DELFEM2_INC}/mshmisc.cpp
      ${DELFEM2_INC}/mshio.h                ${DELFEM2_INC}/mshio.cpp
      ${DELFEM2_INC}/dtri.h                 ${DELFEM2_INC}/dtri.cpp
      ${DELFEM2_INC}/funcs.h                ${DELFEM2_INC}/funcs.cpp
      ${DELFEM2_INC}/mats.h                 ${DELFEM2_INC}/mats.cpp
      ${DELFEM2_INC}/vecxitrsol.h           ${DELFEM2_INC}/vecxitrsol.cpp
      ${DELFEM2_INC}/emat.h                 ${DELFEM2_INC}/emat.cpp
      ${DELFEM2_INC}/bvh.h                  ${DELFEM2_INC}/bvh.cpp
      ${DELFEM2_INC}/slice.h                ${DELFEM2_INC}/slice.cpp
      ${DELFEM2_INC}/bv.h
      
      ${DELFEM2_INC}/geo3_v23m34q.h         ${DELFEM2_INC}/geo3_v23m34q.cpp
      ${DELFEM2_INC}/fem_emats.h            ${DELFEM2_INC}/fem_emats.cpp
      ${DELFEM2_INC}/dtri2_v2dtri.h         ${DELFEM2_INC}/dtri2_v2dtri.cpp
      ${DELFEM2_INC}/ilu_mats.h             ${DELFEM2_INC}/ilu_mats.cpp
      ${DELFEM2_INC}/rig_geo3.h             ${DELFEM2_INC}/rig_geo3.cpp
      ${DELFEM2_INC}/objf_geo3.h            ${DELFEM2_INC}/objf_geo3.cpp
      ${DELFEM2_INC}/femrod.h               ${DELFEM2_INC}/femrod.cpp
      ${DELFEM2_INC}/srchuni_v3.h           ${DELFEM2_INC}/srchuni_v3.cpp
      ${DELFEM2_INC}/gridvoxel.h            ${DELFEM2_INC}/gridvoxel.cpp
      ${DELFEM2_INC}/gridcube.h             ${DELFEM2_INC}/gridcube.cpp            
      ${DELFEM2_INC}/srchbi_v3bvh.h
      ${DELFEM2_INC}/srch_v3bvhmshtopo.h
      )
ELSE()
  add_definitions(-DDFM2_HEADER_ONLY=ON)
  message("compile as HEADER-ONLY library: ${PROJECT_NAME} ${USE_HEADERONLY}")
ENDIF()

add_executable(${MY_BINARY_NAME}
  ${DELFEM2_INC}/evalmathexp.h          ${DELFEM2_INC}/evalmathexp.cpp
  ${DELFEM2_INC}/lp.h                   ${DELFEM2_INC}/lp.cpp
  
 ${3RD_PARTY}/io_gltf.h                ${3RD_PARTY}/io_gltf.cpp
  
  ${SRC_DFM2}
  test_bvh.cpp
  test_lp.cpp
  test_gltf.cpp
  test_fem.cpp
  main.cpp  
)

if(MSVC)
  if(CMAKE_CL_64)
    set(GTEST_MAIN_LIBRARY ${GTEST_DIR}/buildVS64/lib/Release/gtest_main.lib)
    set(GTEST_LIBRARY      ${GTEST_DIR}/buildVS64/lib/Release/gtest.lib)
  else()
    set(GTEST_MAIN_LIBRARY ${GTEST_DIR}/buildVS32/lib/Release/gtest_main.lib)
    set(GTEST_MAIN_LIBRARY ${GTEST_DIR}/buildVS32/lib/Release/gtest.lib)
  endif()
  target_link_libraries(${MY_BINARY_NAME}
    ${GTEST_LIBRARY}
    ${GTEST_MAIN_LIBRARY})  
else()
  target_link_libraries(${MY_BINARY_NAME}
      ${GTEST_DIR}/lib/libgtest.a
      ${GTEST_DIR}/lib/libgtest_main.a
      -pthread)
endif()    

add_test(
  NAME ${MY_BINARY_NAME}
  COMMAND ${MY_BINARY_NAME}
)
