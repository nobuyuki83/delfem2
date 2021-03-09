cmake_minimum_required(VERSION 3.12)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

################################

set(CMAKE_CXX_STANDARD 11)
if(MSVC)
  # GTest is static compiled by default
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")  
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456")
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

find_package( GTest REQUIRED )
find_package(Threads REQUIRED)

enable_testing()

include_directories(
  ${GTEST_INCLUDE_DIR}
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
      ${DELFEM2_INC}/geo3_v23m34q.h         ${DELFEM2_INC}/geo3_v23m34q.cpp
      ${DELFEM2_INC}/cam3_m4q.h             ${DELFEM2_INC}/cam3_m4q.cpp
      ${DELFEM2_INC}/geoplygn2_v2.h         ${DELFEM2_INC}/geoplygn2_v2.cpp
      ${DELFEM2_INC}/geoconvhull3_v3.h      ${DELFEM2_INC}/geoconvhull3_v3.cpp
      ${DELFEM2_INC}/geoproximity3_v3.h     ${DELFEM2_INC}/geoproximity3_v3.cpp
      ${DELFEM2_INC}/geosolidelm_v3.h       ${DELFEM2_INC}/geosolidelm_v3.cpp

      ${DELFEM2_INC}/lp.h                   ${DELFEM2_INC}/lp.cpp
      ${DELFEM2_INC}/evalmathexp.h          ${DELFEM2_INC}/evalmathexp.cpp
      ${DELFEM2_INC}/specialfuncs.h         ${DELFEM2_INC}/specialfuncs.cpp      
      ${DELFEM2_INC}/jagarray.h             ${DELFEM2_INC}/jagarray.cpp
      ${DELFEM2_INC}/file.h                 ${DELFEM2_INC}/file.cpp
      ${DELFEM2_INC}/str.h                  ${DELFEM2_INC}/str.cpp
      ${DELFEM2_INC}/filenpy_str.h          ${DELFEM2_INC}/filenpy_str.cpp

      ${DELFEM2_INC}/dtri.h                 ${DELFEM2_INC}/dtri.cpp
      ${DELFEM2_INC}/dtri2_v2dtri.h         ${DELFEM2_INC}/dtri2_v2dtri.cpp
      ${DELFEM2_INC}/cad2_dtri2.h           ${DELFEM2_INC}/cad2_dtri2.cpp

      ${DELFEM2_INC}/mshprimitive.h         ${DELFEM2_INC}/mshprimitive.cpp
      ${DELFEM2_INC}/mshuni.h               ${DELFEM2_INC}/mshuni.cpp
      ${DELFEM2_INC}/mshmisc.h              ${DELFEM2_INC}/mshmisc.cpp
      ${DELFEM2_INC}/points.h               ${DELFEM2_INC}/points.cpp
      ${DELFEM2_INC}/mshio.h                ${DELFEM2_INC}/mshio.cpp
      ${DELFEM2_INC}/slice.h                ${DELFEM2_INC}/slice.cpp
      
      ${DELFEM2_INC}/lsmats.h               ${DELFEM2_INC}/lsmats.cpp
      ${DELFEM2_INC}/lsilu_mats.h           ${DELFEM2_INC}/lsilu_mats.cpp
      ${DELFEM2_INC}/vecxitrsol.h           ${DELFEM2_INC}/vecxitrsol.cpp
      
      ${DELFEM2_INC}/femutil.h              ${DELFEM2_INC}/femutil.cpp      
      ${DELFEM2_INC}/fempoisson.h           ${DELFEM2_INC}/fempoisson.cpp
      ${DELFEM2_INC}/femdkt.h               ${DELFEM2_INC}/femdkt.cpp
      ${DELFEM2_INC}/femrod.h               ${DELFEM2_INC}/femrod.cpp
      ${DELFEM2_INC}/femmitc3.h             ${DELFEM2_INC}/femmitc3.cpp
      ${DELFEM2_INC}/femhelmholtz.h         ${DELFEM2_INC}/femhelmholtz.cpp
      ${DELFEM2_INC}/femstokes.h            ${DELFEM2_INC}/femstokes.cpp
      ${DELFEM2_INC}/femnavierstokes.h      ${DELFEM2_INC}/femnavierstokes.cpp
      ${DELFEM2_INC}/femsolidlinear.h       ${DELFEM2_INC}/femsolidlinear.cpp
      ${DELFEM2_INC}/femmips_geo3.h         ${DELFEM2_INC}/femmips_geo3.cpp

      ${DELFEM2_INC}/pbd_geo3.h             ${DELFEM2_INC}/pbd_geo3.cpp

      ${DELFEM2_INC}/defarapenergy_geo3.h   ${DELFEM2_INC}/defarapenergy_geo3.cpp

      ${DELFEM2_INC}/srchbvh.h              ${DELFEM2_INC}/srchbvh.cpp
      ${DELFEM2_INC}/srchbv3aabb.h
      ${DELFEM2_INC}/srchbv3sphere.h
      ${DELFEM2_INC}/srchbi_v3bvh.h
      ${DELFEM2_INC}/srch_v3bvhmshtopo.h
      ${DELFEM2_INC}/srchuni_v3.h           ${DELFEM2_INC}/srchuni_v3.cpp

      ${DELFEM2_INC}/rig_geo3.h             ${DELFEM2_INC}/rig_geo3.cpp
      ${DELFEM2_INC}/rigopt.h

      ${DELFEM2_INC}/gridvoxel.h            ${DELFEM2_INC}/gridvoxel.cpp
      ${DELFEM2_INC}/gridcube.h             ${DELFEM2_INC}/gridcube.cpp
      )
ELSE()
  add_definitions(-DDFM2_HEADER_ONLY=ON)
  message("compile as HEADER-ONLY library: ${PROJECT_NAME} ${USE_HEADERONLY}")
ENDIF()

add_executable(${MY_BINARY_NAME}
  ${DELFEM2_INC}/tinygltf/io_gltf.h     ${DELFEM2_INC}/tinygltf/io_gltf.cpp
  ${SRC_DFM2}
  test_bvh.cpp
  test_lp.cpp
  test_gltf.cpp
  test_fem.cpp
  test_geo.cpp
  test_cad.cpp
  test_msh.cpp
  main.cpp  
)

target_link_libraries(${MY_BINARY_NAME}
  ${GTEST_LIBRARY}
  ${GTEST_MAIN_LIBRARY}
  Threads::Threads
)

add_test(
  NAME ${MY_BINARY_NAME}
  COMMAND ${MY_BINARY_NAME}
)
