cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
project(testsDelFEM2)
set(MY_BINARY_NAME runUnitTests)

if(MSVC)
  # GTest is static compiled by default
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")    
  set(CMAKE_CXX_FLAGS "/wd4267")
endif()

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../test_inputs")

set(DELFEM2_INCLUDE_DIR "../include")
set(DELFEM2_INC         "../include/delfem2")

enable_testing()

if(MSVC)
  if(CMAKE_CL_64)
    set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS64/lib/Release/gtest_main.lib)
    set(GTEST_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS64/lib/Release/gtest.lib)
  else()
    set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS32/lib/Release/gtest_main.lib)
    set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/buildVS32/lib/Release/gtest.lib)
  endif()
else()
  set(GTEST_MAIN_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/lib/libgtest_main.a)
  set(GTEST_LIBRARY ${PROJECT_SOURCE_DIR}/googletest/lib/libgtest.a)
endif()
set(GTEST_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/googletest/googletest/include)

include_directories(
  ${GTEST_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

# Add test cpp file
add_executable(${MY_BINARY_NAME}
  ${DELFEM2_INC}/vec3.h           ${DELFEM2_INC}/vec3.cpp
  ${DELFEM2_INC}/vec2.h           ${DELFEM2_INC}/vec2.cpp
  ${DELFEM2_INC}/mat3.h           ${DELFEM2_INC}/mat3.cpp
  ${DELFEM2_INC}/quat.h           ${DELFEM2_INC}/quat.cpp
  ${DELFEM2_INC}/funcs.h          ${DELFEM2_INC}/funcs.cpp
  ${DELFEM2_INC}/specialfuncs.h   ${DELFEM2_INC}/specialfuncs.cpp
  ${DELFEM2_INC}/mshmisc.h        ${DELFEM2_INC}/mshmisc.cpp
  ${DELFEM2_INC}/mshtopo.h        ${DELFEM2_INC}/mshtopo.cpp
  ${DELFEM2_INC}/mshio.h          ${DELFEM2_INC}/mshio.cpp
  ${DELFEM2_INC}/voxel.h          ${DELFEM2_INC}/voxel.cpp
  ${DELFEM2_INC}/evalmathexp.h    ${DELFEM2_INC}/evalmathexp.cpp
  ${DELFEM2_INC}/lp.h             ${DELFEM2_INC}/lp.cpp
  ${DELFEM2_INC}/bvh.h            ${DELFEM2_INC}/bvh.cpp
  ${DELFEM2_INC}/primitive.h      ${DELFEM2_INC}/primitive.cpp
  ${DELFEM2_INC}/slice.h          ${DELFEM2_INC}/slice.cpp
  ${DELFEM2_INC}/emat.h           ${DELFEM2_INC}/emat.cpp
  ${DELFEM2_INC}/dtri.h           ${DELFEM2_INC}/dtri.cpp
  ${DELFEM2_INC}/mats.h           ${DELFEM2_INC}/mats.cpp
  ${DELFEM2_INC}/bv.h  

  ${DELFEM2_INC}/fem_emats.h      ${DELFEM2_INC}/fem_emats.cpp
  ${DELFEM2_INC}/ilu_mats.h       ${DELFEM2_INC}/ilu_mats.cpp
  ${DELFEM2_INC}/dtri_v2.h        ${DELFEM2_INC}/dtri_v2.cpp
  ${DELFEM2_INC}/objfunc_v23.h    ${DELFEM2_INC}/objfunc_v23.cpp
  ${DELFEM2_INC}/srchuni_v3.h     ${DELFEM2_INC}/srchuni_v3.cpp
  ${DELFEM2_INC}/srchbi_v3bvh.h
  ${DELFEM2_INC}/srch_v3bvhmshtopo.h

  test_bvh.cpp
  test_lp.cpp
  test_gltf.cpp
  test_fem.cpp
  main.cpp  
)


if(MSVC)
  target_link_libraries(${MY_BINARY_NAME}
    ${GTEST_LIBRARY}
    ${GTEST_MAIN_LIBRARY})  
else()    
  target_link_libraries(${MY_BINARY_NAME}
    ${GTEST_LIBRARY}
    ${GTEST_MAIN_LIBRARY}
    -pthread)    
endif()    


add_test(
  NAME ${MY_BINARY_NAME}
  COMMAND ${MY_BINARY_NAME}
)
