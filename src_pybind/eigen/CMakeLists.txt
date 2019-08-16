cmake_minimum_required(VERSION 3.0)

set(PROJ_NAME c_eigen)

project(${PROJ_NAME})
set(PYBIND11_CPP_STANDARD -std=c++14)
if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4701 /wd4267 /wd4305 /wd4389 /wd4456 /wd4702 /wd4996 /wd4505 /wd4127 -O2")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -O2")
endif()		

set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

set(INPUT_SRC
  ${DELFEM2_SRC}/vec2.cpp          ${DELFEM2_INC}/vec2.h
  ${DELFEM2_SRC}/vec3.cpp          ${DELFEM2_INC}/vec3.h
  ${DELFEM2_SRC}/mat3.cpp          ${DELFEM2_INC}/mat3.h
  ${DELFEM2_SRC}/quat.cpp          ${DELFEM2_INC}/quat.h

  ${DELFEM2_SRC}/v23m3q.cpp        ${DELFEM2_INC}/v23m3q.h
  ${DELFEM2_SRC}/eigen_rigidbody.cpp     ${DELFEM2_INC}/eigen_rigidbody.h

  ${DELFEM2_SRC}/gl_v23.cpp       ${DELFEM2_INC}/gl_v23.h

  main.cpp
)

find_package(OpenGL REQUIRED)

set(INPUT_INCLUDE_DIR
  ${OPENGL_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  "../../src_cpp/external/eigen-git-mirror"
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
  COMMAND copy "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\c_eigen.*.pyd"   ${NATIVE_PROJECT_SOURCE_DIR}\\..\\..\\PyDelFEM2\\eigen\\
  COMMENT "Post Event"
)
else()
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_eigen.*
  ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/eigen/
  COMMENT "Post Event"
)
endif()
