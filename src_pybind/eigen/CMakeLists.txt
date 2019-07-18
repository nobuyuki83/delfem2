cmake_minimum_required(VERSION 3.0)

set(PROJ_NAME c_eigen)

project(${PROJ_NAME})
set(PYBIND11_CPP_STANDARD -std=c++14)
if(WIN32)
set(CMAKE_CXX_FLAGS "-Wall -O2")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -O2")
endif()		
add_definitions(-DUSE_GL)

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

  ${DELFEM2_SRC}/gl_v23q.cpp       ${DELFEM2_INC}/gl_v23q.h

  main.cpp
)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)

set(INPUT_INCLUDE_DIR
  ${OPENGL_INCLUDE_DIR}
  ${GLEW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  "../../src_cpp/external/eigen-git-mirror"
)
set(INPUT_LIBRARY
  ${OPENGL_LIBRARY}
  ${GLEW_LIBRARY}
)


########################### start use FBX
if(USE_FBX) 
add_definitions(-DUSE_FBX)
IF(APPLE)
  MESSAGE("platform apple os x")
  SET(FBX_ROOT "/Applications/Autodesk/FBX SDK/2019.0")
  SET(FBX_INCLUDE_DIR ${INCLUDE_DIR} ${FBX_ROOT}/include)
  SET(FBX_LINK_LIBS ${FBX_ROOT}/lib/clang/release/libfbxsdk.a)
ENDIF(APPLE)
set(INPUT_SRC
  ${INPUT_SRC}
  ${DELFEM2_SRC}/../external/io_fbx.cpp  
  ${DELFEM2_INC}/../../external/io_fbx.h
  py_fbx.cpp
)
set(INPUT_INCLUDE_DIR
  ${INPUT_INCLUDE_DIR}
  ${FBX_INCLUDE_DIR}  
)
set(INPUT_LIBRARY
  ${INPUT_LIBRARY}
  ${FBX_LINK_LIBS}
)
endif()
########################### end use FBS

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
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_eigen.*
  ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/eigen/
  COMMENT "Post Event"
)
endif()
