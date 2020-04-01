cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

####################################

project(303_RigReadGltf)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(GLFW_DIR "${3RD_PARTY}/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

find_package(OpenGL REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}
)

IF(NOT USE_HEADERONLY)
  message("compile as STATIC libary: ${PROJECT_NAME} ${USE_HEADERONLY}")
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  set(SRC_DFM2
      ${DELFEM2_INC}/quat.h                 ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/camera.h               ${DELFEM2_INC}/camera.cpp
      ${DELFEM2_INC}/mat4.h                 ${DELFEM2_INC}/mat4.cpp
      ${DELFEM2_INC}/color.h                     ${DELFEM2_INC}/color.cpp
      ${DELFEM2_INC}/opengl/glfw/viewer_glfw.h   ${DELFEM2_INC}/opengl/glfw/viewer_glfw.cpp
      )
ELSE()
  message("compile as HEADER-ONLY library: ${PROJECT_NAME} ${USE_HEADERONLY}")
ENDIF()

add_executable(${PROJECT_NAME}
 ${DELFEM2_INC}/vec2.h                       ${DELFEM2_INC}/vec2.cpp
  ${DELFEM2_INC}/vec3.h                      ${DELFEM2_INC}/vec3.cpp
  ${DELFEM2_INC}/mat3.h                      ${DELFEM2_INC}/mat3.cpp
  ${DELFEM2_INC}/mshmisc.h                   ${DELFEM2_INC}/mshmisc.cpp
  ${DELFEM2_INC}/mshio.h                     ${DELFEM2_INC}/mshio.cpp
  ${DELFEM2_INC}/funcs.h                     ${DELFEM2_INC}/funcs.cpp
               
  ${DELFEM2_INC}/v23m3q.h                    ${DELFEM2_INC}/v23m3q.cpp
  
  ${DELFEM2_INC}/opengl/glold_funcs.h        ${DELFEM2_INC}/opengl/glold_funcs.cpp
  ${DELFEM2_INC}/opengl/glold_color.h        ${DELFEM2_INC}/opengl/glold_color.cpp
  ${DELFEM2_INC}/opengl/v3q_glold.h          ${DELFEM2_INC}/opengl/v3q_glold.cpp
  
  ${DELFEM2_INC}/rig_v3q.h                   ${DELFEM2_INC}/rig_v3q.cpp
  ${DELFEM2_INC}/opengl/rig_v3m3q_glold.h    ${DELFEM2_INC}/opengl/rig_v3m3q_glold.cpp
  

  ${3RD_PARTY}/io_gltf.h
  ${3RD_PARTY}/io_gltf.cpp
  ${SRC_DFM2}

  main.cpp
)

###############################################

IF(APPLE)
  find_library( COCOA_LIBRARY Cocoa )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    ${COCOA_LIBRARY}
    ${IOKIT_LIBRARY}
    ${COREVIDEO_LIBRARY}
  )
ELSEIF(UNIX)
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
  )
ELSEIF(MSVC)
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    target_link_libraries(${PROJECT_NAME} 
      ${GLFW_DIR}/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  else()
    target_link_libraries(${PROJECT_NAME} 
      ${GLFW_DIR}/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  endif()
ENDIF()




