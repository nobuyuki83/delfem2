cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

######################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
set(CMAKE_CXX_FLAGS "-O2 /W4\
  /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 \
  /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 ")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations  -std=c++11")
endif()

####################################

project(616_PbdClothCadGltf)

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
      ${DELFEM2_INC}/quat.h                      ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/camera.h               ${DELFEM2_INC}/camera.cpp
      ${DELFEM2_INC}/mat4.h                     ${DELFEM2_INC}/mat4.cpp
      )
ELSE()
  message("compile as HEADER-ONLY library: ${PROJECT_NAME} ${USE_HEADERONLY}")
ENDIF()

# add_definitions(-DDFM2_STATIC_LIBRARY)
add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h                     ${DELFEM2_INC}/vec2.cpp
  ${DELFEM2_INC}/vec3.h                     ${DELFEM2_INC}/vec3.cpp
  ${DELFEM2_INC}/mat3.h                     ${DELFEM2_INC}/mat3.cpp
  ${DELFEM2_INC}/dtri.h                     ${DELFEM2_INC}/dtri.cpp
  ${DELFEM2_INC}/mshtopo.h                  ${DELFEM2_INC}/mshtopo.cpp
  ${DELFEM2_INC}/mshmisc.h                  ${DELFEM2_INC}/mshmisc.cpp
  ${DELFEM2_INC}/bvh.h                      ${DELFEM2_INC}/bvh.cpp
  ${DELFEM2_INC}/funcs.h                    ${DELFEM2_INC}/funcs.cpp
  ${DELFEM2_INC}/vecxitrsol.h               ${DELFEM2_INC}/vecxitrsol.cpp
  ${DELFEM2_INC}/bv.h        
      
  ${DELFEM2_INC}/dtri_v2.h                  ${DELFEM2_INC}/dtri_v2.cpp
  ${DELFEM2_INC}/objfunc_v23.h              ${DELFEM2_INC}/objfunc_v23.cpp
  ${DELFEM2_INC}/objfunc_v23dtri.h          ${DELFEM2_INC}/objfunc_v23dtri.cpp
  ${DELFEM2_INC}/cad2d_v2dtri.h             ${DELFEM2_INC}/cad2d_v2dtri.cpp
  ${DELFEM2_INC}/v23m3q.h                   ${DELFEM2_INC}/v23m3q.cpp
  ${DELFEM2_INC}/srchuni_v3.h               ${DELFEM2_INC}/srchuni_v3.cpp
  ${DELFEM2_INC}/srch_v3bvhmshtopo.h
  
  ${DELFEM2_INC}/opengl/glold_funcs.h       ${DELFEM2_INC}/opengl/glold_funcs.cpp
  ${DELFEM2_INC}/opengl/v3q_glold.h         ${DELFEM2_INC}/opengl/v3q_glold.cpp
  ${DELFEM2_INC}/opengl/caddtri_v3_glold.h  ${DELFEM2_INC}/opengl/caddtri_v3_glold.cpp
  
  ${DELFEM2_INC}/rig_v3q.h                  ${DELFEM2_INC}/rig_v3q.cpp
  ${DELFEM2_INC}/opengl/rig_v3m3q_glold.h   ${DELFEM2_INC}/opengl/rig_v3m3q_glold.cpp
  
  ${DELFEM2_INC}/opengl/glfw/viewer_glfw.h  ${DELFEM2_INC}/opengl/glfw/viewer_glfw.cpp

  ${3RD_PARTY}/io_gltf.h                    ${3RD_PARTY}/io_gltf.cpp
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

