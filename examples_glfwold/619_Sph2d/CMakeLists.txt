cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

######################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
set(CMAKE_CXX_FLAGS "-O2 /W4 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996\
    /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )
endif()

#############################

project(619_Sph2d)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")
set(GLFW_DIR "${3RD_PARTY}/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)  

find_package(OpenGL REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR}
)

IF(NOT USE_HEADERONLY)
  message("compile as STATIC libary: ${PROJECT_NAME} ${USE_HEADERONLY}")
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  set(SRC_DFM2
      ${DELFEM2_INC}/quat.h                      ${DELFEM2_INC}/quat.cpp
      ${DELFEM2_INC}/camera.h               ${DELFEM2_INC}/camera.cpp
      )
ELSE()
  message("compile as HEADER-ONLY library: ${PROJECT_NAME} ${USE_HEADERONLY}")
ENDIF()

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec3.h       ${DELFEM2_INC}/vec3.cpp
  
  ${DELFEM2_INC}/opengl/glold_funcs.h   ${DELFEM2_INC}/opengl/glold_funcs.cpp
  
  ${DELFEM2_INC}/opengl/glfw/viewer_glfw.h      ${DELFEM2_INC}/opengl/glfw/viewer_glfw.cpp
  ${SRC_DFM2}
  main.cpp
)

#####################################

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
