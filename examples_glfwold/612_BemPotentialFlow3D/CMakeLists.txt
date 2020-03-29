cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

##########################################

project(612_BemPotentialFlow3D)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(GLFW_DIR "${3RD_PARTY}/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

add_definitions(-DUSE_MKL)

set(MKL_ROOT
	C:/MKL
	/opt/intel/mkl
	$ENV{MKLROOT})

find_path(MKL_INCLUDE_DIR mkl.h 
	PATHS ${MKL_ROOT} 
	PATH_SUFFIXES include)

find_library(MKL_LIBRARY_CORE
  libmkl_core.a
  PATHS ${MKL_ROOT}/lib ${MKL_ROOT}/lib/intel64)
get_filename_component(MKL_LIBRARY 
  ${MKL_LIBRARY_CORE} DIRECTORY)

find_package(OpenGL REQUIRED)

#######################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR})

add_executable(${PROJECT_NAME}
		${DELFEM2_INC}/vec3.cpp        ${DELFEM2_INC}/vec3.h
		${DELFEM2_INC}/vec2.cpp        ${DELFEM2_INC}/vec2.h
		${DELFEM2_INC}/mat3.cpp        ${DELFEM2_INC}/mat3.h
    ${DELFEM2_INC}/mat4.cpp        ${DELFEM2_INC}/mat4.h
    ${DELFEM2_INC}/quat.cpp        ${DELFEM2_INC}/quat.h
    ${DELFEM2_INC}/bem.cpp         ${DELFEM2_INC}/bem.h
    ${DELFEM2_INC}/primitive.cpp   ${DELFEM2_INC}/primitive.h
    ${DELFEM2_INC}/camera.h        ${DELFEM2_INC}/camera.cpp
    
		${DELFEM2_INC}/v23m3q.cpp      ${DELFEM2_INC}/v23m3q.h
  
    ${DELFEM2_INC}/opengl/glold_funcs.h       ${DELFEM2_INC}/opengl/glold_funcs.cpp
    ${DELFEM2_INC}/opengl/v3q_glold.h         ${DELFEM2_INC}/opengl/v3q_glold.cpp
#    ${DELFEM2_INC}/opengl/glold_v23dtricad.h  ${DELFEM2_INC}/opengl/glold_v23dtricad.cpp
    
    ${DELFEM2_INC}/opengl/glfw/viewer_glfw.h      ${DELFEM2_INC}/opengl/glfw/viewer_glfw.cpp
    
		main.cpp
)

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
	-Wl--start-group
	${MKL_LIBRARY}/libmkl_core.a
	${MKL_LIBRARY}/libmkl_intel_lp64.a
	${MKL_LIBRARY}/libmkl_sequential.a
	-Wl--end-group
	-lpthread
	-lm
	-ldl    
  )
ELSEIF(UNIX)
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
  )
ELSEIF(MSVC)
  if(CMAKE_CL_64)
    target_link_libraries(${PROJECT_NAME} 
      ${GLFW_DIR}/buildVS64/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  else()
    target_link_libraries(${PROJECT_NAME} 
      ${GLFW_DIR}/buildVS32/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  endif()
ENDIF()
