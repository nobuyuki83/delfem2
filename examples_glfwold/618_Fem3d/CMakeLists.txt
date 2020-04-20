cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

######################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
set(CMAKE_CXX_FLAGS "-O2 -std=c++11 -W4\
    /wd4100 /wd4458 /wd4577 /wd4267\
    /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505\
    /wd4838 /wd4305 /wd4505 /wd4456")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -Wno-unused-variable  -Wno-unused-function" )
endif()

####################################

project(618_Fem3d)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")
set(GLFW_DIR "${3RD_PARTY}/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)  

find_package(OpenGL REQUIRED)
include_directories(
	${OPENGL_INCLUDE_DIR}
	${DELFEM2_INCLUDE_DIR}
  ${GLFW_DIR}/include    
)

IF(USE_HEADERONLY)
  add_definitions(-DDFM2_HEADER_ONLY=ON)            
else() 
  if(MSVC)
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/dfm2_static_glfwold.lib )
  else()
    set(DFM2_LIBRARY
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/libdfm2_static_glfwold.a )
  endif()
ENDIF()

add_executable(${PROJECT_NAME}
  main.cpp
)


###############################################

IF(APPLE)
  find_library( COCOA_LIBRARY Cocoa )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    ${COCOA_LIBRARY}
    ${IOKIT_LIBRARY}
    ${COREVIDEO_LIBRARY}
  )
ELSEIF(UNIX)
  target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
  )
ELSEIF(MSVC)
  target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${GLFW_DIR}/src/Release/glfw3.lib
    ${OPENGL_LIBRARY}
  )
ENDIF()
