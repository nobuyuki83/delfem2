cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -g\
  -Wno-unused-variable \
  -Wno-deprecated-declarations" )
ENDIF()

####################################

project(20_PoseBoneBlendshape)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
IF(MSVC)
  set(GLFW_LIBRARY    ../../3rd_party/glfw/src/Release/glfw3.lib)
ELSE()
  set(GLFW_LIBRARY    ../../3rd_party/glfw/src/libglfw3.a)
ENDIF()
get_filename_component(GLFW_LIBRARY ${GLFW_LIBRARY} ABSOLUTE)
get_filename_component(GLFW_INCLUDE_DIR ../../3rd_party/glfw/include ABSOLUTE)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
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
set(DELFEM2_INCLUDE_DIR "../../include")

# zlib
find_package(ZLIB REQUIRED)
message(STATUS "ZLIB_INCLUDE_DIRS:" ${ZLIB_INCLUDE_DIRS})
message(STATUS "ZLIB_LIBRARIES:" ${ZLIB_LIBRARIES})

#########################################

set(3RD_PARTY           "../../3rd_party")
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}
  ${ZLIB_INCLUDE_DIR}
)

set(DELFEM2_INC         "../../include/delfem2")
add_executable(${PROJECT_NAME}
  ${3RD_PARTY}/cnpy/cnpy.h                  
  ${3RD_PARTY}/cnpy/cnpy.cpp
  main.cpp
)

################################################

IF(APPLE)
  find_library( COCOA_LIBRARY Cocoa )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  target_link_libraries(${PROJECT_NAME}
      ${DFM2_LIBRARY}
      ${OPENGL_LIBRARY}
      ${GLFW_LIBRARY}
      ${COCOA_LIBRARY}
      ${IOKIT_LIBRARY}
      ${COREVIDEO_LIBRARY}
      ${ZLIB_LIBRARIES}
      )
ELSEIF(UNIX)
  target_link_libraries(${PROJECT_NAME}
      ${DFM2_LIBRARY}
      ${OPENGL_LIBRARY}
      ${GLFW_LIBRARY}
      -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp -lz
  )
ELSEIF(MSVC)
  target_link_libraries(${PROJECT_NAME}
      ${DFM2_LIBRARY}
      ${GLFW_LIBRARY}
      ${OPENGL_LIBRARY}
      ${ZLIB_LIBRARIES}
  )  
ENDIF()
