cmake_minimum_required(VERSION 3.0)

# compile as "cmake -DEMSCRIPTEN=1 ."  to build with Emscripten
option(EMSCRIPTEN "Use Emscripten:" OFF)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

##########################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(NOT EMSCRIPTEN)
  IF(MSVC)
  ELSE()
    set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
  ENDIF()
ENDIF()

##########################################

project(00_openwin)

##########################################

# glfw
set(GLFW_DIR "../../3rd_party/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE) # get absolute path

# opengl
find_package(OpenGL REQUIRED)

#########################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR}
  #/Users/nobuyuki/emsdk/fastcomp/emscripten/system/include/compat
)

add_executable(${PROJECT_NAME}
  main.cpp
)

########################################

if (NOT EMSCRIPTEN)
  IF(APPLE)
    find_library( COCOA_LIBRARY Cocoa )
    find_library( OpenGL_LIBRARY OpenGL )
    find_library( IOKIT_LIBRARY IOKit )
    find_library( COREVIDEO_LIBRARY CoreVideo )
    target_link_libraries(${PROJECT_NAME} 
      ${OPENGL_LIBRARY}
      ${GLFW_DIR}/src/libglfw3.a
      ${COCOA_LIBRARY} ${OpenGL_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} 
    )
  ELSEIF(UNIX)
    target_link_libraries(${PROJECT_NAME} 
      ${OPENGL_LIBRARY}
      ${GLFW_DIR}/src/libglfw3.a
      -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp 
    )  
  ELSEIF(MSVC)
    target_link_libraries(${PROJECT_NAME}
      ${OPENGL_LIBRARY}
      ${GLFW_DIR}/src/Release/glfw3.lib
    ) 
  ENDIF()
else()
#    set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE BOTH)
#    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".html")
#    set(CMAKE_C_COMPILER "emcc")
    set(CMAKE_CXX_COMPILER "em++")
#    target_compile_definitions(${PROJECT_NAME} PUBLIC -DEMSCRIPTEN=1)
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-o index.html -s USE_GLFW=3 -s USE_WEBGL2=1 -s WASM=1 -std=c++11")
endif()
