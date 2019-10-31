cmake_minimum_required(VERSION 3.0)

project(00_openwin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
add_executable(${PROJECT_NAME} main.cpp)


# compile as "cmake -DEMSCRIPTEN=1 ."  to build with Emscripten

option(EMSCRIPTEN "Use Emscripten:" OFF)
if (NOT EMSCRIPTEN)
  set(GLFW_DIR "../../src_cpp/external/glfw")
  get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE) # get absolute path
  find_package(OpenGL REQUIRED)
  include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLFW_DIR}/include
    ${DELFEM2_INCLUDE_DIR}
  )
  IF(APPLE)
    set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
    find_library( COCOA_LIBRARY Cocoa )
    find_library( OpenGL_LIBRARY OpenGL )
    find_library( IOKIT_LIBRARY IOKit )
    find_library( COREVIDEO_LIBRARY CoreVideo )
    SET( EXTRA_LIBS ${COCOA_LIBRARY} ${OpenGL_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} )
  ELSEIF(UNIX)
    set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
    SET( EXTRA_LIBS   -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp )
  ENDIF()
  target_link_libraries(${PROJECT_NAME} 
    ${OPENGL_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    ${EXTRA_LIBS}
  )
else()
    set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE BOTH)
    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".html")
    set(CMAKE_C_COMPILER "emcc")
    set(CMAKE_CXX_COMPILER "em++")    
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DEMSCRIPTEN=1)
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-s USE_GLFW=3 -s USE_WEBGL2=1 -s FULL_ES3=1")
endif()
