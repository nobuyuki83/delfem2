cmake_minimum_required(VERSION 2.8)

# compile as "cmake -DEMSCRIPTEN=1 ."  to build with Emscripten
option(EMSCRIPTEN "Use Emscripten:" OFF)

####################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(NOT EMSCRIPTEN)
  IF(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
  ELSE()
    set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
  ENDIF()
ENDIF()

####################################

project(41_fem2d_poisson)


set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h                    ${DELFEM2_INC}/vec2.cpp     
  ${DELFEM2_INC}/camera.h                  ${DELFEM2_INC}/camera.cpp
  ${DELFEM2_INC}/color.h                   ${DELFEM2_INC}/color.cpp           
  ${DELFEM2_INC}/mshmisc.h                 ${DELFEM2_INC}/mshmisc.cpp              
  ${DELFEM2_INC}/mshtopo.h                 ${DELFEM2_INC}/mshtopo.cpp        
  ${DELFEM2_INC}/dtri.h                    ${DELFEM2_INC}/dtri.cpp         
  ${DELFEM2_INC}/mats.h                    ${DELFEM2_INC}/mats.cpp    
  ${DELFEM2_INC}/emat.h                    ${DELFEM2_INC}/emat.cpp
  ${DELFEM2_INC}/vecxitrsol.h              ${DELFEM2_INC}/vecxitrsol.cpp     
        
  ${DELFEM2_INC}/dtri_v2.h                 ${DELFEM2_INC}/dtri_v2.cpp
  ${DELFEM2_INC}/ilu_mats.h                ${DELFEM2_INC}/ilu_mats.cpp       
  ${DELFEM2_INC}/fem_emats.h               ${DELFEM2_INC}/fem_emats.cpp                      

  ${DELFEM2_INC}/opengl/gl_funcs.h         ${DELFEM2_INC}/opengl/gl_funcs.cpp
  ${DELFEM2_INC}/opengl/glnew_funcs.h      ${DELFEM2_INC}/opengl/glnew_funcs.cpp
  ${DELFEM2_INC}/opengl/glnew_mshcolor.h   ${DELFEM2_INC}/opengl/glnew_mshcolor.cpp
  ${DELFEM2_INC}/opengl/glfw_viewer.h      ${DELFEM2_INC}/opengl/glfw_viewer.cpp
  ${DELFEM2_INC}/opengl/glfw_cam.h

  ${3RD_PARTY}/glad4/include/glad/glad.h
  ${3RD_PARTY}/glad4/src/glad.c

  main.cpp
)

set(GLAD_INCLUDE_DIR "${3RD_PARTY}/glad4/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

if (NOT EMSCRIPTEN)
  set(GLFW_DIR "${3RD_PARTY}/glfw")
  get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)
  find_package(OpenGL REQUIRED)
  include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLAD_INCLUDE_DIR}
    ${GLFW_DIR}/include
    ${DELFEM2_INCLUDE_DIR}
  )
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
    if(CMAKE_CL_64)
      target_link_libraries(${PROJECT_NAME}
        ${OPENGL_LIBRARY}
        ${GLFW_DIR}/buildVS64/src/Release/glfw3.lib
        ${EXTRA_LIBS}
      )
    else()
        target_link_libraries(${PROJECT_NAME}
        ${OPENGL_LIBRARY}
        ${GLFW_DIR}/buildVS32/src/Release/glfw3.lib
        ${EXTRA_LIBS}
      )
    endif()
  endif()
else()
  set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE BOTH)
  set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".html")
  set(CMAKE_C_COMPILER "emcc")
  set(CMAKE_CXX_COMPILER "em++")
  target_compile_definitions(${PROJECT_NAME} PUBLIC -DEMSCRIPTEN=1)
  set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-s USE_GLFW=3 -s USE_WEBGL2=1 -s FULL_ES3=1")
#  set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-s USE_GLFW=3 -s FULL_ES3=1")
  include_directories(
    ${GLAD_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
  )
endif()

