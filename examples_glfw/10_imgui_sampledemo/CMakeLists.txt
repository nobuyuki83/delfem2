cmake_minimum_required(VERSION 2.8.12)

#####################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(NOT EMSCRIPTEN)
  IF(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
  ELSE()
    set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
  ENDIF()
ENDIF()

######################################

project(10_imgui_sampledemo)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED CONFIG)

set(GLFW_DIR "../../external/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLAD)
set(IMGUI_DIR "../imgui")

set(GLAD_INCLUDE_DIR "../glad4/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${IMGUI_DIR}
  ${IMGUI_DIR}/examples
  ${GLAD_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
  
  ${IMGUI_DIR}/imgui.cpp  
  ${IMGUI_DIR}/imgui.h    
  ${IMGUI_DIR}/imgui_draw.cpp
  ${IMGUI_DIR}/imgui_widgets.cpp
  ${IMGUI_DIR}/imgui_demo.cpp
  ${IMGUI_DIR}/examples/imgui_impl_opengl3.cpp
  ${IMGUI_DIR}/examples/imgui_impl_opengl3.h
  ${IMGUI_DIR}/examples/imgui_impl_glfw.cpp
  ${IMGUI_DIR}/examples/imgui_impl_glfw.h
  
  ../glad4/include/glad/glad.h    ../glad4/src/glad.c
)

IF(APPLE)
  find_library( COCOA_LIBRARY Cocoa )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLEW_LIBRARY}
    ${GLFW_DIR}/src/libglfw3.a
    ${COCOA_LIBRARY}
    ${IOKIT_LIBRARY}
    ${COREVIDEO_LIBRARY}
  )
ELSEIF(UNIX)
  target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARY}
    ${GLEW_LIBRARY}
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
ENDIF()
