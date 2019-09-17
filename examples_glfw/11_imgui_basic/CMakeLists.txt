project(11_imgui_basic)
cmake_minimum_required(VERSION 2.8.12)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED CONFIG)
find_package(GLEW REQUIRED)

set(GLFW_DIR "../glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLEW)

set(IMGUI_DIR "../imgui")

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLEW_INCLUDE_DIR}    
  ${GLFW_DIR}/include
  ${IMGUI_DIR}
  ${IMGUI_DIR}/examples
)

add_executable(${PROJECT_NAME}
  main.cpp
  file_manager.cpp
  file_manager.h
  opengl_shader.cpp
  opengl_shader.h
  ${IMGUI_DIR}/imgui.cpp  
  ${IMGUI_DIR}/imgui.h    
  ${IMGUI_DIR}/imgui_draw.cpp
  ${IMGUI_DIR}/imgui_widgets.cpp
  ${IMGUI_DIR}/examples/imgui_impl_opengl3.cpp
  ${IMGUI_DIR}/examples/imgui_impl_opengl3.h
  ${IMGUI_DIR}/examples/imgui_impl_glfw.cpp
  ${IMGUI_DIR}/examples/imgui_impl_glfw.h
)

IF(APPLE)
  find_library( COCOA_LIBRARY Cocoa )
  find_library( OpenGL_LIBRARY OpenGL )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  SET( EXTRA_LIBS ${COCOA_LIBRARY} ${OpenGL_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} )
ELSEIF(UNIX)
  SET( EXTRA_LIBS   -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp )
ENDIF()

if(MSVC)
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
else()
  target_link_libraries(${PROJECT_NAME} 
    ${OPENGL_LIBRARY}
    ${GLEW_LIBRARY}  
    ${GLFW_DIR}/src/libglfw3.a
    ${EXTRA_LIBS}
  )
endif()
