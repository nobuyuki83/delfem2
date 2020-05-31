cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

#################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
  set(CMAKE_CXX_FLAGS "-O2 /W4 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 \
    /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )
endif()

#################################

project(626_RdgRotation)

# opengl
find_package(OpenGL REQUIRED)

# glfw
IF(MSVC)
  set(GLFW_LIBRARY    ../../3rd_party/glfw/src/Release/glfw3.lib)
ELSE()
  set(GLFW_LIBRARY    ../../3rd_party/glfw/src/libglfw3.a)
ENDIF()
get_filename_component(GLFW_LIBRARY ${GLFW_LIBRARY} ABSOLUTE)
get_filename_component(GLFW_INCLUDE_DIR ../../3rd_party/glfw/include ABSOLUTE)

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

#################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
)

#####################################

IF(APPLE)
  find_library( COCOA_LIBRARY Cocoa )
  find_library( IOKIT_LIBRARY IOKit )
  find_library( COREVIDEO_LIBRARY CoreVideo )
  target_link_libraries(${PROJECT_NAME}
      ${OPENGL_LIBRARY}
      ${GLFW_LIBRARY}
      ${COCOA_LIBRARY}
      ${IOKIT_LIBRARY}
      ${COREVIDEO_LIBRARY}
      ${DFM2_LIBRARY}
      )
ELSEIF(UNIX)
  target_link_libraries(${PROJECT_NAME}
      ${OPENGL_LIBRARY}
      ${GLFW_LIBRARY}
      ${DFM2_LIBRARY}
      -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
      )
ELSEIF(MSVC)
  target_link_libraries(${PROJECT_NAME}
      ${GLFW_LIBRARY}
      ${OPENGL_LIBRARY}
      ${DFM2_LIBRARY})
ENDIF()
