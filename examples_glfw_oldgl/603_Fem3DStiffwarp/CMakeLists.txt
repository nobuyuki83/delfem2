cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
set(CMAKE_CXX_FLAGS "-O2 /W4\
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996\
    /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations" )
endif()

##########################################

project(603_Fem3StiffWarp)

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

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h               ${DELFEM2_INC}/vec2.cpp              
  ${DELFEM2_INC}/vec3.h               ${DELFEM2_INC}/vec3.cpp              
  ${DELFEM2_INC}/mat3.h               ${DELFEM2_INC}/mat3.cpp              
  ${DELFEM2_INC}/quat.h               ${DELFEM2_INC}/quat.cpp
  ${DELFEM2_INC}/camera.h             ${DELFEM2_INC}/camera.cpp
  ${DELFEM2_INC}/color.h              ${DELFEM2_INC}/color.cpp              
  ${DELFEM2_INC}/mshmisc.h            ${DELFEM2_INC}/mshmisc.cpp              
  ${DELFEM2_INC}/mshtopo.h            ${DELFEM2_INC}/mshtopo.cpp                
  ${DELFEM2_INC}/dtri.h               ${DELFEM2_INC}/dtri.cpp
  ${DELFEM2_INC}/emat.h               ${DELFEM2_INC}/emat.cpp 
  ${DELFEM2_INC}/mats.h               ${DELFEM2_INC}/mats.cpp
  ${DELFEM2_INC}/vecxitrsol.h         ${DELFEM2_INC}/vecxitrsol.cpp

  ${DELFEM2_INC}/v23m3q.h             ${DELFEM2_INC}/v23m3q.cpp                
  ${DELFEM2_INC}/dtri_v2.h            ${DELFEM2_INC}/dtri_v2.cpp                
  ${DELFEM2_INC}/fem_emats.h          ${DELFEM2_INC}/fem_emats.cpp              
  ${DELFEM2_INC}/ilu_mats.h           ${DELFEM2_INC}/ilu_mats.cpp              
  
  ${DELFEM2_INC}/opengl/glold_v23.h     ${DELFEM2_INC}/opengl/glold_v23.cpp
  ${DELFEM2_INC}/opengl/glold_color.h   ${DELFEM2_INC}/opengl/glold_color.cpp
  ${DELFEM2_INC}/opengl/glold_funcs.h   ${DELFEM2_INC}/opengl/glold_funcs.cpp
  ${DELFEM2_INC}/opengl/glfw_viewer.h   ${DELFEM2_INC}/opengl/glfw_viewer.cpp                     

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
      ${GLFW_DIR}/buildVS64/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  else()
    target_link_libraries(${PROJECT_NAME} 
      ${GLFW_DIR}/buildVS32/src/Release/glfw3.lib
      ${OPENGL_LIBRARY})
  endif()
ENDIF()

