project(fem2d)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
else()
set( CMAKE_CXX_FLAGS "-Wall  -Wno-deprecated-declarations -std=c++11" )
endif()

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h               ${DELFEM2_INC}/vec2.cpp     
  ${DELFEM2_INC}/camera.h             ${DELFEM2_INC}/camera.cpp
  ${DELFEM2_INC}/color.h              ${DELFEM2_INC}/color.cpp           
  ${DELFEM2_INC}/mshmisc.h            ${DELFEM2_INC}/mshmisc.cpp              
  ${DELFEM2_INC}/mshtopo.h            ${DELFEM2_INC}/mshtopo.cpp        
  ${DELFEM2_INC}/dtri.h               ${DELFEM2_INC}/dtri.cpp         
  ${DELFEM2_INC}/mats.h               ${DELFEM2_INC}/mats.cpp    
  ${DELFEM2_INC}/emat.h               ${DELFEM2_INC}/emat.cpp     
     
  ${DELFEM2_INC}/dtri_v2.h            ${DELFEM2_INC}/dtri_v2.cpp
  ${DELFEM2_INC}/ilu_mats.h           ${DELFEM2_INC}/ilu_mats.cpp       
  ${DELFEM2_INC}/fem_emats.h          ${DELFEM2_INC}/fem_emats.cpp                      

  ${DELFEM2_INC}/opengl/gl2_color.h   ${DELFEM2_INC}/opengl/gl2_color.cpp
  ${DELFEM2_INC}/opengl/gl2_funcs.h   ${DELFEM2_INC}/opengl/gl2_funcs.cpp

  main.cpp
)

target_link_libraries(${PROJECT_NAME}
 ${GLUT_LIBRARY} 
 ${OPENGL_LIBRARY}
)
