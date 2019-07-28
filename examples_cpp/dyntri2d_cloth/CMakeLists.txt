project(dyntri2d_cloth)
cmake_minimum_required(VERSION 2.8)

if(MSVC)
set(CMAKE_CXX_FLAGS "/W4 /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456 -O2 -std=c++11")
else()
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
endif()

set(DELFEM2_INCLUDE_DIR "../../src_cpp/include")
set(DELFEM2_SRC         "../../src_cpp/cpp")
set(DELFEM2_INC         "../../src_cpp/include/delfem2")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  ${DELFEM2_INC}/vec2.h              ${DELFEM2_SRC}/vec2.cpp                                
  ${DELFEM2_INC}/dyntri.h            ${DELFEM2_SRC}/dyntri.cpp                                
  ${DELFEM2_INC}/mshtopo.h           ${DELFEM2_SRC}/mshtopo.cpp                                
  ${DELFEM2_INC}/fem_ematrix.h       ${DELFEM2_SRC}/fem_ematrix.cpp                                
  ${DELFEM2_INC}/matrix_sparse.h     ${DELFEM2_SRC}/matrix_sparse.cpp                                
  ${DELFEM2_INC}/ilu_sparse.h        ${DELFEM2_SRC}/ilu_sparse.cpp                                
  ${DELFEM2_INC}/fem.h               ${DELFEM2_SRC}/fem.cpp                                      
      
  ${DELFEM2_INC}/dyntri_v2.h         ${DELFEM2_SRC}/dyntri_v2.cpp                                
                                     
  ${DELFEM2_INC}/cloth_internal.h     
  ${DELFEM2_INC}/gl_funcs.h          ${DELFEM2_SRC}/gl_funcs.cpp

  main.cpp  
)

target_link_libraries(${PROJECT_NAME}  
 ${GLUT_LIBRARY} 
 ${OPENGL_LIBRARY}
)

