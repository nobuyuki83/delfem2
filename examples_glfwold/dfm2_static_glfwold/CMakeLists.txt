cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -g \
    -Wno-deprecated-declarations \
    -Wno-unused-const-variable ")
ENDIF()

####################################

set(PROJ_NAME dfm2_static_glfwold)
project(${PROJ_NAME})

add_definitions(-DDFM2_STATIC_LIBRARY=ON)

# glfw
get_filename_component(GLFW_INCLUDE_DIR ../../3rd_party/glfw/include ABSOLUTE)

# glad
get_filename_component(GLAD_INCLUDE_DIR ../../3rd_party/glad3/include ABSOLUTE)

# opengl
find_package(OpenGL REQUIRED)

set(DELFEM2_INCLUDE_DIR "../../include")
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
)

set(DELFEM2_INC         "../../include/delfem2")
add_library(${PROJ_NAME} STATIC
    ${DELFEM2_INC}/vec2.h                      ${DELFEM2_INC}/vec2.cpp
    ${DELFEM2_INC}/vec3.h                      ${DELFEM2_INC}/vec3.cpp
    ${DELFEM2_INC}/mat3.h                      ${DELFEM2_INC}/mat3.cpp
    ${DELFEM2_INC}/mat4.h                      ${DELFEM2_INC}/mat4.cpp
    ${DELFEM2_INC}/quat.h                      ${DELFEM2_INC}/quat.cpp
    ${DELFEM2_INC}/camera.h                    ${DELFEM2_INC}/camera.cpp
    ${DELFEM2_INC}/primitive.h                 ${DELFEM2_INC}/primitive.cpp
    ${DELFEM2_INC}/specialfuncs.h              ${DELFEM2_INC}/specialfuncs.cpp
    ${DELFEM2_INC}/mshtopo.h                   ${DELFEM2_INC}/mshtopo.cpp
    ${DELFEM2_INC}/mshmisc.h                   ${DELFEM2_INC}/mshmisc.cpp
    ${DELFEM2_INC}/mshio.h                     ${DELFEM2_INC}/mshio.cpp
    ${DELFEM2_INC}/color.h                     ${DELFEM2_INC}/color.cpp
    ${DELFEM2_INC}/imgio.h                     ${DELFEM2_INC}/imgio.cpp
    ${DELFEM2_INC}/dtri.h                      ${DELFEM2_INC}/dtri.cpp
    ${DELFEM2_INC}/funcs.h                     ${DELFEM2_INC}/funcs.cpp
    ${DELFEM2_INC}/slice.h                     ${DELFEM2_INC}/slice.cpp
    ${DELFEM2_INC}/vecxitrsol.h                ${DELFEM2_INC}/vecxitrsol.cpp
    ${DELFEM2_INC}/emat.h                      ${DELFEM2_INC}/emat.cpp
    ${DELFEM2_INC}/mats.h                      ${DELFEM2_INC}/mats.cpp
    ${DELFEM2_INC}/iss.h                       ${DELFEM2_INC}/iss.cpp
    ${DELFEM2_INC}/adf.h                       ${DELFEM2_INC}/adf.cpp
    ${DELFEM2_INC}/bem.h                       ${DELFEM2_INC}/bem.cpp
    ${DELFEM2_INC}/bvh.h                       ${DELFEM2_INC}/bvh.cpp
    ${DELFEM2_INC}/bv.h
    ${DELFEM2_INC}/pgeo.h
    ${DELFEM2_INC}/sdf.h
        
    ${DELFEM2_INC}/geo3_v23m34q.h              ${DELFEM2_INC}/geo3_v23m34q.cpp
    ${DELFEM2_INC}/gizmo_geo3.h                ${DELFEM2_INC}/gizmo_geo3.cpp
    ${DELFEM2_INC}/fem_emats.h                 ${DELFEM2_INC}/fem_emats.cpp
    ${DELFEM2_INC}/dtri2_v2dtri.h              ${DELFEM2_INC}/dtri2_v2dtri.cpp
    ${DELFEM2_INC}/dtri3_v3dtri.h              ${DELFEM2_INC}/dtri3_v3dtri.cpp
    ${DELFEM2_INC}/dtet_v3.h                   ${DELFEM2_INC}/dtet_v3.cpp
    ${DELFEM2_INC}/ilu_mats.h                  ${DELFEM2_INC}/ilu_mats.cpp
    ${DELFEM2_INC}/cad2_dtri2.h                ${DELFEM2_INC}/cad2_dtri2.cpp
    ${DELFEM2_INC}/cad3d.h                     ${DELFEM2_INC}/cad3d.cpp
    ${DELFEM2_INC}/objf_geo3.h                 ${DELFEM2_INC}/objf_geo3.cpp
    ${DELFEM2_INC}/objfdtri_objfdtri23.h       ${DELFEM2_INC}/objfdtri_objfdtri23.cpp
    ${DELFEM2_INC}/def.h                       ${DELFEM2_INC}/def.cpp    
    ${DELFEM2_INC}/rig_geo3.h                  ${DELFEM2_INC}/rig_geo3.cpp
    ${DELFEM2_INC}/srchuni_v3.h                ${DELFEM2_INC}/srchuni_v3.cpp
    ${DELFEM2_INC}/gridvoxel.h                 ${DELFEM2_INC}/gridvoxel.cpp
    ${DELFEM2_INC}/gridcube.h                  ${DELFEM2_INC}/gridcube.cpp
    ${DELFEM2_INC}/femrod.h                    ${DELFEM2_INC}/femrod.cpp    
    ${DELFEM2_INC}/srchbi_v3bvh.h
    ${DELFEM2_INC}/srch_v3bvhmshtopo.h
    ${DELFEM2_INC}/cloth_internal.h
    ${DELFEM2_INC}/rigopt.h
    
    ${DELFEM2_INC}/opengl/tex_gl.h             ${DELFEM2_INC}/opengl/tex_gl.cpp
    ${DELFEM2_INC}/opengl/v2_glold.h           ${DELFEM2_INC}/opengl/v2_glold.cpp
    ${DELFEM2_INC}/opengl/v3q_glold.h          ${DELFEM2_INC}/opengl/v3q_glold.cpp
    ${DELFEM2_INC}/opengl/color_glold.h        ${DELFEM2_INC}/opengl/color_glold.cpp
    ${DELFEM2_INC}/opengl/funcs_glold.h        ${DELFEM2_INC}/opengl/funcs_glold.cpp
    ${DELFEM2_INC}/opengl/cad2dtriv2_glold.h   ${DELFEM2_INC}/opengl/cad2dtriv2_glold.cpp
    ${DELFEM2_INC}/opengl/caddtri_v3_glold.h   ${DELFEM2_INC}/opengl/caddtri_v3_glold.cpp
    ${DELFEM2_INC}/opengl/rigv3_glold.h        ${DELFEM2_INC}/opengl/rigv3_glold.cpp
    ${DELFEM2_INC}/opengl/colorv3_glold.h      ${DELFEM2_INC}/opengl/colorv3_glold.cpp
    ${DELFEM2_INC}/opengl/r2t_gl.h             ${DELFEM2_INC}/opengl/r2t_gl.cpp
    ${DELFEM2_INC}/opengl/r2tglo_glold.h       ${DELFEM2_INC}/opengl/r2tglo_glold.cpp
    ${DELFEM2_INC}/opengl/gridcube_glold.h     ${DELFEM2_INC}/opengl/gridcube_glold.cpp
    ${DELFEM2_INC}/opengl/glfw/viewer_glfw.h   ${DELFEM2_INC}/opengl/glfw/viewer_glfw.cpp
    ${DELFEM2_INC}/opengl/glfw/cam_glfw.h
)

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(WIN32)
    message("copy lib from: ${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\dfm2_static_glfwold.lib")
    message("copy lib to: ${NATIVE_PROJECT_SOURCE_DIR}")
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND copy
        "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\dfm2_static_glfwold.lib"
        ${NATIVE_PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
        )
else()
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/libdfm2_static_glfwold.*
        ${PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
        )
endif()
