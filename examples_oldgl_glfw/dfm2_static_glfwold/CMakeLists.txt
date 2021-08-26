cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -g")
ENDIF()

####################################

set(PROJ_NAME dfm2_static_glfwold)
project(${PROJ_NAME})

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# glad
get_filename_component(GLAD_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad3/include ABSOLUTE)

# opengl
find_package(OpenGL REQUIRED)

add_definitions(-DDFM2_STATIC_LIBRARY=ON)
set(DELFEM2_INCLUDE_DIR "../../include")

#####################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_INCLUDE_DIR}
  ${DELFEM2_INCLUDE_DIR}
  ${GLAD_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party
)

set(DELFEM2_INC         "../../include/delfem2")
add_library(${PROJ_NAME} STATIC
    ${DELFEM2_INC}/mat2.h                      ${DELFEM2_INC}/mat2.cpp
    ${DELFEM2_INC}/vec2.h                      ${DELFEM2_INC}/vec2.cpp
    ${DELFEM2_INC}/vec3.h                      ${DELFEM2_INC}/vec3.cpp
    ${DELFEM2_INC}/mat3.h                      ${DELFEM2_INC}/mat3.cpp
    ${DELFEM2_INC}/mat4.h                      ${DELFEM2_INC}/mat4.cpp
    ${DELFEM2_INC}/quat.h                      ${DELFEM2_INC}/quat.cpp
    ${DELFEM2_INC}/cam3_m4q.h                  ${DELFEM2_INC}/cam3_m4q.cpp
    ${DELFEM2_INC}/geo3_v23m34q.h              ${DELFEM2_INC}/geo3_v23m34q.cpp
    ${DELFEM2_INC}/gizmo_geo3.h                ${DELFEM2_INC}/gizmo_geo3.cpp
    ${DELFEM2_INC}/geoplygn2_v2.h              ${DELFEM2_INC}/geoplygn2_v2.cpp
    ${DELFEM2_INC}/geoproximity3_v3.h          ${DELFEM2_INC}/geoproximity3_v3.cpp
    ${DELFEM2_INC}/geodelaunay3_v3.h           ${DELFEM2_INC}/geodelaunay3_v3.cpp
    ${DELFEM2_INC}/geosolidelm_v3.h            ${DELFEM2_INC}/geosolidelm_v3.cpp
    ${DELFEM2_INC}/geoconvhull3.h 

    ${DELFEM2_INC}/specialfuncs.h              ${DELFEM2_INC}/specialfuncs.cpp
    ${DELFEM2_INC}/jagarray.h                  ${DELFEM2_INC}/jagarray.cpp
    
    ${DELFEM2_INC}/mshuni.h                    ${DELFEM2_INC}/mshuni.cpp
    ${DELFEM2_INC}/mshmix.h                    ${DELFEM2_INC}/mshmix.cpp
    ${DELFEM2_INC}/mshsubdiv.h                 ${DELFEM2_INC}/mshsubdiv.cpp
    ${DELFEM2_INC}/mshmisc.h                   ${DELFEM2_INC}/mshmisc.cpp
    ${DELFEM2_INC}/points.h                    ${DELFEM2_INC}/points.cpp
    ${DELFEM2_INC}/mshio.h                     ${DELFEM2_INC}/mshio.cpp
    ${DELFEM2_INC}/mshprimitive.h              ${DELFEM2_INC}/mshprimitive.cpp
    ${DELFEM2_INC}/slice.h                     ${DELFEM2_INC}/slice.cpp
    ${DELFEM2_INC}/mshtopoio.h                 ${DELFEM2_INC}/mshtopoio.cpp

    ${DELFEM2_INC}/imgio.h                     ${DELFEM2_INC}/imgio.cpp
    ${DELFEM2_INC}/file.h                      ${DELFEM2_INC}/file.cpp
    ${DELFEM2_INC}/filenpy_str.h               ${DELFEM2_INC}/filenpy_str.cpp
    ${DELFEM2_INC}/str.h                       ${DELFEM2_INC}/str.cpp

    ${DELFEM2_INC}/color.h                     ${DELFEM2_INC}/color.cpp

    ${DELFEM2_INC}/femutil.h                   ${DELFEM2_INC}/femutil.cpp
    ${DELFEM2_INC}/femhelmholtz.h              ${DELFEM2_INC}/femhelmholtz.cpp
    ${DELFEM2_INC}/femmitc3.h                  ${DELFEM2_INC}/femmitc3.cpp
    ${DELFEM2_INC}/femrod.h                    ${DELFEM2_INC}/femrod.cpp
    ${DELFEM2_INC}/femcloth.h                  ${DELFEM2_INC}/femcloth.cpp
    ${DELFEM2_INC}/fempoisson.h                ${DELFEM2_INC}/fempoisson.cpp
    ${DELFEM2_INC}/femsolidlinear.h            ${DELFEM2_INC}/femsolidlinear.cpp
    ${DELFEM2_INC}/femstokes.h                 ${DELFEM2_INC}/femstokes.cpp
    ${DELFEM2_INC}/femnavierstokes.h           ${DELFEM2_INC}/femnavierstokes.cpp
    ${DELFEM2_INC}/femmips_geo3.h              ${DELFEM2_INC}/femmips_geo3.cpp
    ${DELFEM2_INC}/femsolidhyper.h             ${DELFEM2_INC}/femsolidhyper.cpp
    ${DELFEM2_INC}/cloth_internal.h

    ${DELFEM2_INC}/pbd_geo3.h                  ${DELFEM2_INC}/pbd_geo3.cpp
    ${DELFEM2_INC}/pbd_geo3dtri23.h            ${DELFEM2_INC}/pbd_geo3dtri23.cpp

    ${DELFEM2_INC}/rgd_v2m3.h                  ${DELFEM2_INC}/rgd_v2m3.cpp

    ${DELFEM2_INC}/isrf_iss.h                  ${DELFEM2_INC}/isrf_iss.cpp
    ${DELFEM2_INC}/isrf_adf.h                  ${DELFEM2_INC}/isrf_adf.cpp
    ${DELFEM2_INC}/isrf_sdf.h

    ${DELFEM2_INC}/bem.h                       ${DELFEM2_INC}/bem.cpp

    ${DELFEM2_INC}/srchbvh.h                   ${DELFEM2_INC}/srchbvh.cpp
    ${DELFEM2_INC}/srchbv3sphere.h
    ${DELFEM2_INC}/srchbv3aabb.h
    ${DELFEM2_INC}/srchuni_v3.h                ${DELFEM2_INC}/srchuni_v3.cpp
    ${DELFEM2_INC}/srchbi_v3bvh.h
    ${DELFEM2_INC}/srch_v3bvhmshtopo.h

    ${DELFEM2_INC}/dtri.h                      ${DELFEM2_INC}/dtri.cpp
    ${DELFEM2_INC}/pgeo.h
    ${DELFEM2_INC}/dtri2_v2dtri.h              ${DELFEM2_INC}/dtri2_v2dtri.cpp
    ${DELFEM2_INC}/dtri3_v3dtri.h              ${DELFEM2_INC}/dtri3_v3dtri.cpp
    ${DELFEM2_INC}/dtet_v3.h                   ${DELFEM2_INC}/dtet_v3.cpp
    ${DELFEM2_INC}/cad2_dtri2.h                ${DELFEM2_INC}/cad2_dtri2.cpp
    ${DELFEM2_INC}/cad3d.h                     ${DELFEM2_INC}/cad3d.cpp
    
    ${DELFEM2_INC}/vecxitrsol.h                ${DELFEM2_INC}/vecxitrsol.cpp
    ${DELFEM2_INC}/lsmats.h                    ${DELFEM2_INC}/lsmats.cpp
    ${DELFEM2_INC}/vecxitrsol.h                ${DELFEM2_INC}/vecxitrsol.cpp
    ${DELFEM2_INC}/lsilu_mats.h                ${DELFEM2_INC}/lsilu_mats.cpp
    ${DELFEM2_INC}/ls_masterslave.h            ${DELFEM2_INC}/ls_masterslave.cpp
    ${DELFEM2_INC}/lsitrsol.h
    ${DELFEM2_INC}/lsvecx.h

    ${DELFEM2_INC}/deflap.h                    ${DELFEM2_INC}/deflap.cpp    
    ${DELFEM2_INC}/defarap.h                   ${DELFEM2_INC}/defarap.cpp
    ${DELFEM2_INC}/gridvoxel.h                 ${DELFEM2_INC}/gridvoxel.cpp
    ${DELFEM2_INC}/gridcube.h                  ${DELFEM2_INC}/gridcube.cpp

    ${DELFEM2_INC}/clusterpoints.h             ${DELFEM2_INC}/clusterpoints.cpp

    ${DELFEM2_INC}/rig_geo3.h                  ${DELFEM2_INC}/rig_geo3.cpp
    ${DELFEM2_INC}/rigopt.h
    
    ${DELFEM2_INC}/opengl/funcs.h           ${DELFEM2_INC}/opengl/funcs.cpp
    ${DELFEM2_INC}/opengl/tex.h             ${DELFEM2_INC}/opengl/tex.cpp
    ${DELFEM2_INC}/opengl/r2t.h             ${DELFEM2_INC}/opengl/r2t.cpp
    
    ${DELFEM2_INC}/opengl/old/v2.h           ${DELFEM2_INC}/opengl/old/v2.cpp
    ${DELFEM2_INC}/opengl/old/v3q.h          ${DELFEM2_INC}/opengl/old/v3q.cpp
    ${DELFEM2_INC}/opengl/old/color.h        ${DELFEM2_INC}/opengl/old/color.cpp
    ${DELFEM2_INC}/opengl/old/funcs.h        ${DELFEM2_INC}/opengl/old/funcs.cpp
    ${DELFEM2_INC}/opengl/old/mshuni.h       ${DELFEM2_INC}/opengl/old/mshuni.cpp
    ${DELFEM2_INC}/opengl/old/cad2dtriv2.h   ${DELFEM2_INC}/opengl/old/cad2dtriv2.cpp
    ${DELFEM2_INC}/opengl/old/caddtri_v3.h   ${DELFEM2_INC}/opengl/old/caddtri_v3.cpp
    ${DELFEM2_INC}/opengl/old/rigv3.h        ${DELFEM2_INC}/opengl/old/rigv3.cpp
    ${DELFEM2_INC}/opengl/old/colorv3.h      ${DELFEM2_INC}/opengl/old/colorv3.cpp
    ${DELFEM2_INC}/opengl/old/r2tglo.h       ${DELFEM2_INC}/opengl/old/r2tglo.cpp
    ${DELFEM2_INC}/opengl/old/gridcube.h     ${DELFEM2_INC}/opengl/old/gridcube.cpp
    ${DELFEM2_INC}/opengl/old/gizmo.h        ${DELFEM2_INC}/opengl/old/gizmo.cpp
    
    ${DELFEM2_INC}/glfw/util.h     
    ${DELFEM2_INC}/glfw/viewer2.h     
    ${DELFEM2_INC}/glfw/viewer2.cpp
    ${DELFEM2_INC}/glfw/viewer3.h     
    ${DELFEM2_INC}/glfw/viewer3.cpp

    ${DELFEM2_INC}/openglstb/glyph.h         ${DELFEM2_INC}/openglstb/glyph.cpp
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
