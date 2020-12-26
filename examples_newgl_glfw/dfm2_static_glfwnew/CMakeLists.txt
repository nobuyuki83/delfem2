cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

####################################

set(PROJ_NAME dfm2_static_glfwnew)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")

set(GLFW_DIR "${3RD_PARTY}/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

add_definitions(-DDFM2_STATIC_LIBRARY=ON)

find_package(OpenGL REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR}
  ${3RD_PARTY}/glad4/include
)

add_library(${PROJ_NAME} STATIC
    ${DELFEM2_INC}/vec2.h                      ${DELFEM2_INC}/vec2.cpp
    ${DELFEM2_INC}/vec3.h                      ${DELFEM2_INC}/vec3.cpp
    ${DELFEM2_INC}/mat3.h                      ${DELFEM2_INC}/mat3.cpp
    ${DELFEM2_INC}/mat4.h                      ${DELFEM2_INC}/mat4.cpp
    ${DELFEM2_INC}/quat.h                      ${DELFEM2_INC}/quat.cpp
    ${DELFEM2_INC}/cam3_m4q.h                  ${DELFEM2_INC}/cam3_m4q.cpp
    ${DELFEM2_INC}/geo3_v23m34q.h              ${DELFEM2_INC}/geo3_v23m34q.cpp

    ${DELFEM2_INC}/specialfuncs.h              ${DELFEM2_INC}/specialfuncs.cpp
    ${DELFEM2_INC}/jagarray.h                  ${DELFEM2_INC}/jagarray.cpp

    ${DELFEM2_INC}/points.h                    ${DELFEM2_INC}/points.cpp
    ${DELFEM2_INC}/mshuni.h                    ${DELFEM2_INC}/mshuni.cpp
    ${DELFEM2_INC}/mshmisc.h                   ${DELFEM2_INC}/mshmisc.cpp
    ${DELFEM2_INC}/mshio.h                     ${DELFEM2_INC}/mshio.cpp
    ${DELFEM2_INC}/mshprimitive.h              ${DELFEM2_INC}/mshprimitive.cpp

    ${DELFEM2_INC}/color.h                     ${DELFEM2_INC}/color.cpp
    ${DELFEM2_INC}/imgio.h                     ${DELFEM2_INC}/imgio.cpp
    ${DELFEM2_INC}/str.h                       ${DELFEM2_INC}/str.cpp
    ${DELFEM2_INC}/file.h                      ${DELFEM2_INC}/file.cpp
    ${DELFEM2_INC}/filenpy_str.h               ${DELFEM2_INC}/filenpy_str.cpp
    ${DELFEM2_INC}/slice.h                     ${DELFEM2_INC}/slice.cpp
    ${DELFEM2_INC}/iss.h                       ${DELFEM2_INC}/iss.cpp
    ${DELFEM2_INC}/adf.h                       ${DELFEM2_INC}/adf.cpp

    ${DELFEM2_INC}/srchbvh.h                   ${DELFEM2_INC}/srchbvh.cpp
    ${DELFEM2_INC}/srchbv3sphere.h
    ${DELFEM2_INC}/srchbv3aabb.h

    ${DELFEM2_INC}/pgeo.h
    
    ${DELFEM2_INC}/lsmats.h                    ${DELFEM2_INC}/lsmats.cpp
    ${DELFEM2_INC}/lsilu_mats.h                ${DELFEM2_INC}/lsilu_mats.cpp
    ${DELFEM2_INC}/vecxitrsol.h                ${DELFEM2_INC}/vecxitrsol.cpp

    ${DELFEM2_INC}/fempoisson.h                ${DELFEM2_INC}/fempoisson.cpp
    ${DELFEM2_INC}/femdkt.h                    ${DELFEM2_INC}/femdkt.cpp
    ${DELFEM2_INC}/femcloth.h                  ${DELFEM2_INC}/femcloth.cpp
    ${DELFEM2_INC}/femhelmholtz.h              ${DELFEM2_INC}/femhelmholtz.cpp
    ${DELFEM2_INC}/femsolidlinear.h            ${DELFEM2_INC}/femsolidlinear.cpp
    ${DELFEM2_INC}/femstokes.h                 ${DELFEM2_INC}/femstokes.cpp
    ${DELFEM2_INC}/femcloth.h                  ${DELFEM2_INC}/femcloth.cpp
    ${DELFEM2_INC}/femnavierstokes.h           ${DELFEM2_INC}/femnavierstokes.cpp
    ${DELFEM2_INC}/femutil.h                   ${DELFEM2_INC}/femutil.cpp
    ${DELFEM2_INC}/femmips_geo3.h              ${DELFEM2_INC}/femmips_geo3.cpp

    ${DELFEM2_INC}/pbd_geo3dtri23.h            ${DELFEM2_INC}/pbd_geo3dtri23.cpp
    ${DELFEM2_INC}/pbd_geo3.h                  ${DELFEM2_INC}/pbd_geo3.cpp

    ${DELFEM2_INC}/dtri.h                      ${DELFEM2_INC}/dtri.cpp
    ${DELFEM2_INC}/dtri2_v2dtri.h              ${DELFEM2_INC}/dtri2_v2dtri.cpp
    ${DELFEM2_INC}/dtri3_v3dtri.h              ${DELFEM2_INC}/dtri3_v3dtri.cpp
    ${DELFEM2_INC}/dtet_v3.h                   ${DELFEM2_INC}/dtet_v3.cpp
    ${DELFEM2_INC}/cad2_dtri2.h                ${DELFEM2_INC}/cad2_dtri2.cpp

    ${DELFEM2_INC}/defarapenergy_geo3.h        ${DELFEM2_INC}/defarapenergy_geo3.cpp
    ${DELFEM2_INC}/deflap.h                    ${DELFEM2_INC}/deflap.cpp
    ${DELFEM2_INC}/defarap.h                   ${DELFEM2_INC}/defarap.cpp

    ${DELFEM2_INC}/rig_geo3.h                  ${DELFEM2_INC}/rig_geo3.cpp
    
    ${DELFEM2_INC}/srchuni_v3.h                ${DELFEM2_INC}/srchuni_v3.cpp
    ${DELFEM2_INC}/srchbi_v3bvh.h
    ${DELFEM2_INC}/srch_v3bvhmshtopo.h
    ${DELFEM2_INC}/cloth_internal.h
    ${DELFEM2_INC}/rigopt.h
    
    ${DELFEM2_INC}/opengl/tex.h             ${DELFEM2_INC}/opengl/tex.cpp
    ${DELFEM2_INC}/opengl/funcs.h           ${DELFEM2_INC}/opengl/funcs.cpp
    ${DELFEM2_INC}/opengl/r2t.h             ${DELFEM2_INC}/opengl/r2t.cpp

    ${DELFEM2_INC}/opengl/new/funcs.h        ${DELFEM2_INC}/opengl/new/funcs.cpp
    ${DELFEM2_INC}/opengl/new/mshcolor.h     ${DELFEM2_INC}/opengl/new/mshcolor.cpp
    ${DELFEM2_INC}/opengl/new/r2tgln.h       ${DELFEM2_INC}/opengl/new/r2tgln.cpp
    ${DELFEM2_INC}/opengl/new/v23dtricad.h   ${DELFEM2_INC}/opengl/new/v23dtricad.cpp

    ${DELFEM2_INC}/opengl/glfw/viewer_glfw.h   ${DELFEM2_INC}/opengl/glfw/viewer_glfw.cpp
)

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(WIN32)
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND copy
        ${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\dfm2_static_glfwnew.lib
        ${NATIVE_PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
    )
else()
    add_custom_command(TARGET ${PROJ_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/libdfm2_static_glfwnew.*
        ${PROJECT_SOURCE_DIR}
        COMMENT "Post Event"
    )
endif()
