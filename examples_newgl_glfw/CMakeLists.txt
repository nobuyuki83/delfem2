cmake_minimum_required(VERSION 3.12)

option(USE_STATIC_LIB "compile delfem2 as a static library" OFF)

if(USE_STATIC_LIB)
  message("$$$$$$$$$$$$$$$$$$$$$")
  message("generate project as STATIC library: USE_STATIC_LIB=${USE_STATIC_LIB}")
  project(examples_glfwnew_static)
  add_subdirectory(dfm2_static_glfwnew)
ELSE()
  message("$$$$$$$$$$$$$$$$$$$$$")
  message("generate project as HEADER-ONLY library:  USE_STATIC_LIB=${USE_STATIC_LIB}")
  project(examples_glfnew_hdronly)
ENDIF()

add_subdirectory(00_openwin)
add_subdirectory(01_drawrect)
add_subdirectory(02_nav3d)
add_subdirectory(03_texture)
add_subdirectory(04_render2texture)
add_subdirectory(05_OffScreenProjection)
add_subdirectory(06_VolumeRenderingSlice)
add_subdirectory(07_RayCasting)
add_subdirectory(08_SphereTracing)
add_subdirectory(09_VolumeRenderingRayCasting)
add_subdirectory(10_FlatColorTriangles)

add_subdirectory(20_cad2d)

add_subdirectory(30_cadmsh2d)

add_subdirectory(40_femcloth)
add_subdirectory(41_fem2d_poisson)
add_subdirectory(42_fem2d_linearsolid)
