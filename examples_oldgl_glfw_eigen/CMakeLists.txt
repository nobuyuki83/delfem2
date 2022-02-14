cmake_minimum_required(VERSION 3.12)

IF (USE_STATIC_LIB)
  message("$$$$$$$$$$$$$$$$$$$$$")
  message("generate project as STATIC library: USE_STATIC_LIB=${USE_STATIC_LIB}")
  project(examples_static)
  add_subdirectory(libdfm2_oldgl_glfw_eigen)
ELSE ()
  message("$$$$$$$$$$$$$$$$$$$$$")
  message("generate project as HEADER-ONLY library: USE_STATIC_LIB=${USE_STATIC_LIB}")
  project(examples_hdronly)
ENDIF ()

project(examples_oldgl_glfw_hdronly)

# simple examples without delfem2 
add_subdirectory(00_EigenModes)
add_subdirectory(01_ConvexHull)
add_subdirectory(02_FemSolidLinear2)
add_subdirectory(03_FemSolidLinear3)
add_subdirectory(04_ReadTriangleMesh)
add_subdirectory(05_SketchContour2)