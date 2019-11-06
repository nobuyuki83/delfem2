cmake_minimum_required(VERSION 2.8)

project(examples_glfw_oldgl)


# simple examples without delfem2 
add_subdirectory(000_OpenWindowWithGLFW)


# basic funcs
add_subdirectory(100_ViewNavigationIn3D)
add_subdirectory(101_subdiv)
add_subdirectory(102_convhull)
add_subdirectory(103_noise2d)
add_subdirectory(104_noise3d)
add_subdirectory(105_meshslice)
add_subdirectory(106_sphericalharmonics)
add_subdirectory(107_paramgeo2d)
add_subdirectory(108_paramgeo3d)
add_subdirectory(109_ContourLine2D)


# meshing
add_subdirectory(201_dyntri2d_triangulation)
add_subdirectory(202_dyntri2d_remesh)
add_subdirectory(203_dyntri3d_edgecollapse)
add_subdirectory(204_dyntet_tetrahedralization)


# immplicit shape
add_subdirectory(300_AdaptiveDistanceField)
add_subdirectory(301_IsoSurfaceStuffing)


# skeleton & rigging
add_subdirectory(400_ReadFileBVH)
add_subdirectory(401_ReadFileGLTF)


# cad
add_subdirectory(500_2DCad)
add_subdirectory(501_2DCadEdit)
add_subdirectory(502_2DCadMeshEdit)
add_subdirectory(503_3DCad)