cmake_minimum_required(VERSION 2.8)

project(examples_glfw_oldgl)


# simple examples without delfem2 
add_subdirectory(000_OpenWindowWithGLFW)


# basic funcs
add_subdirectory(100_ViewNavigationIn3D)
add_subdirectory(101_SubdivisionSurfaceCatmull)
add_subdirectory(102_ConvexHull)
add_subdirectory(103_noise2d)
add_subdirectory(104_noise3d)
add_subdirectory(105_meshslice)
add_subdirectory(106_sphericalharmonics)
add_subdirectory(107_paramgeo2d)
add_subdirectory(108_paramgeo3d)
add_subdirectory(109_ContourLine2D)
add_subdirectory(110_ReadFileBVH)
add_subdirectory(111_ReadFileGLTF)
add_subdirectory(112_OffscreenRendering)


# meshing
add_subdirectory(201_dyntri2d_triangulation)
add_subdirectory(202_dyntri2d_remesh)
add_subdirectory(203_dyntri3d_edgecollapse)
add_subdirectory(204_dyntet_tetrahedralization)


# immplicit shape
add_subdirectory(300_AdaptiveDistanceField)
add_subdirectory(301_IsoSurfaceStuffing)


add_subdirectory(400_LaplaceMeshDeform)
add_subdirectory(401_ExponentialMap)
add_subdirectory(402_MapOnMesh)


# cad
add_subdirectory(500_2DCad)
add_subdirectory(501_2DCadEdit)
add_subdirectory(502_2DCadMeshEdit)
add_subdirectory(503_3DCad)

# fem
add_subdirectory(600_Fem2DHelmholtz)
add_subdirectory(601_FemPlateBendingMitc3)
add_subdirectory(602_Fem3DEigen)
add_subdirectory(603_Fem3DStiffwarp)
add_subdirectory(604_FemCloth_DynamicTriangle)

# spatial hash 
add_subdirectory(700_HashBvhSelfCollision)
add_subdirectory(701_HashLbvh3D)

# pbd
add_subdirectory(800_PbdGrid2D)
add_subdirectory(801_PbdCloth)

#bem
add_subdirectory(900_BemPotentialFlow3D)

