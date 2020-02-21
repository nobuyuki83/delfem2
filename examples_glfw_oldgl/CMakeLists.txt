cmake_minimum_required(VERSION 2.8)

project(examples_glfw_oldgl)


# simple examples without delfem2 
add_subdirectory(000_OpenWindowWithGLFW)


# basic funcs
add_subdirectory(100_ViewNavigationIn3D)
add_subdirectory(101_SubdivCatmull)
add_subdirectory(102_ConvexHull)
add_subdirectory(103_Noise2)
add_subdirectory(104_Noise3)
add_subdirectory(105_MeshSlices)
add_subdirectory(106_SphericalHarmonics)
add_subdirectory(107_ParamGeo2)
add_subdirectory(108_ParamGeo3)
add_subdirectory(109_ContourLine2D)
add_subdirectory(110_ReadFileBVH)
add_subdirectory(111_ReadFileGLTF)
add_subdirectory(112_OffscreenRendering)
add_subdirectory(113_ProjectionBox)


# meshing
add_subdirectory(201_DynTri2dTriangulation)
add_subdirectory(202_DynTri2dRemesh)
add_subdirectory(203_DynTri3dEdgeCollapses)
add_subdirectory(204_DynTetTetrahedralization)


# immplicit shape
add_subdirectory(300_AdaptiveDistanceField)
add_subdirectory(301_IsoSurfaceStuffing)
add_subdirectory(302_HashBvhSelfCollision)
add_subdirectory(303_HashLbvh3D)

add_subdirectory(400_ExponentialMap)
add_subdirectory(401_GridOnMesh)
add_subdirectory(402_QuadSubdivOnMesh)
add_subdirectory(403_LaplacianMeshDeformation)
add_subdirectory(404_AsRigidAsPossibleDeformation)


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
add_subdirectory(605_FemClothInternal)
add_subdirectory(606_PbdGrid2D)
add_subdirectory(607_PbdCloth)
add_subdirectory(608_FdmFluid2D)

#bem
add_subdirectory(900_BemPotentialFlow3D)

