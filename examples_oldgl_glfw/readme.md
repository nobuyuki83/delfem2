# DelFEM2 C++ Examples using Legacy OpenGL

These demos use OpenGL version 2.1 and GLSL shaer version 1.2 which are depricated in many environment. But still it is convenient to use legacy functions such as glBegin(), glEnd(). We will eventually consider porting these demo into newer OpenGL >= 3.3 in the [examples_newgl_glfw](../examples_newgl_glfw) folder.



## How To Build

These demos depend on the GLFW library to open an OpenGL window. If you haven't install `glfw` in your computer, please read following Document to set up GLFW.

- [How to Set Up GLFW Library](../docs/setup_glfw.md)
- [GLFWライブラリの設定方法](../docs/setup_glfw_jp.md)

With `glfw` installed, you can build the demos simply by 

```bash
mkdir build && cd build
cmake ..
cmake --build .
```



## Simple Example without DelFEM2

### [000_OpenWin](000_OpenWindowWithGLFW)
<img src="000_OpenWindowWithGLFW/thumbnail.png" width=200px>



## Basic Demo

### [100_Nav3D](100_ViewNavigationIn3D)
<img src="100_ViewNavigationIn3D/thumbnail.png" width=200px>

### [101_SubdivCatmull](101_SubdivCatmull)
<img src="101_SubdivCatmull/thumbnail.png" width=200px>

### [102_ConvexHull3](102_ConvexHull3)
<img src="102_ConvexHull3/thumbnail.png" width=200px>

### [103_Noise2d](103_Noise2)
<img src="103_Noise2/thumbnail.png" width=200px>

### [104_CageDef3Mvc](104_CageDef3Mvc)
<img src="104_CageDef3Mvc/thumbnail.png" width=200px>

### [105_MeshSlice](105_MeshSlices)
<img src="105_MeshSlices/thumbnail.png" width=200px>

### [106_SphericalHarmonics](106_SphericalHarmonics)
<img src="106_SphericalHarmonics/thumbnail.png" width=200px>

### [107_ParamGeo2d](107_ParamGeo2)
<img src="107_ParamGeo2/thumbnail.png" width=200px>

### [108_ParamGeo3d](108_ParamGeo3)
<img src="108_ParamGeo3/thumbnail.png" width=200px>

### [109_ContourLine2D](109_ContourLine2D)
<img src="109_ContourLine2D/thumbnail.png" width=200px>

### [110_ClusteringMesh](110_ClusteringMesh)
<img src="110_ClusteringMesh/thumbnail.png" width=200px>

### [111_ExponentialMapElem](111_ExponentialMapElem)
<img src="111_ExponentialMapElem/thumbnail.png" width=200px>

### [112_SamplingMeshTri3](112_SamplingMeshTri3)
<img src="112_SamplingMeshTri3/thumbnail.png" width=200px>

### [113_KnitGeometry](113_KnitGeometry)
<img src="113_KnitGeometry/thumbnail.png" width=200px>

### [114_ConvexHull2](114_ConvexHull2)
<img src="114_ConvexHull2/thumbnail.png" width=200px>

### [115_ExponentialMapPoint](115_ExponentialMapPoint)
<img src="115_ExponentialMapPoint/thumbnail.png" width=200px>

### [116_GridOnMesh](116_GridOnMesh)
<img src="116_GridOnMesh/thumbnail.png" width=200px>

### [117_QuadSubdivOnMesh](117_QuadSubdivOnMesh)
<img src="117_QuadSubdivOnMesh/thumbnail.png" width=200>

### [118_FastMarchingMethod](118_FastMarchingMethod)
<img src="118_FastMarchingMethod/thumbnail.png" width=200>

### [119_Glyph](119_Glyph)
<img src="119_Glyph/thumbnail.png" width=200>

### [120_Adf3](120_AdaptiveDistanceField)
<img src="120_AdaptiveDistanceField/thumbnail.png" width=200>

### [121_IsoSurfaceStuffing](121_IsoSurfaceStuffing)
<img src="121_IsoSurfaceStuffing/thumbnail.png" width=200>

### [122_HashBvhSelfCollision](122_HashBvhSelfCollision)
<img src="122_HashBvhSelfCollision/thumbnail.png" width=200>

### [123_HashLbvh3D](123_HashLbvh3D)
<img src="123_HashLbvh3D/thumbnail.png" width=200>

### [124_4RotSymField](124_4RotSymField)
<img src="124_4RotSymField/thumbnail.png" width=200>

### [125_GizmoRot](125_GizmoRot)
<img src="125_GizmoRot/thumbnail.png" width=200>

### [126_GizmoTransl](126_GizmoTransl)
<img src="126_GizmoTransl/thumbnail.png" width=200>

### [127_GizmoAffine](127_GizmoAffine)
<img src="127_GizmoAffine/thumbnail.png" width=200>

### [128_GjkSat2](128_GjkSat2)
<img src="128_GjkSat2/thumbnail.png" width=200>

### [129_CubeGridEdit](129_CubeGridEdit)
<img src="129_CubeGridEdit/thumbnail.png" width=200 alt="129_CubeVoxelEdit">

### [134_Primitives](134_Primitives)
<img src="134_Primitives/thumbnail.png" width=200 alt="134_Primitives">

### [137_VoxelLineIntersection](137_VoxelLineIntersection)
<img src="137_VoxelLineIntersection/thumbnail.png" width=200 alt="137_VoxleLineIntersection">

### [138_RigReadFileBiovision](138_RigReadFileBiovision)
<img src="138_RigReadFileBiovision/thumbnail.png" width=200>

### [139_Pick](139_Pick)
<img src="139_Pick/thumbnail.png" width=200>

### [140_PickPropagate](140_PickPropagate)
<img src="140_PickPropagate/thumbnail.png" width=200>

### [141_RigBendCapsule](141_RigBendCapsule)
<img src="141_RigBendCapsule/thumbnail.png" width=200>

### [142_BinaryClusteringPoints3D](142_BinaryClusteringPoints3D)
<img src="142_BinaryClusteringPoints3D/thumbnail.png" width=200>

### [143_Image2MeshInterpolation](143_Image2MeshInterpolation)
<img src="143_Image2MeshInterpolation/thumbnail.png" width=200>

### [144_BinaryClusteringPoints2D](144_BinaryClusteringPoints2D)
<img src="144_BinaryClusteringPoints2D/thumbnail.png" width=200>

### [145_Lloyd](145_Lloyd)
<img src="145_Lloyd/thumbnail.png" width=200>

### [146_4RotSymFieldMultilevel](146_4RotSymFieldMultilevel)
<img src="146_4RotSymFieldMultilevel/thumbnail.png" width=200>

### [147_HashLbvh3RayTri](147_HashLbvh3RayTri)
<img src="147_HashLbvh3RayTri/thumbnail.png" width=200>

### [148_RayCasting](148_RayCasting)
<img src="148_RayCasting/thumbnail.png" width=200>


## Dynamic Triangle

### [201_DynTri2Triangulation](201_DynTri2dTriangulation)
<img src="201_DynTri2dTriangulation/thumbnail.png" width=200>

### [202_DynTri2Remesh](202_DynTri2dRemesh)
<img src="202_DynTri2dRemesh/thumbnail.png" width=200>

### [203_DynTri3_EdgeCollapse](203_DynTri3dEdgeCollapses)
<img src="203_DynTri3dEdgeCollapses/thumbnail.png" width=200>

### [204_DynTet_Tetrahedralization](204_DynTetTetrahedralization)
<img src="204_DynTetTetrahedralization/thumbnail.png" width=200>

### [205_DynTri3SimpilfyQuadratic](205_DynTri3SimpilfyQuadratic)
<img src="205_DynTri3SimpilfyQuadratic/thumbnail.png" width=200>



## CAD

### [500_Cad2d](500_Cad2d)
<img src="500_Cad2d/thumbnail.png" width=200>

### [501_Cad2dEdit](501_Cad2dEdit)
<img src="501_Cad2dEdit/thumbnail.png" width=200>

### [502_Cad2dMeshEdit](502_Cad2dMeshEdit)
<img src="502_Cad2dMeshEdit/thumbnail.png" width=200>

### [503_Cad3d](503_Cad3d)
<img src="503_Cad3d/thumbnail.png" width=200>

### [504_Cad2dSvg](504_Cad2dSvg)
<img src="504_Cad2dSvg/thumbnail.png" width=200>

### [505_Cad3ReadStep](505_Cad3ReadStep)
<img src="505_Cad3ReadStep/thumbnail.png" width=200>



## Simulation (FEM,FDM,PBD,BEM,SPH,Arap)

### [600_Fem2Helmholtz](600_Fem2DHelmholtz)
<img src="600_Fem2DHelmholtz/thumbnail.png" width=200>

### [601_FemPlateBendingMitc3](601_FemPlateBendingMitc3)
<img src="601_FemPlateBendingMitc3/thumbnail.png" width=200>

### [602_Fem3Eigen](602_Fem3DEigen)
<img src="602_Fem3DEigen/thumbnail.png" width=200>

### [603_Fem3StiffWarp](603_Fem3DStiffwarp)
<img src="603_Fem3DStiffwarp/thumbnail.png" width=200>

### [604_FemClothDTri](604_FemCloth_DynamicTriangle)
<img src="604_FemCloth_DynamicTriangle/thumbnail.png" width=200>

### [605_FemClothInternal](605_FemClothInternal)
<img src="605_FemClothInternal/thumbnail.png" width=200>

### [606_PbdDeform2D](606_PbdGrid2D)
<img src="606_PbdGrid2D/thumbnail.png" width=200>

### [607_PbdCloth](607_PbdCloth)
<img src="607_PbdCloth/thumbnail.png" width=200>

### [608_StableFluids2Staggered](608_StableFluids2Staggered)
<img src="608_StableFluids2Staggered/thumbnail.png" width=200>

### [609_DefLaplacianMesh](609_DefLaplacian)
<img src="609_DefLaplacian/thumbnail.png" width=200>

### [610_DefArapEdge](610_DefArapEdge)
<img src="610_DefArapEdge/thumbnail.png" width=200>

### [611_DefArap](611_DefArap)
<img src="611_DefArap/thumbnail.png" width=200>

### [612_BemPotentialFlow3D](612_BemPotentialFlow3D)
<img src="612_BemPotentialFlow3D/thumbnail.png" width=200>

### [613_FemRod3Darboux](613_FemRod3Darboux)
<img src="613_FemRod3Darboux/thumbnail.png" width=200>

### [614_FemMitc3Eigen](614_FemMitc3Eigen)
<img src="614_FemMitc3Eigen/thumbnail.png" width=200>

### [615_PbdClothCad](615_PbdClothCad)
<img src="615_PbdClothCad/thumbnail.png" width=200>

### [616_FemRod3StraightStatic](616_FemRod3StraightStatic)
<img src="616_FemRod3StraightStatic/thumbnail.png" width=200>

### [617_FemRod3HairDarbouxSelfcollision](617_FemRod3HairDarbouxSelfcollision)
<img src="617_FemRod3HairDarbouxSelfcollision/thumbnail.png" width=200>

### [618_Fem3Tet](618_Fem3Tet)
<img src="618_Fem3Tet/thumbnail.png" width=200>

### [619_Sph3](619_Sph3)
<img src="619_Sph3/thumbnail.png" width=200>

### [620_FemClothSelfCollision](620_FemClothSelfCollision)
<img src="620_FemClothSelfCollision/thumbnail.png" width=200>

### [621_Fem2](621_Fem2)
<img src="621_Fem2/thumbnail.png" width=200 alt="621_Fem2">

### [622_FemHairDarbouxStatic](622_FemHairDarbouxStatic)
<img src="622_FemHairDarbouxStatic/thumbnail.png" width=200 alt="622_FemRodHairStatic">

### [623_FemHairDarbouxDynamic](623_FemHairDarbouxDynamic)
<img src="623_FemHairDarbouxDynamic/thumbnail.png" width=200 alt="622_FemRodHairDynamic">

### [624_DefArapUi](624_DefArapUi)
<img src="624_DefArapUi/thumbnail.png" width=200>

### [625_Fem2MasterSlave](625_Fem2MasterSlave)
<img src="625_Fem2MasterSlave/thumbnail.png" width=200 alt="625_Fem2MasterSlave">

### [626_RgdRotation](626_RgdRotation)
<img src="626_RgdRotation/thumbnail.png" width=200>

### [627_Rgd2dCollision](627_Rgd2dCollision)
<img src="627_Rgd2dCollision/thumbnail.png" width=200>

### [628_FemHyper3](628_FemHyper3)
<img src="628_FemHyper3/thumbnail.png" width=200>

### [629_Rgd2Caterpillar](629_Rgd2Caterpillar)
<img src="629_Rgd2Caterpillar/thumbnail.png" width=200>

### [630_FemRod2Static](630_FemRod2Static)
<img src="630_FemRod2Static/thumbnail.png" width=200>

### [631_FemRod2Dynamic](631_FemRod2Dynamic)
<img src="631_FemRod2Dynamic/thumbnail.png" width=200>

### [632_PbdSpring2](632_PbdSpring2)
<img src="632_PbdSpring2/thumbnail.png" width=200>

### [633_Mpm2Snow](633_Mpm2Snow)
<img src="633_Mpm2Snow/thumbnail.png" width=200>

### [634_FemLinearSolid3](634_FemLinearSolid3)
<img src="634_FemLinearSolid3/thumbnail.png" width=200>

### [635_FemRod2Edit](635_FemRod2Edit)
<img src="635_FemRod2Edit/thumbnail.png" width=200>

### [636_InvertibleFem](636_InvertibleFem)
<img src="636_InvertibleFem/thumbnail.png" width=200>

### [637_StableFluids2CellCentered](637_StableFluids2CellCentered)
<img src="637_StableFluids2CellCentered/thumbnail.png" width=200>

### [639_AccousticFTDT2](639_AccousticFTDT2)
<img src="639_AccousticFTDT2/thumbnail.png" width=200>
