# DelFEM2 C++ Examples using Legacy OpenGL and Glad

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



## Simple Demo

### [00_Noise3](00_Noise3)
<img src="00_Noise3/thumbnail.png" width=200px>

### [01_GlBuffer](01_GlBuffer)
<img src="01_GlBuffer/thumbnail.png" width=200>



## VBO

### [11_OffscreenRendering](11_OffscreenRendering)
<img src="11_OffscreenRendering/thumbnail.png" width=200px>

### [12_ProjectionBox](12_ProjectionBox)
<img src="12_ProjectionBox/thumbnail.png" width=200px>

### [13_CollisionLineHightfield](13_CollisionLineHightfield)
<img src="13_CollisionLineHightfield/thumbnail.png" width=200px>

### [14_Voxelize](14_Voxelize)
<img src="14_Voxelize/thumbnail.png" width=200>

### [15_VoxelMorph](15_VoxelMorph)
<img src="15_VoxelMorph/thumbnail.png" width=200>

### [16_VoxelGeodesic](16_VoxelGeodesic)
<img src="16_VoxelGeodesic/thumbnail.png" width=200>

### [17_RigVoxelGeodesic](17_RigVoxelGeodesic)
<img src="17_RigVoxelGeodesic/thumbnail.png" width=200>

### [18_DefLaplacianFitProj](18_DefLaplacianFitProj)
<img src="18_DefLaplacianFitProj/thumbnail.png" width=200>



## Shader

### [20_Shader](20_Shader)
<img src="20_Shader/thumbnail.png" width=200>

### [21_ShaderTexLaplace](21_ShaderTexLaplace)
<img src="21_ShaderTexLaplace/thumbnail.png" width=200>

### [22_ShaderContour](22_ShaderContour)
<img src="22_ShaderContour/thumbnail.png" width=200>

