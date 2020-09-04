![](docs/imgs/social_preview.png)


<a href="http://doge.mit-license.org"><img src="http://img.shields.io/:license-mit-blue.svg"></a> 

| Linux (Travis-CI) | Windows (GitHub Action) |
|----|----|
| [![travis_status](https://travis-ci.org/nobuyuki83/delfem2.svg?branch=master)](https://travis-ci.org/nobuyuki83/delfem2) | ![](https://github.com/nobuyuki83/delfem2/workflows/CI_Windows/badge.svg) |



# DelFEM2

DelFEM2 is a end-to-end framework for geometry processing and FEM simulation covering wide range of components including shape editing, meshing, FEM simulation, linear solver, variational mesh deformer, and visualization. DelFEM2 is aiming to be an interactive digital engineering and authoring tool.


The implementation is based on my old open source project [DelFEM](https://github.com/nobuyuki83/DelFEM) library


***
# Manual &  Tutorial

There are currently no tutorial available for this library. To understand the code, please look at the exaxmples and tests  under following directoris.

+ C++ examples
  + [delfem2/examples_glfwold](examples_glfwold): examples with the legacy OpenGL
  + [delfem2/examples_glfwnew](examples_glfwnew):  examples with the modern OpenGL
  + [delfem2/examples_smpl](delfem2/examples_smpl): example using SMPL model
  + [delfem2/examples_cuda](examples_cuda): examples using cuda
+ C++ test:
  + [delfem2/test_cpp](test_cpp): tests using C++
  + [delfem2/test_cuda](test_cuda) : test using cuda



***
# Install

No installation is necessary to use DelFEM2. All the source code of DelFEM2 can be download using git using the following command:
```
git clone https://github.com/nobuyuki83/delfem2.git
```

DelFEM2 can be compiled either as a header-only library or as a static library. Nothing complicated is necessary if DelFEM2 is used as a header only library -- just by include header files and compile the code with option ```DFM2_HEADER_ONLY```. To use DelFEM2 as a static library, you may compiles  several dependent DelFEM2 source files and link them manually (this is not very complicated too).

Most of the source code of DelFEM2 does not have any external dependency. However, for some advancd functionality OpenGL or Unit Test, you many need to download dependent repositories and compile them manually. One can download all the dependent C++ repositories with

```
git submodle update --init --recursive
```

This command downloads all the external third party codes into the directory ```delfem2/3rd_party```. Currently DelFEM have binding to the following C++ open source projects.

- glfw
- glad
- cereal
- googletest
- cnpy
- imgui
- pybind11
- tinygltf
- stb_image.h

These projects are awesome and I would like to express huge  appreciation for contributers of these projects.






***

# Acknowledgement

For the testing, DelFEM2 used following data, images and 3D models:

- ```test_piputs/RiggedFigure.glb```from https://github.com/cx20/gltf-test/tree/master/sampleModels/RiggedFigure
- ```test_inputs/CesiumMan.glb``` from https://github.com/cx20/gltf-test/tree/master/sampleModels/CesiumMan 
- ```test_inputs/rollsRoyce.obj``` from https://poly.google.com/view/3DtJTlxgO_U
- ```test_inputs/walk.bvh``` from https://sites.google.com/a/cgspeed.com/cgspeed/motion-capture/cmu-bvh-conversion
- ```test_inputs/jump.bvh``` from https://sites.google.com/a/cgspeed.com/cgspeed/motion-capture/cmu-bvh-conversion
- ```test_input/uglysweater.jpg``` from https://www.flickr.com/photos/74812695@N00/6516495819
- ```test_input/bolt.STEP``` from https://grabcad.com/library/anchor-bolt-2
- ```test_inputs/epping_forest_01_1k.hdr``` from https://hdrihaven.com/hdri/?h=epping_forest_01


