![](docs/imgs/social_preview.png)


<a href="http://doge.mit-license.org"><img src="http://img.shields.io/:license-mit-blue.svg"></a> 

| Linux | Linux | Windows |
|----|----|----|
| [![wercker status](https://app.wercker.com/status/03b6d924ec82270e22a04c3584fbf4de/s/master "wercker status")](https://app.wercker.com/project/byKey/03b6d924ec82270e22a04c3584fbf4de) | [![travis_status](https://travis-ci.org/nobuyuki83/delfem2.svg?branch=master)](https://travis-ci.org/nobuyuki83/delfem2) | ![](https://github.com/nobuyuki83/delfem2/workflows/CI/badge.svg) |



# DelFEM2

DelFEM2 is a end-to-end framework for geometry processing and FEM simulation covering wide range of components including shape editing, meshing, FEM simulation, linear solver, variational mesh deformer, and visualization. DelFEM2 is aiming to be an interactive digital engineering and authoring tool.

Aside from the C++ implementation, python wrapper called PyDelFEM2 is provided. PyDelFEM2 can run various types of FEM simulation just a 10-20 lines of codes. Here is the example of solving the Poisson's equation in a square domain.

```
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

cad = dfm2.Cad2D()
cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,0, +0,+0, 0,+1, -1,+1.0])
mesh,map_cad2mesh = cad.mesh(0.05)
fem = dfm2.FEM_Poisson(source=1.0)
fem.updated_topology(mesh)
npIdP = cad.points_edge([0,1,2,3], mesh.np_pos)
fem.ls.bc[npIdP] = 1
fem.solve()
field = dfm2.VisFEM_ColorContour(fem,"value")
dfm2.gl._glfw.winDraw3d([field])
 ```
The result of this code woud be the following window

![Poisson](docs/imgs/poisson.png)


The implementation is based on my old open source project [DelFEM](https://github.com/nobuyuki83/DelFEM) library

Please find out more detail in this [project document](https://nobuyuki83.github.io/delfem2/)


***
# Manual &  Tutorial

There are currently no tutorial available for this library. To understand the code, please look at the exaxmples and tests  under following directoris.

+ C++ examples
  + [delfem2/examples_glfwold](examples_glfwold): examples with the legacy OpenGL
  + [delfem2/examples_glfwnew](examples_glfwnew):  examples with the modern OpenGL
  + [delfem2/examples_cuda](examples_cuda): examples using cuda
+ C++ test:
  + [delfem2/test_cpp](test_cpp): tests using C++
  + [delfem2/test_cuda](test_cuda) : test using cuda
+ Python examples
  + [delfem2/examples_py](examples_py) : examples using python
  + [delfem2/examples_pyqt](examples_pyqt) examples using PyQt gui library
  + [delfem2/examples_jupyter](examples_jupyter) : examples using Jupyter
  + [delfem2/examples_blender](examples_blender) : examples with Blender python scripting
+  Python tests
  + [delfem2/test_py](test_py) : test using python




***
# Install

## C++

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



## Python

The most recommended way to install PyDelFEM2, which is the python binding of DelFEM2, is to build it from the source code. The following command down load the source code and its C++ dependencies and build python modules and download its python dependencies.

```
git clone https://github.com/nobuyuki83/delfem2.git
git submodle update --init --recursive
pip3 install -e .
```

Here are some trouble shooting tips: 
- For Ubuntu if you don't have git install it with ```sudo apt-get install git```

- For Ubuntu, if you don't have pip installed, get it with ```sudo apt-get install python3-pip```

- The installation fails if OpenGL packages are missing. For Ubuntu, install them ```sudo apt-get install freeglut3-dev libglfw3-dev libglew-dev```




Alternatively, PyDelFEM2 can be also installed from the GitHub repository can be done with the command:
```
pip3 install git+https://github.com/nobuyuki83/delfem2
```

PyDelFEM2 can be installed from PyPL simply with the following command, however this is not recommended as the version in PyPL is not updated frequently.

```
pip3 install PyDelFEM2
```


PyDelFEM runs on Python3. Python2 is not supported. PyDelFEM2 depends on following awesome python packages:
- numpy
- glfw
- PyOpenGL
- PySide2

These dependencies are written in ```REQUIRED_PACKAGES``` in the setup.py, so they are automatically installed when installing the PyDelFEM2 pakage using the ```setup.py``` or ```pip3```.


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


