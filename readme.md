[![wercker status](https://app.wercker.com/status/03b6d924ec82270e22a04c3584fbf4de/s/master "wercker status")](https://app.wercker.com/project/byKey/03b6d924ec82270e22a04c3584fbf4de)  [![travis_status](https://travis-ci.org/nobuyuki83/delfem2.svg?branch=master)](https://travis-ci.org/nobuyuki83/delfem2)


# DelFEM2

DelFEM2 is a end-to-end framework for geometry processing and FEM simulation covering wide range of components including shape editing, meshing, FEM simulation, linear solver, and visualization. DelFEM2 is aiming an interactive digital engineering and authoring tool. 

Aside from the C++ implementation, python wrapper called PyDelFEM2 is provided. PyDelFEM2 can run various types of FEM simulation just a 10-20 lines of codes. Here is the example of solving the Poisson's equation in a square domain.

```
import PyDelFEM2 as dfm2
import PyDelFEM2.gl._glfw

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
 

The implementation is based on the [DelFEM](https://github.com/nobuyuki83/DelFEM) library

Please find out more detail in this [project document](https://nobuyuki83.github.io/delfem2/)

# Manual Tutorial

Please look at the python exaxmples under  ```examples_cpp```, ```examples_py```, ```examples_pyqt```, and ```examples_jupyter```.

Reading the test code under ```tests\``` also help understanding the code behaviour


# Dependency

PyDelFEM runs on Python3. Python2 is not supported.

PyDelFEM2 depends on following python packages:
- numpy
- glfw
- PyOpenGL  
- PySide2

These dependency is written in ```REQUIRED_PACKAGES``` in the setup.py, so they are automatically installed when installing the PyDelFEM2 pakage using the ```setup.py``` or ```pip3```.


# Install

## from PyPl

PyDelFEM2 can be installed simply with 

```
pip3 install PyDelFEM2
```

For Ubuntu, if you don't have pip installed, get it with ```sudo apt-get install python3-pip```

The installation fails if OpenGL packages are missing. For Ubuntu, install them ```sudo apt-get install freeglut3-dev libglfw3-dev libglew-dev```



## from GitHub

Installation from the GitHub repository can be done with the command:
```
pip3 install git+https://github.com/nobuyuki83/delfem2
```


## building from the source code

```
git clone https://github.com/nobuyuki83/delfem2.git
git submodle update --init --recursive
python3 setup.py install
```

For Ubuntu if you don't have git install it with ```sudo apt-get install git```



