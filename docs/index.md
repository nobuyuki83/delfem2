# DelFEM2

[![wercker status](https://app.wercker.com/status/03b6d924ec82270e22a04c3584fbf4de/s/master "wercker status")](https://app.wercker.com/project/byKey/03b6d924ec82270e22a04c3584fbf4de) [![travis_status](https://travis-ci.org/nobuyuki83/delfem2.svg?branch=master)](https://travis-ci.org/nobuyuki83/delfem2)


DelFEM2 is a toolset for Finite Element Analsys with Python/C++

- Seemless integration with Numpy.
- Integrated CAD / Mesher / Solver / Visualizer


## Features

- Mesh Generation
	- 2D Triangle Mesh (Delaunay Triangulation)
	- 3D Triangle Mesh (IsoSurface Stuffing)
- Poission's Equation
	- 2D Triangle Mesh
	- 3D Tetrahedra Mesh
- Diffusion Equation
	- 2D Triangle Mesh
	- 3D Tetrahedra Mesh
- Linear Solid Equation
	- 2D Trianlge Mesh
	- 3D Tetrahedra Mesh
- Simulation of Cloth
	- 3D Triangle Mesh
- Storkes Equation
	- 2D Trianlge Mesh
	- 3D Tetrahedra Mesh
- Navier Storkes Equation
	- 2D Trianlge Mesh


## Download

Please download the latest project by cloning the [GitHub repository](https://github.com/nobuyuki83/delfem2)

## Project layout

    docs           # This document
    src_cpp/
        cpp        # cpp files
        include    # include files for cpp
        external   # external dependencies
    module_py
        dfm2       # python module folder
    example_cpp    # Cpp examples
    example_py     # Python examples

## Python

To build the python module pleese see [building instruction](install).

For the examples see [python examples](example_py).


## C++

See the [cpp examples](example_cpp)


## Licence

DelFEM2 is distributed under the [MIT licence](https://github.com/nobuyuki83/delfem2/blob/master/LICENSE)



## Contact

DelFEM2 is currently developed and maintained by [Nobuyuki Umetani](http://www.nobuyuki-umetani.com/). If you have questions or comments please [contact](mailto:n.umetani@gmail.com).


## Copyright

Nobuyuki Umetani -- All rights reserved.



