# DelFEM2

[![wercker status](https://app.wercker.com/status/03b6d924ec82270e22a04c3584fbf4de/s/master "wercker status")](https://app.wercker.com/project/byKey/03b6d924ec82270e22a04c3584fbf4de) [![travis_status](https://travis-ci.org/nobuyuki83/delfem2.svg?branch=master)](https://travis-ci.org/nobuyuki83/delfem2)


DelFEM2 is a toolset for Finite Element Analsys with Python/C++

- Seemless integration with Numpy.
- Integrated CAD / Mesher / Solver / Visualizer


## Features

- Built in mesh generation
	- 2D triangle mesh (Delaunay triangulation)
	- 3D tetrahedra mesh (Isosurface stuffing)
	- Dynamic remesh
- Built in sparse linar solver
	- Congugate gradient method
	- BiCGStab method
	- ILU(0) preconditioner
- Supports many equations on many types of mesh
	- Poission's equation (Tri P1, Tet P1)
	- Diffusion equation (Tri P1, Tet P1)
	- Linear solid equation (Tri P1, Tet P1)
	- Simulation of cloth (Tri P1)
	- Storkes equation (Tri P1, Tet P1)
	- Navier-Storkes equation with PSPG-SUPG stabilization (Tri P1)
- Visualization
	- color contour
	- mesh displacement
- IO
	- VTK file output





## Download

Please download the latest project by cloning the [GitHub repository](https://github.com/nobuyuki83/delfem2)

## Project layout

    docs          # This document
    src_cpp/
        cpp       # cpp files
        include   # include files for cpp
        external  # external dependencies
    module_py
        dfm2      # python module folder
    example_cpp   # cpp examples
    example_py    # python examples
    test          # Unit tests for C++/Python
    test_inputs   # Input files for test

## Use DelFEM2 from Python

To build the python module pleese see [building instruction](install).

For the examples see [python examples](example_py).


## Use DelFEM2 from C++

See the [cpp examples](example_cpp)


## Licence & Copyright

DelFEM2 is distributed under the [MIT licence](https://github.com/nobuyuki83/delfem2/blob/master/LICENSE). 

In a nut shell, it is free to use/distribute this software commercially for any purpose. But please always keep the copyright notice below and the MIT lisence statement.


	Copyright (C) Nobuyuki Umetani 2019


## Contact

DelFEM2 is currently developed and maintained by [Nobuyuki Umetani](http://www.nobuyuki-umetani.com/). If you have questions or comments please [contact](mailto:n.umetani@gmail.com).


