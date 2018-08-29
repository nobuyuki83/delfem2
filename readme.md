[![wercker status](https://app.wercker.com/status/03b6d924ec82270e22a04c3584fbf4de/s/master "wercker status")](https://app.wercker.com/project/byKey/03b6d924ec82270e22a04c3584fbf4de)     [![Documentation Status](https://readthedocs.org/projects/delfem2/badge/?version=latest)](https://delfem2.readthedocs.io/en/latest/?badge=latest)


# DelFEM2

A handy toolset for coding geometry processing and fem simulation

The implementation is based on the [DelFEM](https://github.com/nobuyuki83/DelFEM) library


## How to Build

### Build Examples
```
cd examples
mkdir buildMake
cd buildMake
cmake ..
make
```

### Build Python Binding
```
cd python
mkdir buildMake
cd buildMake
cmake ..
make
```


## Examples

- working directory is in the folder `examples/`
- binary files are put in the folder `examples/bin/`
- input files are put in the folder `test_inputs/`
- The scripts to run the executable is put in `examples/script/`


| Name | Screen Shot |
| ------------- | ------------- |
| triangulation | ![triangulation](docs/screenshot_triangulation.png) |
| transform_handler  | ![handler](docs/screenshot_handler.png) |
| cloth_internal | ![cloth_internal](docs/screenshot_clothinternal.png) |
| subdiv | ![subdiv](docs/screenshot_subdiv.png) |
| read_bvh | ![read_bvh](docs/screenshot_readbvh.png) |
| exponential_map | ![expmap](docs/screenshot_expmap.png) |
| selfcollision_bvh | ![selfcollisionbvh](docs/screenshot_selfcollisionbvh.png) |
| edge_collapse | ![edgecollapse](docs/screenshot_edgecollapse.png) |


## Coding

### Philosophy
- remove dependency between the codes as much as possible
- low level code must compile C++98


### Rule
- use double space for a tab


### Element Indexing Rule

Based on the elemnt index rule in VTK( https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf )

![element index](./docs/element_index.png)


### Coordinate for the Depth Computation

![depth_coord](./docs/depth.png)


### Naming Convention (cpp)
* use the extension ".h" instead of ".hpp"
* filename should be in lower case. The underscore represent depndency. For example "aa_bb.cpp" means this is a implementation of class "aaa" and it depends on a class "bbb"
* The function name should be written in camel case notation that sarts with upper case letter (e.g., Verb_Object_Adverb)
* Geometric Operator
  * Nearest
  * Intersection
  * IsInside
  * Volume
  * Area
* Use Read <-> Write for the file io. Don't use Load <-> Save
* Naming Point
  * aXY -> general 2D point
  * aXYZ-> general 3D point
  * CAD vertex -> Vertex
  * Mesh corner point -> Point
  * FEM points it can be inside an element (may be on edge or on face) -> Node
* How to call a mesh?
  * MeshTri3D
  * MeshQuad2D
  * MeshHex3D
  * MeshMix3D
  * MeshElem3D
  * Point3D
  * MeshHex3DSurface