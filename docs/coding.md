# Coding Convention



## Overall

- remove dependency between the codes as much as possible
- put implementations in different code based on its dependency 
- use double space for a tab





## Element Indexing Rule

Based on the elemnt index rule in [VTK](https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf)

![element_index](imgs/element_index.png)


## Coordinate for the Depth Computation


![depth](imgs/depth.png)



## Naming Convention 



### file name

use the extension ".h" instead of ".hpp"

There are roughly three types of source code and naming is different

1. stand alone source code -> short lowercase name (e.g, vec2.h, mats.h)

2. source code that depends on few other codes -> short lowercase name and its dependency connected by "_". (e.g, objf_geo3.cpp)

3. source code that depeds on many other codes -> Uppercase name

   

### function & class name

The function name should be written in camel case notation that sarts with upper case letter (e.g., Verb_Object_Adverb)

- Geometric Operator
  - Nearest
  - Intersection
  - IsInside
  - Volume
  - Area

- Use Read <-> Write for the file io. Don't use Load <-> Save

- Naming Various Types of Points
    - aXY -> general 2D point
    - aXYZ-> general 3D point
    - CAD vertex -> Vertex
    - Mesh corner point -> Point
    - FEM points it can be inside an element (may be on edge or on face) -> Node  

How to call a mesh?
- MeshTri3D
- MeshQuad2D
- MeshHex3D
- MeshMix3D
- MeshElem3D
- Point3D
- MeshHex3DSurface



### which function goes to which file?

Matrix, Vector & Quaternion

- The order of priority between objects is "Vec2 < Mat2 < Vec3 < Quat < Mat3 < Mat4". Higher priority object contains the functions that use lower priority object.



