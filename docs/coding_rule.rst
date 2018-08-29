Coding Rule
===========


Overall
-------

- remove dependency between the codes as much as possible
- low level code must compile C++98
- use double space for a tab


Element Indexing Rule
---------------------

Based on the elemnt index rule in VTK(`VTK Manual <https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf>`_)

.. image:: element_index.png


Coordinate for the Depth Computation
------------------------------------

.. image:: depth.png


Naming Convention (cpp)
------------------------
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