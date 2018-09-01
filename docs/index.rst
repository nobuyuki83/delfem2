.. delfem2 documentation master file, created by
   sphinx-quickstart on Wed Aug 29 15:21:09 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DelFEM2 toolset for geometry processing and FEM
===============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   examples
   coding_rule   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Building Python Module
======================

:: 
  
  cd module_py
  mkdir buildMake
  cd buildMake
  cmake ..
  make

You can also use ``pip`` to　install the module_py using the ``setup.py``

::

  pip install . -e



Examples
========
The examples of the project built with DelFEM2 can be find in :doc:`examples`

Rule
====
Coding and naming rule can be found here: :doc:`coding_rule`


Refence
=======

.. autofunction:: dfm2.winDraw3d

.. autofunction:: dfm2.imgDraw3d

.. autofunction:: dfm2.setSomeLighting

.. autoclass::  dfm2.AABB3
   :undoc-members:

.. autoclass:: dfm2.VoxelGrid
   :undoc-members:

.. autofunction:: dfm2.meshQuad3d_voxelGrid