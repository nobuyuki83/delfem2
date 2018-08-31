#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include<pybind11/numpy.h>

#include <vector>

#include "delfem2/voxel.h"
//#include "delfem2/funcs_glut.h"

#include "delfem2/mshtopoio_gl.h"
#include "delfem2/bv.h"

// TODO:Make a wrapper class of the VoxelGrid?
CMeshElem MeshQuad3D_VoxelGrid(const CVoxelGrid& vg){
  CMeshElem me;
  vg.GetQuad(me.aPos, me.aElem);
  me.elem_type = MESHELEM_QUAD;
  me.ndim = 3;
  return me;
}

//////////////////////////////////////////////////////////////////////////////////////////

namespace py = pybind11;

void init_mshtopoio_gl(py::module &m);

PYBIND11_MODULE(dfm2, m) {
  m.doc() = "pybind11 delfem2 binding";
  // aabb
  py::class_<CBV3D_AABB>(m,"AABB3",
                         "three dimensional axis aligned bounding box class")
  .def(py::init<>())
  .def("minMaxLocXY",
       &CBV3D_AABB::MinMaxLocXY,
       "get range of the x and y local coordiante");
  
  // mesh
  init_mshtopoio_gl(m);
  
  // voxel
  py::class_<CVoxelGrid>(m, "VoxelGrid",
                         "voxel grid class")
  .def(py::init<>())
  .def("add",&CVoxelGrid::Add,"add voxel at the integer coordinate");
  
  m.def("meshQuad3d_voxelGrid",
        &MeshQuad3D_VoxelGrid,
        "get quad mesh from voxel grid");
  
  // gl misc
  m.def("setSomeLighting",
        &setSomeLighting,
        "set some lighting that looks good for me");
}

