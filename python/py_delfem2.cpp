#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include<pybind11/numpy.h>

#include <vector>

#include "delfem2/voxel.h"
//#include "delfem2/funcs_glut.h"

#include "delfem2/mshtopoio_gl.h"

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

py::array_t<double> add_arrays
(py::array_t<double> input1,
 py::array_t<double> input2)
{
  auto buf1 = input1.request(), buf2 = input2.request();
  
  if (buf1.size != buf2.size)
    throw std::runtime_error("Input shapes must match");
  
  /*  allocate the buffer */
  py::array_t<double> result = py::array_t<double>(buf1.size);
  
  auto buf3 = result.request();
  
  double *ptr1 = (double *) buf1.ptr,
  *ptr2 = (double *) buf2.ptr,
  *ptr3 = (double *) buf3.ptr;
  int X = buf1.shape[0];
  int Y = buf1.shape[1];
  
  for (size_t idx = 0; idx < X; idx++)
    for (size_t idy = 0; idy < Y; idy++)
      ptr3[idx*Y + idy] = ptr1[idx*Y+ idy] + ptr2[idx*Y+ idy];
  
  // reshape array to match input shape
  result.resize({X,Y});
  
  return result;
}


void init_mshtopoio_gl(py::module &m);

PYBIND11_MODULE(_dfm2, m) {
  m.doc() = "pybind11 delfem2 binding";
  // mesh
  init_mshtopoio_gl(m);
  
  // voxel
  py::class_<CVoxelGrid>(m, "VoxelGrid")
  .def(py::init<>())
  .def("add",&CVoxelGrid::Add);
  m.def("meshQuad3d_voxelGrid", &MeshQuad3D_VoxelGrid);
  
  // gl misc
  m.def("setSomeLighting", &setSomeLighting);
}

