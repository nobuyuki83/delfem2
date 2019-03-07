#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/rigidbody.h"

namespace py = pybind11;

void init_rigidbody(py::module &m){
  py::class_<CRigidBodyAssembly_Static>(m,"RigidBodyAssembly_Static")
  .def(py::init<>())
  .def("draw",&CRigidBodyAssembly_Static::Draw);
//  .def("init_gl",&CTexture::LoadTex)
//  .def("minmax_xyz",&CTexture::MinMaxXYZ);
}
