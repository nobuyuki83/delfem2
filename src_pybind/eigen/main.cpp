#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/eigen_rigidbody.h"

namespace py = pybind11;

void init_rigidbody(py::module &m){
  py::class_<CRigidBodyAssembly_Static>(m,"RigidBodyAssembly_Static")
  .def(py::init<>())
  .def(py::init<std::vector<CRigidBody>,std::vector<CJoint>>())
  .def("solve", &CRigidBodyAssembly_Static::Solve)
  .def("draw",&CRigidBodyAssembly_Static::Draw)
  .def("minmax_xyz",&CRigidBodyAssembly_Static::MinMaxXYZ);
//  .def("init_gl",&CTexture::LoadTex)
//  .def("minmax_xyz",&CTexture::MinMaxXYZ);
  
  py::class_<CRigidBody>(m,"RigidBody")
  .def(py::init<double, std::vector<double>>())
  .def("add_contact_point",&CRigidBody::addCP);
  
  py::class_<CJoint>(m,"Joint")
  .def(py::init<int,int, std::vector<double>>());
}


PYBIND11_MODULE(c_eigen, m) {
  m.doc() = "pybind11 delfem2 binding";
  ///////////////////////////////////
  init_rigidbody(m);
}
