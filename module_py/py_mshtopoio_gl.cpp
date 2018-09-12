#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "delfem2/mshtopoio_gl.h"

namespace py = pybind11;

void init_mshtopoio_gl(py::module &m){
  py::enum_<MESHELEM_TYPE>(m, "MeshElemType")
  .value("Tri",     MESHELEM_TYPE::MESHELEM_TRI)
  .value("Quad",    MESHELEM_TYPE::MESHELEM_QUAD)
  .value("Tet",     MESHELEM_TYPE::MESHELEM_TET)
  .value("Pyramid", MESHELEM_TYPE::MESHELEM_PYRAMID)
  .value("Wedge",   MESHELEM_TYPE::MESHELEM_WEDGE)
  .value("Hex",     MESHELEM_TYPE::MESHELEM_HEX)
  .export_values();
  
  py::class_<CMeshElem>(m, "MeshElem")
  .def(py::init<>())
  .def(py::init<const std::string&>(),"open the file in the path")
  .def("draw", &CMeshElem::Draw) //
  .def("minmax_xyz", &CMeshElem::AABB3_MinMax)
  .def("read", &CMeshElem::Read)
  .def("drawFace_elemWiseNorm", &CMeshElem::DrawFace_ElemWiseNorm)
  .def("drawEdge", &CMeshElem::DrawEdge)
  .def("scaleXYZ", &CMeshElem::ScaleXYZ)
  .def("subdiv",   &CMeshElem::Subdiv)
  .def_readonly("listElem", &CMeshElem::aElem)
  .def_readonly("listPos",  &CMeshElem::aPos)
  .def_readonly("elemType", &CMeshElem::elem_type)
  .def_readonly("nDim",     &CMeshElem::ndim)
  .def_readwrite("color_face",  &CMeshElem::color_face);
  
  py::class_<CTriangulationOutput>(m, "TriangulationOutput")
  .def(py::init<>())
  .def_readonly("meshElem", &CTriangulationOutput::me);
  
  m.def("triangulation",&Triangulation,
        py::arg("aXY"),
        py::arg("edge_length")=0.03);
  
}
