#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/mshtopoio_gl.h"

namespace py = pybind11;

CMeshElem HightMap(py::array_t<double> input0){
  py::buffer_info buf1 = input0.request();
  if (buf1.ndim != 2 ){
    throw std::runtime_error("Number of dimensions must be one");
  }
  const int nx = buf1.shape[0];
  const int ny = buf1.shape[1];
  CMeshElem msh;
  msh.elem_type = MESHELEM_QUAD;
  msh.ndim = 3;
  {
    double *ptr1 = (double *) buf1.ptr;
    std::vector<double> aXY;
    MeshTri2D_Grid(aXY, msh.aElem,
                   nx-1, ny-1);
    const int np = aXY.size()/2;
    msh.aPos.resize(np*3);
    for(int iy=0;iy<ny;++iy){
      for(int ix=0;ix<nx;++ix){
        int ip = iy*nx+ix;
        msh.aPos[ip*3+0] = aXY[ip*2+0];
        msh.aPos[ip*3+1] = aXY[ip*2+1];
        msh.aPos[ip*3+2] = ptr1[ip];
      }
    }
  }
  return msh;
}

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
  
  py::class_<CMeshMultiElem>(m,"MeshMultiElem")
  .def(py::init<>())
  .def("read_obj", &CMeshMultiElem::ReadObj)
  .def("minmax_xyz", &CMeshMultiElem::AABB3_MinMax)
  .def("draw",&CMeshMultiElem::Draw)
  .def("scaleXYZ",&CMeshMultiElem::ScaleXYZ)
  .def("translateXYZ",&CMeshMultiElem::TranslateXYZ);
  
  py::class_<CTriangulationOutput>(m, "TriangulationOutput")
  .def(py::init<>())
  .def_readonly("meshElem", &CTriangulationOutput::me);
  
  m.def("triangulation",&Triangulation,
        py::arg("aXY"),
        py::arg("edge_length")=0.03);
  
  m.def("hight_map", &HightMap);
  
}
