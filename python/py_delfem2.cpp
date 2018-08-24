#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"

#include "delfem2/funcs_gl.h"
//#include "delfem2/funcs_glut.h"

class CMeshElem{
public:
  void Read(const std::string& fname){
    Read_Ply(fname, aPos, aTri);
    elem_type = FEMELEM_TRI;
    ndim = 3;
  }
  void Draw(){
    if( elem_type == FEMELEM_TRI ){
      if( ndim == 3 ){ DrawMeshTri3D_FaceNorm(aPos, aTri); }
    }
  }
  void ScaleXYZ(double s){
    Scale(s,aPos);
  }
public:
  MESHELEM_TYPE elem_type;
  std::vector<int> aTri;
  /////
  int ndim;
  std::vector<double> aPos;
};

namespace py = pybind11;

PYBIND11_MODULE(dfm2, m) {
  m.doc() = "pybind11 example plugin";
  //////
  py::class_<CMeshElem>(m, "MeshTri")
  .def(py::init<>())
  .def("read", &CMeshElem::Read)
  .def("draw", &CMeshElem::Draw)
  .def("scale_xyz",&CMeshElem::ScaleXYZ)
  .def_readonly("array_tri", &CMeshElem::aTri)
  .def_readonly("array_xyz", &CMeshElem::aPos)
  .def_readonly("elem_type", &CMeshElem::elem_type);
  
  py::enum_<MESHELEM_TYPE>(m, "FemElemType")
  .value("Tri",     MESHELEM_TYPE::FEMELEM_TRI)
  .value("Quad",    MESHELEM_TYPE::FEMELEM_QUAD)
  .value("Tet",     MESHELEM_TYPE::FEMELEM_TET)
  .value("Pyramid", MESHELEM_TYPE::FEMELEM_PYRAMID)
  .value("Wedge",   MESHELEM_TYPE::FEMELEM_WEDGE)
  .value("Hex",     MESHELEM_TYPE::FEMELEM_HEX)
  .export_values();
  
  m.def("set_some_lighting", &setSomeLighting);
}

