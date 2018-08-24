#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/voxel.h"

#include "delfem2/funcs_gl.h"
//#include "delfem2/funcs_glut.h"

class CMeshElem{
public:
  CMeshElem(){}
  CMeshElem(const std::string& fpath){
    this->Read(fpath);
  }
  void Read(const std::string& fname){
    Read_Ply(fname, aPos, aTri);
    elem_type = MESHELEM_TRI;
    ndim = 3;
  }
  void Draw(){
    if( elem_type == MESHELEM_TRI ){
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
  py::class_<CMeshElem>(m, "MeshElem")
  .def(py::init<>())
  .def(py::init<const std::string&>())
  .def("read", &CMeshElem::Read)
  .def("draw", &CMeshElem::Draw)
  .def("scale_xyz",&CMeshElem::ScaleXYZ)
  .def_readonly("array_tri", &CMeshElem::aTri)
  .def_readonly("array_pos", &CMeshElem::aPos)
  .def_readonly("elem_type", &CMeshElem::elem_type);
  
  py::class_<CVoxelGrid>(m, "VoxelGrid")
  .def(py::init<>());
  
  py::enum_<MESHELEM_TYPE>(m, "FemElemType")
  .value("Tri",     MESHELEM_TYPE::MESHELEM_TRI)
  .value("Quad",    MESHELEM_TYPE::MESHELEM_QUAD)
  .value("Tet",     MESHELEM_TYPE::MESHELEM_TET)
  .value("Pyramid", MESHELEM_TYPE::MESHELEM_PYRAMID)
  .value("Wedge",   MESHELEM_TYPE::MESHELEM_WEDGE)
  .value("Hex",     MESHELEM_TYPE::MESHELEM_HEX)
  .export_values();
  
  m.def("set_some_lighting", &setSomeLighting);
}

