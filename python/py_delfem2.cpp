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
    Read_Ply(fname, aPos, aElem);
    elem_type = MESHELEM_TRI;
    ndim = 3;
  }
  void DrawFace_ElemWiseNorm(){
    if( elem_type == MESHELEM_TRI ){
      if( ndim == 3 ){ DrawMeshTri3D_FaceNorm(aPos, aElem); }
    }
    else if( elem_type == MESHELEM_QUAD ){
      if( ndim == 3 ){ DrawMeshQuad3D_FaceNorm(aPos, aElem); }
    }
  }
  void ScaleXYZ(double s){
    Scale(s,aPos);
  }
public:
  MESHELEM_TYPE elem_type;
  std::vector<int> aElem;
  /////
  int ndim;
  std::vector<double> aPos;
};

CMeshElem MeshQuad3D_VoxelGrid(const CVoxelGrid& vg){
  CMeshElem me;
  vg.GetQuad(me.aPos, me.aElem);
  me.elem_type = MESHELEM_QUAD;
  me.ndim = 3;
  return me;
}

namespace py = pybind11;

PYBIND11_MODULE(dfm2, m) {
  m.doc() = "pybind11 example plugin";
  //////
  py::class_<CMeshElem>(m, "MeshElem")
  .def(py::init<>())
  .def(py::init<const std::string&>())
  .def("read", &CMeshElem::Read)
  .def("drawFace_elemWiseNorm", &CMeshElem::DrawFace_ElemWiseNorm)
  .def("scaleXYZ",&CMeshElem::ScaleXYZ)
  .def_readonly("listElem", &CMeshElem::aElem)
  .def_readonly("listPos",  &CMeshElem::aPos)
  .def_readonly("elemType", &CMeshElem::elem_type);
  
  py::class_<CVoxelGrid>(m, "VoxelGrid")
  .def(py::init<>())
  .def("add",&CVoxelGrid::Add);
  
  py::enum_<MESHELEM_TYPE>(m, "FemElemType")
  .value("Tri",     MESHELEM_TYPE::MESHELEM_TRI)
  .value("Quad",    MESHELEM_TYPE::MESHELEM_QUAD)
  .value("Tet",     MESHELEM_TYPE::MESHELEM_TET)
  .value("Pyramid", MESHELEM_TYPE::MESHELEM_PYRAMID)
  .value("Wedge",   MESHELEM_TYPE::MESHELEM_WEDGE)
  .value("Hex",     MESHELEM_TYPE::MESHELEM_HEX)
  .export_values();
  
  m.def("setSomeLighting", &setSomeLighting);
  m.def("meshQuad3d_voxelGrid", &MeshQuad3D_VoxelGrid);
}

