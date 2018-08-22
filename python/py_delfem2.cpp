#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "delfem2/mshio.h"
#include "delfem2/msh.h"

#include "delfem2/funcs_gl.h"

class CMeshTri{
public:
  void Read(const std::string& fname){
    Load_Ply(fname, aPos, aTri);
  }
  void Draw(){
    DrawMeshTri3D_FaceNorm(aPos, aTri);
  }
  void ScaleXYZ(double s){
    Scale(s,aPos);
  }
public:
  std::vector<int> aTri;
  std::vector<double> aPos;
};

namespace py = pybind11;

PYBIND11_MODULE(dfm2, m) {
  m.doc() = "pybind11 example plugin";
  //////
  py::class_<CMeshTri>(m, "MeshTri")
  .def(py::init<>())
  .def("read", &CMeshTri::Read)
  .def("draw", &CMeshTri::Draw)
  .def("scale_xyz",&CMeshTri::ScaleXYZ)
  .def_readonly("array_tri", &CMeshTri::aTri)
  .def_readonly("array_xyz", &CMeshTri::aPos);
  
  m.def("set_some_lighitng", &setSomeLighting);
}

