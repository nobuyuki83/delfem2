#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/funcs_gl.h"
#include "delfem2/mshtopoio_gl.h"

namespace py = pybind11;


void DrawField_ColorMap
(const py::array_t<double>& pos,
 const py::array_t<int>& elm,
 const py::array_t<double>& val,
 const CColorMap& color_map)
{
//  DrawMeshTri2D_ScalarP1(me.aPos,me.aElem,a.data(),1,0,colorMap);
  const int np = pos.shape()[0];
  const int ndim = pos.shape()[1];
  const int nelm = elm.shape()[0];
  const int nnoelm = elm.shape()[1];
  const int ndimv = val.shape()[1];
  assert( val.shape()[0] == np);
  if( ndim == 3 ){
    DrawMeshTri3D_ScalarP1(pos.data(), np,
                           elm.data(), nelm,
                           val.data(),
                           color_map.aColor);
  }
  else if( ndim == 2 ){
    DrawMeshTri2D_ScalarP1(pos.data(), np,
                           elm.data(), nelm,
                           val.data(), 1,0,
                           color_map.aColor);
  }
}

void DrawField_Disp
(const py::array_t<double>& pos,
 const py::array_t<int>& elm,
 const py::array_t<double>& disp)
{
  //  DrawMeshTri2D_ScalarP1(me.aPos,me.aElem,a.data(),1,0,colorMap);
  const int np = pos.shape()[0];
  const int ndim = pos.shape()[1];
  const int nelm = elm.shape()[0];
  const int nnoelm = elm.shape()[1];
  const int ndimv = disp.shape()[1];
  assert( disp.shape()[0] == np );
  assert( disp.shape()[1] == ndim );
  if( ndim == 3 ){
  }
  else if( ndim == 2 ){
    DrawMeshTri2D_FaceDisp2D(pos.data(), np,
                             elm.data(), nelm,
                             disp.data());
  }
}

void init_field(py::module &m){
  m.def("drawField_colorMap",&DrawField_ColorMap);
  m.def("drawField_disp",&DrawField_Disp);
}
