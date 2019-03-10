#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/funcs_gl.h"
#include "delfem2/mshtopoio_gl.h"

namespace py = pybind11;


void DrawField(const py::array_t<double>& pos,
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
  }
}

void init_field(py::module &m){
  m.def("draw_field",&DrawField);
}
