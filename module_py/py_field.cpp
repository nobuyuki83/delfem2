#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/funcs_gl.h"
#include "delfem2/mshtopoio_gl.h"

namespace py = pybind11;


void DrawField(const CMeshElem& me,
               const py::array_t<double>& a,
               const CColorMap& color_map)
{
//  DrawMeshTri2D_ScalarP1(me.aPos,me.aElem,a.data(),1,0,colorMap);
    DrawMeshTri3D_ScalarP1(me.aPos,me.aElem,a.data(),
                           color_map.aColor);
}

void init_field(py::module &m){
  m.def("draw_field",&DrawField);
}
