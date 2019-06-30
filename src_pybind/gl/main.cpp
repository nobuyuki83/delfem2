#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <map>
#include <deque>

#include "delfem2/funcs_glew.h" // have to be included in the beginning
#include "delfem2/funcs_gl.h"
#include "delfem2/color_gl.h"
//#include "delfem2/mshtopoio_gl.h"

//#include "delfem2/voxel.h"
//#include "delfem2/bv.h"    // include gl
//#include "delfem2/cad2d.h"
//#include "delfem2/sdf.h"
//#include "delfem2/isosurface_stuffing.h"
//#include "delfem2/mathexpeval.h"

namespace py = pybind11;

//////////////////////////////////////////////////////////////////////////////////////////

void init_sampler(py::module &m);
void init_texture(py::module &m);

PYBIND11_MODULE(c_gl, m) {
  m.doc() = "pybind11 delfem2 binding";
  ///////////////////////////////////
  init_sampler(m);
  init_texture(m);
   
 ////////////////////////////////////

  py::class_<CColorMap>(m,"ColorMap")
  .def(py::init<>())
  .def(py::init<double, double, const std::string&>());
  
  ///////////////////////////////////
  // gl misc
  m.def("setSomeLighting",  &setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl",       &setUpGLSL, "compile shader program");
  m.def("glew_init",        &glewInit);
  m.def("draw_sphere",      &DrawSphereAt );
}


