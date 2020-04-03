#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/opengl/tex_gl.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

dfm2::opengl::CTexture GetTextureFromNumpy(const py::array_t<unsigned char>& a){
  assert(a.ndim()==3);
  assert(a.shape()[2] == 3);
  const int h = a.shape()[0];
  const int w = a.shape()[1];
  dfm2::opengl::CTexture tex(w,h,a.data(),"bgr");
  return tex;
}

void init_texture(py::module &m) {
  py::class_<dfm2::opengl::CTexture>(m, "Texture")
      .def(py::init<>())
      .def("minmax_xyz", &dfm2::opengl::CTexture::MinMaxXYZ)
      //
      .def_readonly("width", &dfm2::opengl::CTexture::w)
      .def_readonly("height", &dfm2::opengl::CTexture::h)
      .def("draw", &dfm2::opengl::CTexture::Draw_oldGL)
      .def("init_gl", &dfm2::opengl::CTexture::InitGL)
      .def("set_minmax_xy", &dfm2::opengl::CTexture::SetMinMaxXY);


  m.def("get_texture", &GetTextureFromNumpy);
}
