#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/gl.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "delfem2/gl_tex.h"

namespace py = pybind11;


CTexture GetTextureFromNumpy(const py::array_t<unsigned char>& a){
  assert(a.ndim()==3);
  assert(a.shape()[2] == 3);
  const int h = a.shape()[0];
  const int w = a.shape()[1];
  CTexture tex(w,h,a.data(),"bgr");
  return tex;
}

void init_texture(py::module &m){
  py::class_<CTexture>(m,"Texture")
  .def(py::init<>())
  .def("draw",&CTexture::Draw)
  .def("init_gl",&CTexture::LoadTex)
  .def("minmax_xyz",&CTexture::MinMaxXYZ);
  
  m.def("get_texture", &GetTextureFromNumpy);
}
