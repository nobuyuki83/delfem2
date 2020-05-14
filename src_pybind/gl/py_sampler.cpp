#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/opengl/r2tglo_glold.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

// -----------------

py::array_t<float> depth_buffer(dfm2::opengl::CRender2Tex_DrawOldGL& sampler)
{
  std::vector<float> aZ;
  sampler.ExtractFromTexture_Depth(aZ);
  assert(aZ.size()==sampler.nResY*sampler.nResX);
  std::vector<size_t> strides = {sizeof(float)*sampler.nResX,sizeof(float)};
  std::vector<size_t> shape = {(size_t)sampler.nResY,(size_t)sampler.nResX};
  size_t ndim = 2;
  return py::array(py::buffer_info(aZ.data(), sizeof(float),
                                   py::format_descriptor<float>::format(),
                                   ndim, shape, strides));
}

py::array_t<unsigned char> color_buffer_4byte(dfm2::opengl::CRender2Tex_DrawOldGL& sampler)
{
  std::vector<unsigned char> aRGBA;
  sampler.ExtractFromTexture_RGBA8UI(aRGBA);
  assert(aRGBA.size()==sampler.nResY*sampler.nResX*4);
  std::vector<size_t> strides = {sizeof(unsigned char)*sampler.nResX*4,sizeof(unsigned char)*4,sizeof(unsigned char)};
  std::vector<size_t> shape = {(size_t)sampler.nResY,(size_t)sampler.nResX,4};
  size_t ndim = 3;
  return py::array(py::buffer_info(aRGBA.data(), sizeof(unsigned char),
                                   py::format_descriptor<unsigned char>::format(),
                                   ndim, shape, strides));
}

py::array_t<float> color_buffer_4float(dfm2::opengl::CRender2Tex_DrawOldGL& sampler)
{
  std::vector<float> aRGBA;
  sampler.ExtractFromTexture_RGBA32F(aRGBA);
  assert(aRGBA.size()==sampler.nResY*sampler.nResX*4);
  std::vector<size_t> strides = {sizeof(float)*sampler.nResX*4,sizeof(float)*4,sizeof(float)};
  std::vector<size_t> shape = {(size_t)sampler.nResY,(size_t)sampler.nResX,4};
  size_t ndim = 3;
  return py::array(py::buffer_info(aRGBA.data(), sizeof(float),
                                   py::format_descriptor<float>::format(),
                                   ndim, shape, strides));
}


void init_sampler(py::module &m)
{
  // ---------------------------------------
  // Depth&Color Sampler
  py::class_<dfm2::opengl::CRender2Tex_DrawOldGL>(m,
                                                  "CppGPUSampler",
                                                  "sample color and depth in the frame buffer")
  .def(py::init<>())
  .def("draw",       &dfm2::opengl::CRender2Tex_DrawOldGL::Draw)
  .def("init_gl",    &dfm2::opengl::CRender2Tex_DrawOldGL::InitGL)
  .def("minmax_xyz", &dfm2::opengl::CRender2Tex_DrawOldGL::AABBVec3)
  //
  .def("set_texture_property",
       &dfm2::opengl::CRender2Tex_DrawOldGL::SetTextureProperty,
       py::arg("size_res_width"),
       py::arg("size_res_height"),
       py::arg("is_rgba_8ui") )
  .def("set_coordinate",
       &dfm2::opengl::CRender2Tex_DrawOldGL::SetCoord,
       py::arg("len_grid"),
       py::arg("depth_max"),
       py::arg("org"),
       py::arg("dir_prj"),
       py::arg("dir_width"))
  .def("start",      &dfm2::opengl::CRender2Tex_DrawOldGL::Start)
  .def("end",        &dfm2::opengl::CRender2Tex_DrawOldGL::End)
  .def("get_pos_ray_collide",   &dfm2::opengl::CRender2Tex_DrawOldGL::getGPos)
  .def("set_zero_to_depth",     &dfm2::opengl::CRender2Tex_DrawOldGL::SetZeroToDepth)
  .def("get_depth",             &dfm2::opengl::CRender2Tex_DrawOldGL::GetDepth)
  .def("get_color",             &dfm2::opengl::CRender2Tex_DrawOldGL::GetColor)
  .def_readwrite("point_color", &dfm2::opengl::CRender2Tex_DrawOldGL::colorPoint)
  .def_readwrite("len_axis",    &dfm2::opengl::CRender2Tex_DrawOldGL::draw_len_axis)
  .def_readwrite("is_draw_tex", &dfm2::opengl::CRender2Tex_DrawOldGL::isDrawTex)
  .def_readwrite("point_size",  &dfm2::opengl::CRender2Tex_DrawOldGL::pointSize);
    
  m.def("depth_buffer", &depth_buffer);
  m.def("color_buffer_4byte", &color_buffer_4byte);
  m.def("color_buffer_4float", &color_buffer_4float);
}
