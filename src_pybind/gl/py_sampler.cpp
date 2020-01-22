#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/opengl/gl_framebuffer.h"
#include "delfem2/opengl/glold_smplr.h"

namespace py = pybind11;

// -----------------

py::array_t<float> depth_buffer(CGPUSamplerDrawer& sampler)
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

py::array_t<unsigned char> color_buffer_4byte(CGPUSamplerDrawer& sampler)
{
  std::vector<unsigned char> aRGBA;
  sampler.ExtractFromTexture_Color(aRGBA);
  assert(aRGBA.size()==sampler.nResY*sampler.nResX*4);
  std::vector<size_t> strides = {sizeof(unsigned char)*sampler.nResX*4,sizeof(unsigned char)*4,sizeof(unsigned char)};
  std::vector<size_t> shape = {(size_t)sampler.nResY,(size_t)sampler.nResX,4};
  size_t ndim = 3;
  return py::array(py::buffer_info(aRGBA.data(), sizeof(unsigned char),
                                   py::format_descriptor<unsigned char>::format(),
                                   ndim, shape, strides));
}

/*
py::array_t<float> color_buffer_4float(CGPUSampler& sampler)
{
  assert((int)sampler.aF_RGBA.size()==sampler.nResY*sampler.nResX*4);
  std::vector<size_t> strides = {sizeof(float)*sampler.nResX*4,sizeof(float)*4,sizeof(float)};
  std::vector<size_t> shape = {(size_t)sampler.nResY,(size_t)sampler.nResX,4};
  size_t ndim = 3;
  return py::array(py::buffer_info(sampler.aF_RGBA.data(), sizeof(float),
                                   py::format_descriptor<float>::format(),
                                   ndim, shape, strides));
}
 */


void init_sampler(py::module &m)
{
  // ---------------------------------------
  // FrameBuffer
  py::class_<CFrameBufferManager>(m,"CppFrameBufferManager", "Buffer Class for Depth")
  .def(py::init<>())
  .def("set_buffer_size", &CFrameBufferManager::Init)
  .def("start",           &CFrameBufferManager::Start)
  .def("end",             &CFrameBufferManager::End);
  
  // ---------------------------------------
  // Depth&Color Sampler
  py::class_<CGPUSamplerDrawer>(m,"CppGPUSampler", "sample color and depth in the frame buffer")
  .def(py::init<>())
  .def("draw",       &CGPUSamplerDrawer::Draw)
  .def("minmax_xyz", &CGPUSamplerDrawer::MinMaxXYZ)
  .def("init_gl",    &CGPUSamplerDrawer::InitGL)
  // --------------------------
  .def("init",       &CGPUSamplerDrawer::Init,
       py::arg("size_res_width"),
       py::arg("size_res_height") )
  .def("set_coordinate", &CGPUSamplerDrawer::SetCoord,
       py::arg("len_grid"),
       py::arg("depth_max"),
       py::arg("org"),
       py::arg("dir_prj"),
       py::arg("dir_width"))
  // ---------------
  .def("start",      &CGPUSamplerDrawer::Start)
  .def("end",        &CGPUSamplerDrawer::End)
  .def("get_pos_ray_collide",   &CGPUSamplerDrawer::getGPos)
  .def("set_zero_to_depth",     &CGPUSamplerDrawer::SetZeroToDepth)
  .def("get_depth",             &CGPUSamplerDrawer::GetDepth)
  .def("get_color",             &CGPUSamplerDrawer::GetColor)
  .def_readwrite("bgcolor",     &CGPUSamplerDrawer::bgcolor)
  .def_readwrite("point_color", &CGPUSamplerDrawer::colorPoint)
  .def_readwrite("len_axis",    &CGPUSamplerDrawer::draw_len_axis)
  .def_readwrite("is_draw_tex", &CGPUSamplerDrawer::isDrawTex)
  .def_readwrite("point_size",  &CGPUSamplerDrawer::pointSize);
  
  
  m.def("depth_buffer", &depth_buffer);
  m.def("color_buffer_4byte", &color_buffer_4byte);
//  m.def("color_buffer_4float", &color_buffer_4float);
}
