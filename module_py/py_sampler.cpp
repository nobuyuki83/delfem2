#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/funcs_glew.h"
#include "delfem2/gpusampler.h"

namespace py = pybind11;

///////////////////////////////////

py::array_t<float> depth_buffer(CGPUSampler& sampler)
{
  assert((int)sampler.aZ.size()==sampler.nResY*sampler.nResX);
  std::vector<size_t> strides = {sizeof(float)*sampler.nResX,sizeof(float)};
  std::vector<size_t> shape = {(size_t)sampler.nResY,(size_t)sampler.nResX};
  size_t ndim = 2;
  return py::array(py::buffer_info(sampler.aZ.data(), sizeof(float),
                                   py::format_descriptor<float>::format(),
                                   ndim, shape, strides));
}

py::array_t<unsigned char> color_buffer(CGPUSampler& sampler)
{
  assert((int)sampler.aUC_RGBA.size()==sampler.nResY*sampler.nResX*4);
  std::vector<size_t> strides = {sizeof(unsigned char)*sampler.nResX*4,sizeof(unsigned char)*4,sizeof(unsigned char)};
  std::vector<size_t> shape = {(size_t)sampler.nResY,(size_t)sampler.nResX,4};
  size_t ndim = 3;
  return py::array(py::buffer_info(sampler.aUC_RGBA.data(), sizeof(unsigned char),
                                   py::format_descriptor<unsigned char>::format(),
                                   ndim, shape, strides));
}




void init_sampler(py::module &m){
  ///////////////////////////////////
  // FrameBuffer
  py::class_<CFrameBufferManager>(m,"FrameBufferManager", "Buffer Class for Depth")
  .def(py::init<>())
  .def(py::init<const std::vector<int>&,std::string,bool>(),
       py::arg("win_size"),py::arg("is_color"),py::arg("is_depth"))
  .def("set_buffer_size", &CFrameBufferManager::Init)
  .def("start",           &CFrameBufferManager::Start)
  .def("end",             &CFrameBufferManager::End);
  
  //////////////////////////////////
  // Depth&Color Sampler
  py::class_<CGPUSampler>(m,"GPUSampler","sample color and depth in the frame buffer")
  .def(py::init<>())
  .def("init",       &CGPUSampler::Init,
       py::arg("size_res_width"),py::arg("size_res_height"),py::arg("format_color"),py::arg("is_depth"))
  .def("set_coordinate", &CGPUSampler::SetCoord,
       py::arg("len_grid"),py::arg("depth_max"), py::arg("org"), py::arg("dir_prj"), py::arg("dir_width"))
  .def("draw",       &CGPUSampler::Draw)
  .def("minmax_xyz", &CGPUSampler::MinMaxXYZ)
  .def("start",      &CGPUSampler::Start)
  .def("end",        &CGPUSampler::End)
  .def("init_gl",    &CGPUSampler::LoadTex)
  .def_readwrite("bgcolor", &CGPUSampler::bgcolor)
  .def_readwrite("color",  &CGPUSampler::color)
  .def_readwrite("len_axis",  &CGPUSampler::draw_len_axis);
  
  m.def("depth_buffer", &depth_buffer);
  m.def("color_buffer", &color_buffer);
}
