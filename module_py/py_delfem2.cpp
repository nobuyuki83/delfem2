#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>

#include "delfem2/funcs_glew.h" // have to be included in the beginning
#include "delfem2/depth_v3_gl.h"
#include "delfem2/mshtopoio_gl.h"

#include "delfem2/voxel.h"
#include "delfem2/bv.h"    // include gl
#include "delfem2/cad2d.h"

// TODO:Make a wrapper class of the VoxelGrid?
CMeshElem MeshQuad3D_VoxelGrid(const CVoxelGrid& vg){
  CMeshElem me;
  vg.GetQuad(me.aPos, me.aElem);
  me.elem_type = MESHELEM_QUAD;
  me.ndim = 3;
  return me;
}

class CAxisXYZ {
public:
  CAxisXYZ(): len(1.0){
    line_width = 1.0;
  }
  CAxisXYZ(double len): len(len){
    line_width=1.0;
  }
  void Draw() const{
    glLineWidth(line_width);
    DrawAxis(len);
  }
  std::vector<double> MinMaxXYZ() const{
    std::vector<double> mm(6,0);
    mm[1] = len;  mm[3] = len;  mm[5] = len;
    return mm;
  }
public:
  double len;
  double line_width;
};


//////////////////////////////////////////////////////////////////////////////////////////

namespace py = pybind11;

py::array_t<float> depth_buffer(CDepth& depth)
{
  std::vector<size_t> strides = {sizeof(float)*depth.nResX,sizeof(float)};
  std::vector<size_t> shape = {(size_t)depth.nResY,(size_t)depth.nResX};
  size_t ndim = 2;
  return py::array(py::buffer_info(depth.aZ.data(), sizeof(float),
                                   py::format_descriptor<float>::format(),
                                   ndim, shape, strides));
}

void init_mshtopoio_gl(py::module &m);
void init_fbx(py::module &m);

PYBIND11_MODULE(dfm2, m) {
  m.doc() = "pybind11 delfem2 binding";
  ///////////////////////////////////
  
#ifdef USE_FBX
  init_fbx(m);
#endif

  ///////////////////////////////////
  // axis arrigned boudning box
  py::class_<CBV3D_AABB>(m,"AABB3", "3D axis aligned bounding box class")
  .def(py::init<>())
  .def(py::init<const std::vector<double>&>())
  .def("__str__",            &CBV3D_AABB::str, "print x_min,x_max,y_min,y_max,z_min,z_max")
  .def("minmax_xyz",         &CBV3D_AABB::MinMaxXYZ)
  .def("draw",               &CBV3D_AABB::Draw, "draw edge of the bounding box to opengl")
  .def("set_minmax_xyz",     &CBV3D_AABB::SetMinMaxXYZ)
  .def("add_minmax_xyz",     &CBV3D_AABB::Add_AABBMinMax)
  .def("list_xyz",           &CBV3D_AABB::Point3D_Vox, "corner xyz coords in voxel point order")
  .def("diagonal_length",    &CBV3D_AABB::DiagonalLength, "diagonal length of the bounding box")
  .def("max_length",         &CBV3D_AABB::MaxLength, "diagonal length of the bounding box")
  .def_readwrite("isActive", &CBV3D_AABB::is_active);
  
  py::class_<CAxisXYZ>(m,"AxisXYZ","3D axis class")
  .def(py::init<>())
  .def(py::init<double>(), py::arg("len"))
  .def("draw",                 &CAxisXYZ::Draw)
  .def("minmax_xyz",           &CAxisXYZ::MinMaxXYZ)
  .def_readwrite("len",        &CAxisXYZ::len)
  .def_readwrite("line_width", &CAxisXYZ::line_width);
  
  ///////////////////////////////////
  // mesh
  init_mshtopoio_gl(m);
  
  ///////////////////////////////////
  // voxel
  py::class_<CVoxelGrid>(m, "VoxelGrid", "voxel grid class")
  .def(py::init<>())
  .def("add",&CVoxelGrid::Add,"add voxel at the integer coordinate");
  
  m.def("meshQuad3d_voxelGrid",  &MeshQuad3D_VoxelGrid, "get quad mesh from voxel grid");
  
  ///////////////////////////////////
  // cad
  
  py::class_<CCad2D>(m, "Cad2D", "2D CAD class")
  .def(py::init<>())
  .def("add_square", &CCad2D::Add_Square)
  .def("draw",       &CCad2D::Draw)
  .def("minmax_xyz", &CCad2D::MinMaxXYZ);
  
  ///////////////////////////////////
  // depth
  
  py::class_<CFrameBufferManager>(m,"FrameBufferManager", "Buffer Class for Depth")
  .def(py::init<>())
  .def(py::init<const std::vector<int>&>(),py::arg("win_size"))
  .def("set_buffer_size", &CFrameBufferManager::Init)
  .def("start",           &CFrameBufferManager::Start)
  .def("end",             &CFrameBufferManager::End);
  
  py::class_<CDepth>(m,"Depth","Depth projection class")
  .def(py::init<>())
  .def(py::init<int,int,double,double,const std::vector<double>&,std::vector<double>&,std::vector<double>&>(),
       py::arg("size_res_width"),py::arg("size_res_height"),py::arg("len_grid"),py::arg("depth_max"),
       py::arg("org"), py::arg("dir_prj"), py::arg("dir_width"))
  .def("set_coordinate", &CDepth::SetCoord,
       py::arg("size_res_width"),py::arg("size_res_height"),py::arg("len_grid"),py::arg("depth_max"),
       py::arg("org"), py::arg("dir_prj"), py::arg("dir_width"))
  .def("draw",       &CDepth::Draw)
  .def("minmax_xyz", &CDepth::MinMaxXYZ)
  .def("start",      &CDepth::Start)
  .def("end",        &CDepth::End)
  .def_readwrite("color",  &CDepth::color)
  .def_readwrite("len_axis",  &CDepth::draw_len_axis);
  
  m.def("depth_buffer", &depth_buffer);

  ///////////////////////////////////
  // gl misc
  m.def("setSomeLighting",  &setSomeLighting, "set some lighting that looks good for me");
}
