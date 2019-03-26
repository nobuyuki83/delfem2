#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <map>
#include <deque>

#include "delfem2/funcs_glew.h" // have to be included in the beginning
#include "delfem2/funcs_gl.h"
#include "delfem2/mshtopoio_gl.h"

#include "delfem2/voxel.h"
#include "delfem2/bv.h"    // include gl
#include "delfem2/cad2d.h"

#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"
#include "delfem2/fem.h"

namespace py = pybind11;

//////////////////////////////////////////////////////////////////////////////////////////


void init_polyline(py::module &m);
void init_mshtopoio_gl(py::module &m);
void init_sampler(py::module &m);
void init_fbx(py::module &m);
void init_texture(py::module &m);
void init_rigidbody(py::module &m);
void init_field(py::module &m);

std::tuple<std::vector<double>,std::vector<int>> GetMesh_VoxelGrid
(const CVoxelGrid& vg)
{
  std::vector<double> aXYZ;
  std::vector<int> aQuad;
  vg.GetQuad(aXYZ, aQuad);
  return std::tie(aXYZ,aQuad);
}


std::tuple<py::array_t<double>, py::array_t<int>> GetMesh_Cad
(const CCad2D& cad, double len)
{
  std::vector<double> aXY;
  std::vector<int> aTri;
  cad.Meshing(aXY,aTri, len);
  ////
  py::array_t<double> npXY({(int)aXY.size()/2,2}, aXY.data());
  py::array_t<int> npTri({(int)aTri.size()/3,3}, aTri.data());
  return std::tie(npXY,npTri);
}


void MatrixSquareSparse_SetPattern
(CMatrixSquareSparse& mss,
 const py::array_t<int>& psup_ind,
 const py::array_t<int>& psup)
{
  assert( psup_ind.ndim()  == 1 );
  assert( psup.ndim()  == 1 );
  const int np = mss.m_nblk_col;
  assert( psup_ind.shape()[0] == np+1 );
  mss.SetPattern(psup_ind.data(), psup_ind.shape()[0],
                 psup.data(),     psup.shape()[0]);
}

void MatrixSquareSparse_SetFixBC
(CMatrixSquareSparse& mss,
 const py::array_t<int>& flagbc)
{
  mss.SetBoundaryCondition(flagbc.data(),flagbc.shape()[0]);
}


void PrecondILU0
(CPreconditionerILU&  mat_ilu,
 const CMatrixSquareSparse& mss)
{
  mat_ilu.Initialize_ILU0(mss);
}

void PyMergeLinSys_Poission2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double alpha, double source,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_Poission2D(mss, (double*)buff_vecb.ptr,
                         alpha, source,
                         aXY.data(), aXY.shape()[0],
                         aTri.data(), aTri.shape()[0],
                         aVal.data());
}


void PySolve_PCG
(py::array_t<double>& vec_b,
 py::array_t<double>& vec_x,
 double conv_ratio, double iteration,
 const CMatrixSquareSparse& mat_A,
 const CPreconditionerILU& ilu_A)
{
//  std::cout << "solve pcg" << std::endl;
  auto buff_vecb = vec_b.request();
  auto buff_vecx = vec_x.request();
  Solve_PCG((double*)buff_vecb.ptr,
            (double*)buff_vecx.ptr,
            conv_ratio,iteration,
            mat_A,ilu_A);
}

void PySortIndexedArray
(py::array_t<int>& psup_ind,
 py::array_t<int>& psup)
{
//  std::cout << "hoge " << psup_ind.size() << " " << psup.size() << std::endl;
  auto buff_psup = psup.request();
  SortIndexedArray(psup_ind.data(), psup_ind.shape()[0]-1, (int*)buff_psup.ptr);
}

void PyCad_SetBCFlagEdge
(py::array_t<int>& vec_bc,
 const py::array_t<double>& aXY,
 const py::array_t<int>& np_ie,
 const CCad2D& cad,
 int iflag,
 double torelance)
{
//  std::cout << "cad_SetBCFlagEdge" << std::endl;
  auto buff_bc = vec_bc.request();
  std::vector<int> aIE(np_ie.data(),np_ie.data()+np_ie.shape()[0]);
  cad.setBCFlagEdge((int*)buff_bc.ptr,
                    aXY.data(), aXY.shape()[0],
                    aIE,iflag,torelance);
}


PYBIND11_MODULE(dfm2, m) {
  m.doc() = "pybind11 delfem2 binding";
  ///////////////////////////////////
  
#ifdef USE_FBX
  init_fbx(m);
#endif

  ///////////////////////////////////
  init_mshtopoio_gl(m);
  init_sampler(m);
  init_polyline(m);
  init_texture(m);
  init_rigidbody(m);
  init_field(m);
  
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
  .def("center",             &CBV3D_AABB::Center, "center position")
  .def_readwrite("isActive", &CBV3D_AABB::is_active);
  
  py::class_<CAxisXYZ>(m,"AxisXYZ","3D axis class")
  .def(py::init<>())
  .def(py::init<double>(), py::arg("len"))
  .def("draw",                 &CAxisXYZ::Draw)
  .def("minmax_xyz",           &CAxisXYZ::MinMaxXYZ)
  .def_readwrite("len",        &CAxisXYZ::len)
  .def_readwrite("line_width", &CAxisXYZ::line_width);
  
  ///////////////////////////////////
  // voxel
  py::class_<CVoxelGrid>(m, "VoxelGrid", "voxel grid class")
  .def(py::init<>())
  .def("add",&CVoxelGrid::Add,"add voxel at the integer coordinate");
  
  m.def("getmesh_voxelgrid",&GetMesh_VoxelGrid);
  m.def("getMesh_cad",&GetMesh_Cad);
  
  ///////////////////////////////////
  // cad
  py::class_<CCad2D>(m, "Cad2D", "2D CAD class")
  .def(py::init<>())
  .def("add_polygon", &CCad2D::AddPolygon)
  .def("draw",       &CCad2D::Draw)
  .def("mouse",      &CCad2D::Mouse)
  .def("motion",      &CCad2D::Motion)
  .def("minmax_xyz", &CCad2D::MinMaxXYZ)
  .def("meshing",&CCad2D::Meshing);
  
  
  py::class_<CColorMap>(m,"ColorMap")
  .def(py::init<>())
  .def(py::init<double, double, const std::string&>());
  
  py::class_<CMatrixSquareSparse>(m,"MatrixSquareSparse")
  .def(py::init<>())
  .def("initialize", &CMatrixSquareSparse::Initialize)
  .def("setZero",    &CMatrixSquareSparse::SetZero);
  
  py::class_<CPreconditionerILU>(m,"PreconditionerILU")
  .def(py::init<>())
  .def("ilu_decomp", &CPreconditionerILU::DoILUDecomp)
  .def("set_value", &CPreconditionerILU::SetValueILU);
//  .def(py::init<const CPreconditionerILU&>);
  
  m.def("matrixSquareSparse_setPattern", &MatrixSquareSparse_SetPattern);
  m.def("matrixSquareSparse_setFixBC", &MatrixSquareSparse_SetFixBC);
  m.def("precond_ilu0",  &PrecondILU0);
  m.def("mergeLinSys_poission2D", &PyMergeLinSys_Poission2D);
  m.def("linsys_solve_pcg", &PySolve_PCG);
  m.def("sortIndexedArray", &PySortIndexedArray);
  m.def("cad_setBCFlagEdge",&PyCad_SetBCFlagEdge);


  ///////////////////////////////////
  // gl misc
  m.def("setSomeLighting",  &setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl", setUpGLSL, "compile shader program");
  m.def("glew_init", glewInit);
  m.def("draw_sphere", DrawSphereAt );
}
