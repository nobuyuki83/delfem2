#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/mshtopoio_gl.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/dyntri_v3.h"

namespace py = pybind11;

/*
CMeshElem HightMap(py::array_t<double> input0){
  py::buffer_info buf1 = input0.request();
  if (buf1.ndim != 2 ){
    throw std::runtime_error("Number of dimensions must be one");
  }
  const int nx = buf1.shape[0];
  const int ny = buf1.shape[1];
  CMeshElem msh;
  msh.elem_type = MESHELEM_QUAD;
  msh.ndim = 3;
  {
    double *ptr1 = (double *) buf1.ptr;
    std::vector<double> aXY;
    MeshQuad2D_Grid(aXY, msh.aElem,
                   nx-1, ny-1);
    const int np = aXY.size()/2;
    msh.aPos.resize(np*3);
    for(int iy=0;iy<ny;++iy){
      for(int ix=0;ix<nx;++ix){
        int ip = iy*nx+ix;
        msh.aPos[ip*3+0] = aXY[ip*2+0];
        msh.aPos[ip*3+1] = aXY[ip*2+1];
        msh.aPos[ip*3+2] = ptr1[ip];
      }
    }
  }
  return msh;
}
 */

std::tuple<std::vector<double>,std::vector<int>>
PyGetMesh_Grid
(int mx, int my)
{
  std::vector<double> aXY;
  std::vector<int> aQuad;
  MeshQuad2D_Grid(aXY, aQuad,
                  mx-1, my-1);
  return std::tie(aXY,aQuad);
}

std::tuple<std::vector<double>,std::vector<int>> ReadMesh_Ply
(const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<int> aTri;
  Read_Ply(fname, aXYZ, aTri);
  return std::tie(aXYZ,aTri);
}

std::tuple<std::vector<double>,std::vector<int>> ReadMesh_Obj
(const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<int> aTri;
  Read_Obj(fname, aXYZ, aTri);
  return std::tie(aXYZ,aTri);
}

void PyDrawMesh_FaceNorm
(const py::array_t<double>& pos,
 const py::array_t<int>& elm)
{
  assert(pos.ndim()==2);
  assert(elm.ndim()==2);
  const auto shape_pos = pos.shape();
  const auto shape_elm = elm.shape();
  if( shape_pos[1] == 3 ){ // 3D Mesh
    if( shape_elm[1] == 3 ){  DrawMeshTri3D_FaceNorm(pos.data(), elm.data(), shape_elm[0]); }
    if( shape_elm[1] == 4 ){  DrawMeshQuad3D_FaceNorm(pos.data(), elm.data(), shape_elm[0]); }
  }
}

void PyDrawMesh_Edge
(const py::array_t<double>& pos,
 const py::array_t<int>& elm)
{
  assert(pos.ndim()==2);
  assert(elm.ndim()==2);
  const auto shape_pos = pos.shape();
  const auto shape_elm = elm.shape();
  if( shape_pos[1] == 3 ){ // 3D Mesh
    if( shape_elm[1] == 3 ){  DrawMeshTri3D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( shape_elm[1] == 4 ){  DrawMeshQuad3D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
  if( shape_pos[1] == 2 ){ // 3D Mesh
    if( shape_elm[1] == 3 ){  DrawMeshTri2D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( shape_elm[1] == 4 ){  DrawMeshQuad2D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
}

std::tuple<std::vector<double>,std::vector<int>>
PySubviv
(const std::vector<double>& aXYZ0, const std::vector<int>& aQuad0)
{
  std::vector<int> aQuad1;
  std::vector<int> aEdgeFace0;
  std::vector<int> psupIndQuad0, psupQuad0;
  QuadSubdiv(aQuad1,
             psupIndQuad0,psupQuad0, aEdgeFace0,
             aQuad0, aXYZ0.size()/3);
  ///////
  std::vector<double> aXYZ1;
  SubdivisionPoints_QuadCatmullClark(aXYZ1,
                                     aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,aQuad0,aXYZ0);
  return std::tie(aXYZ1,aQuad1);
}

std::tuple<std::vector<double>,std::vector<int>,std::vector<int>,std::vector<int>>
Triangulation
(const std::vector< std::vector<double> >& aaXY,
 double edge_length)
{
//  CTriangulationOutput out;
//  std::vector< std::vector<double> > aaXY;
//  aaXY.push_back(aXY);
  std::vector<int> aPtrVtxInd;
  std::vector<int> aVtxInd;
  std::vector<int> aElm;
  std::vector<double> aPos;
  GenerateTesselation2(aElm, aPos,
                       aPtrVtxInd, aVtxInd,
                       edge_length, true, aaXY);
  return std::tie(aPos,aElm, aPtrVtxInd,aVtxInd);
}

std::tuple<py::array_t<int>, py::array_t<int>>
GetPsup(const py::array_t<int>& elm, int npoint)
{
  std::vector<int> psup_ind, psup;
  makeOneRingNeighborhood(psup_ind, psup,
                          elm.data(), elm.shape()[0], elm.shape()[1], npoint);
  py::array_t<int> np_psup_ind((pybind11::size_t)psup_ind.size(), psup_ind.data());
  py::array_t<int> np_psup((pybind11::size_t)psup.size(), psup.data());
  return std::tie(np_psup_ind, np_psup);
}


///////////////////////////////////////////////////////////////////////////////////////

void init_mshtopoio_gl(py::module &m){
  py::enum_<MESHELEM_TYPE>(m, "MeshElemType")
  .value("Tri",     MESHELEM_TYPE::MESHELEM_TRI)
  .value("Quad",    MESHELEM_TYPE::MESHELEM_QUAD)
  .value("Tet",     MESHELEM_TYPE::MESHELEM_TET)
  .value("Pyramid", MESHELEM_TYPE::MESHELEM_PYRAMID)
  .value("Wedge",   MESHELEM_TYPE::MESHELEM_WEDGE)
  .value("Hex",     MESHELEM_TYPE::MESHELEM_HEX)
  .export_values();
  
  /*
  py::class_<CMeshElem>(m, "MeshElem")
  .def(py::init<>())
  .def(py::init<const std::string&>(),"open the file in the path")
  .def("draw", &CMeshElem::Draw) //
  .def("minmax_xyz", &CMeshElem::AABB3_MinMax)
  .def("read", &CMeshElem::Read)
  .def("write_obj",&CMeshElem::Write_Obj)
  .def("drawFace_elemWiseNorm", &CMeshElem::DrawFace_ElemWiseNorm)
  .def("drawEdge", &CMeshElem::DrawEdge)
  .def("scaleXYZ", &CMeshElem::ScaleXYZ)
  .def("subdiv",   &CMeshElem::Subdiv)
  .def_readonly("listElem", &CMeshElem::aElem)
  .def_readonly("listPos",  &CMeshElem::aPos)
  .def_readonly("elemType", &CMeshElem::elem_type)
  .def_readonly("nDim",     &CMeshElem::ndim)
  .def_readwrite("color_face",  &CMeshElem::color_face)
  .def_readwrite("is_draw_edge", &CMeshElem::is_draw_edge);
   */
  
  
  
  m.def("read_nastran_triangle",&Read_MeshTri3D_Nas_CMeshElem);
  
  py::class_<CMeshMultiElem>(m,"MeshMultiElem")
  .def(py::init<>())
  .def("read_obj", &CMeshMultiElem::ReadObj)
  .def("minmax_xyz", &CMeshMultiElem::AABB3_MinMax)
  .def("draw",&CMeshMultiElem::Draw)
  .def("scaleXYZ",&CMeshMultiElem::ScaleXYZ)
  .def("translateXYZ",&CMeshMultiElem::TranslateXYZ);
  
  m.def("triangulation",&Triangulation);
  
  m.def("triangulation",&Triangulation,
        py::arg("aXY"),
        py::arg("edge_length")=0.03);
  
//  m.def("hight_map", &HightMap);
  
  m.def("get_psup",&GetPsup);

  m.def("read_ply",&ReadMesh_Ply);
  m.def("read_obj",&ReadMesh_Obj);
  
  m.def("draw_mesh_facenorm", &PyDrawMesh_FaceNorm);
  m.def("draw_mesh_edge", &PyDrawMesh_Edge);
  m.def("subdiv",&PySubviv);
  m.def("get_mesh_grid",&PyGetMesh_Grid);
  
}
