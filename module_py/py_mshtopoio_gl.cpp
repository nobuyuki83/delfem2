#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/mshtopoio_gl.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/dyntri.h"
#include "delfem2/dyntri_v2.h"
#include "delfem2/dyntri_v3.h"

namespace py = pybind11;

std::tuple<py::array_t<double>,py::array_t<int>> PyMeshQuad2D_Grid
(int mx, int my)
{
  std::vector<double> aXY;
  std::vector<int> aQuad;
  MeshQuad2D_Grid(aXY, aQuad,
                  mx-1, my-1);
  py::array_t<double> npXY({(int)aXY.size()/2,2}, aXY.data());
  py::array_t<int> npQuad({(int)aQuad.size()/4,4}, aQuad.data());
  return std::forward_as_tuple(npXY,npQuad);
}

std::tuple<py::array_t<double>,py::array_t<int>> PyMeshTri3D_ReadPly
(const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<int> aTri;
  Read_Ply(fname, aXYZ, aTri);
  py::array_t<double> np_XYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<int> np_Tri({(int)aTri.size()/3,3}, aTri.data());
  return std::forward_as_tuple(np_XYZ,np_Tri);
}

std::tuple<py::array_t<double>,py::array_t<int>> PyMeshTri3D_ReadObj
(const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<int> aTri;
  Read_Obj(fname, aXYZ, aTri);
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<int> npTri({(int)aTri.size()/3,3}, aTri.data());
  return std::forward_as_tuple(npXYZ,npTri);
}

void PyMeshTri3D_WriteObj
(const std::string& fname,
 const py::array_t<double>& aXYZ,
 const py::array_t<int>& aTri)
{
  Write_Obj(fname,
            aXYZ.data(), aXYZ.shape()[0],
            aTri.data(), aTri.shape()[0]);
}

std::tuple<py::array_t<double>,py::array_t<int>> PyMeshTri3D_ReadNastran
(const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<int> aTri;
  Read_MeshTri3D_Nas(aXYZ, aTri, fname.c_str());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<int> npTri({(int)aTri.size()/3,3}, aTri.data());
  return std::forward_as_tuple(npXYZ,npTri);
}

std::tuple<py::array_t<double>,py::array_t<int>> PyMeshQuad3D_Subviv
(const py::array_t<double>& aXYZ0, const py::array_t<int>& aQuad0)
{
  std::vector<int> aQuad1;
  std::vector<int> psupIndQuad0, psupQuad0;
  std::vector<int> aEdgeFace0;
  QuadSubdiv(aQuad1,
             psupIndQuad0,psupQuad0, aEdgeFace0,
             aQuad0.data(), aQuad0.shape()[0], aXYZ0.shape()[0]);
  ///////
  std::vector<double> aXYZ1;
  SubdivisionPoints_QuadCatmullClark(aXYZ1,
                                     aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,
                                     aQuad0.data(),aQuad0.shape()[0],
                                     aXYZ0.data(), aXYZ0.shape()[0]);
  py::array_t<double> npXYZ1({(int)aXYZ1.size()/3,3}, aXYZ1.data());
  py::array_t<int> npQuad1({(int)aQuad1.size()/4,4}, aQuad1.data());
  return std::forward_as_tuple(npXYZ1,npQuad1);
}

std::tuple<py::array_t<double>,py::array_t<int>> PyMeshHex3D_Subviv
(const py::array_t<double>& aXYZ0, const py::array_t<int>& aHex0)
{
  std::vector<int> aHex1;
  std::vector<int> psupIndHex0, psupHex0;
  std::vector<int> aQuadHex0;
  HexSubdiv(aHex1,
            psupIndHex0, psupHex0,
            aQuadHex0,
            ///
            aHex0.data(), aHex0.shape()[0], aXYZ0.shape()[0]);
  ///////
  std::vector<double> aXYZ1;
  SubdivisionPoints_Hex(aXYZ1,
                        psupIndHex0,psupHex0,aQuadHex0,
                        aHex0.data(), aHex0.shape()[0],
                        aXYZ0.data(), aXYZ0.shape()[0]);
  py::array_t<double> npXYZ1({(int)aXYZ1.size()/3,3}, aXYZ1.data());
  py::array_t<int> npHex1({(int)aHex1.size()/8,8}, aHex1.data());
  return std::forward_as_tuple(npXYZ1,npHex1);
}

void PyMeshDynTri3D_Initialize
(CMeshDynTri3D& mesh,
 const py::array_t<double>& po,
 const py::array_t<int>& tri)
{
  mesh.Initialize(po.data(), po.shape()[0], po.shape()[1],
                  tri.data(), tri.shape()[0]);
}

void PyMeshDynTri2D_Initialize
(CMeshDynTri2D& mesh,
 const py::array_t<double>& po,
 const py::array_t<int>& tri)
{
  assert(po.shape()[1]==2);
  mesh.Initialize(po.data(), po.shape()[0],
                  tri.data(), tri.shape()[0]);
}


//////////////////////////////////////////////////////////////////////////////////////////

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
    if( shape_elm[1] == 8 ){  DrawMeshHex3D_FaceNorm(pos.data(), elm.data(), shape_elm[0]); }
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
    if( shape_elm[1] == 8 ){  DrawMeshHex3D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]);  }
  }
  if( shape_pos[1] == 2 ){ // 2D Mesh
    if( shape_elm[1] == 3 ){  DrawMeshTri2D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( shape_elm[1] == 4 ){  DrawMeshQuad2D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
}



std::tuple<std::vector<double>,std::vector<int>,std::vector<int>,std::vector<int>>
PyTriangulation
(const std::vector< std::vector<double> >& aaXY,
 double edge_length)
{
  std::vector<int> loop1_ind,loop1;
  std::vector<double> aXY0;
  JArray_FromVecVec_XY(loop1_ind,aXY0,
                       aaXY);
  assert( CheckInputBoundaryForTriangulation(loop1_ind,aXY0) );
  /////
  loop1.resize(aXY0.size()/2);
  for(unsigned int ip=0;ip<aXY0.size()/2;++ip){ loop1[ip] = ip; }
  /////
  FixLoopOrientation(loop1,
                     loop1_ind,aXY0);
  ResamplingLoop(loop1_ind,loop1,aXY0,
                 edge_length );
  std::vector<double> aPos;
  std::vector<int> aElm;
  CInputTriangulation_Uniform param(1.0);
  Triangulation(aElm, aPos,
                loop1_ind, loop1,
                edge_length, param, aXY0);
  return std::forward_as_tuple(aPos,aElm, loop1_ind,loop1);
}

std::tuple<py::array_t<int>, py::array_t<int>>
PyJArray_MeshPsup(const py::array_t<int>& elm, int npoint)
{
  std::vector<int> psup_ind, psup;
  JArray_MeshOneRingNeighborhood(psup_ind, psup,
                                      elm.data(), elm.shape()[0], elm.shape()[1], npoint);
  py::array_t<int> np_psup_ind((pybind11::size_t)psup_ind.size(), psup_ind.data());
  py::array_t<int> np_psup((pybind11::size_t)psup.size(), psup.data());
  return std::forward_as_tuple(np_psup_ind, np_psup);
}

void PyJArray_Sort
(py::array_t<int>& psup_ind,
 py::array_t<int>& psup)
{
  //  std::cout << "hoge " << psup_ind.size() << " " << psup.size() << std::endl;
  auto buff_psup = psup.request();
  JArray_Sort(psup_ind.data(), psup_ind.shape()[0]-1, (int*)buff_psup.ptr);
}

std::tuple<py::array_t<int>, py::array_t<int>>
PyJArray_AddDiagonal(py::array_t<int>& psup_ind0, py::array_t<int>& psup0)
{
  std::vector<int> psup_ind, psup;
  JArray_AddDiagonal(psup_ind,psup,
                          psup_ind0.data(),psup_ind0.shape()[0], psup0.data(),psup0.shape()[0]);
  py::array_t<int> np_psup_ind((pybind11::size_t)psup_ind.size(), psup_ind.data());
  py::array_t<int> np_psup((pybind11::size_t)psup.size(), psup.data());
  return std::forward_as_tuple(np_psup_ind, np_psup);
}


py::array_t<int> GetElemQuad_DihedralTri
(py::array_t<int>& aTri,
 int np)
{
  assert( aTri.shape()[1] == 3 );
  const int nTri = aTri.shape()[0];
  std::vector<int> aElemSurRel;
  makeSurroundingRelationship(aElemSurRel,
                              aTri.data(), nTri,
                              MESHELEM_TRI, np);
  ////
  std::vector<int> aQuad;
  for(int itri=0; itri<nTri; ++itri){
    for(int iedtri=0;iedtri<3;++iedtri){
      int jtri = aElemSurRel[itri*6+iedtri*2+0];
      if( jtri == -1 ) continue;
      if( jtri < itri ) continue;
      int jedtri = aElemSurRel[itri*6+iedtri*2+1];
      assert( itri == aElemSurRel[jtri*6+jedtri*2+0] );
      int ipo0 = aTri.at(itri,iedtri);
      int ipo1 = aTri.at(jtri,jedtri);
      int ipo2 = aTri.at(itri,(iedtri+1)%3);
      int ipo3 = aTri.at(itri,(iedtri+2)%3);
      assert( aTri.at(jtri,(jedtri+2)%3) == ipo2 );
      assert( aTri.at(jtri,(jedtri+1)%3) == ipo3 );
      aQuad.push_back(ipo0);
      aQuad.push_back(ipo1);
      aQuad.push_back(ipo2);
      aQuad.push_back(ipo3);
    }
  }
  ////
  py::array_t<int> npQuad({(int)aQuad.size()/4,4}, aQuad.data());
  return npQuad;
}

std::tuple<double,double>
PyQuality_MeshTri2D
(const py::array_t<double>& np_xy,
 const py::array_t<int>& np_tri)
{
  double max_aspect;
  double min_area;
  Quality_MeshTri2D(max_aspect, min_area,
                    np_xy.data(),
                    np_tri.data(), np_tri.shape()[0]);
   return std::tie(max_aspect, min_area);
}



///////////////////////////////////////////////////////////////////////////////////////

void init_mshtopoio_gl(py::module &m){
  py::enum_<MESHELEM_TYPE>(m, "MESH_ELEM_TYPE")
  .value("TRI",     MESHELEM_TYPE::MESHELEM_TRI)
  .value("QUAD",    MESHELEM_TYPE::MESHELEM_QUAD)
  .value("TET",     MESHELEM_TYPE::MESHELEM_TET)
  .value("PYRAMID", MESHELEM_TYPE::MESHELEM_PYRAMID)
  .value("WEDGE",   MESHELEM_TYPE::MESHELEM_WEDGE)
  .value("HEX",     MESHELEM_TYPE::MESHELEM_HEX)
  .export_values();
  
  py::class_<CMeshMultiElem>(m,"MeshMultiElem")
  .def(py::init<>())
  .def("read_obj", &CMeshMultiElem::ReadObj)
  .def("minmax_xyz", &CMeshMultiElem::AABB3_MinMax)
  .def("draw",&CMeshMultiElem::Draw)
  .def("scaleXYZ",&CMeshMultiElem::ScaleXYZ)
  .def("translateXYZ",&CMeshMultiElem::TranslateXYZ);
  
  py::class_<CMeshDynTri3D>(m, "CppMeshDynTri3D")
  .def(py::init<>())
  .def("draw",              &CMeshDynTri3D::draw)
  .def("draw_face",         &CMeshDynTri3D::Draw_FaceNorm)
  .def("draw_edge",         &CMeshDynTri3D::Draw_Edge)
  .def("check",             &CMeshDynTri3D::Check)
  .def("ntri",              &CMeshDynTri3D::nTri)
  .def("delete_tri_edge",   &CMeshDynTri3D::DeleteTriEdge)
  .def("minmax_xyz",        &CMeshDynTri3D::MinMax_XYZ)
  .def("insert_point_elem", &CMeshDynTri3D::insertPointElem)
  .def("delaunay_around_point", &CMeshDynTri3D::DelaunayAroundPoint);
  
  py::class_<CMeshDynTri2D>(m, "CppMeshDynTri2D")
  .def(py::init<>())
  .def("draw",              &CMeshDynTri2D::draw)
  .def("draw_face",         &CMeshDynTri2D::Draw_FaceNorm)
  .def("draw_edge",         &CMeshDynTri2D::Draw_Edge)
  .def("check",             &CMeshDynTri2D::Check)
  .def("ntri",              &CMeshDynTri2D::nTri)
  .def("delete_tri_edge",   &CMeshDynTri2D::DeleteTriEdge)
  .def("minmax_xyz",        &CMeshDynTri2D::MinMax_XYZ)
  .def("insert_point_elem", &CMeshDynTri2D::insertPointElem)
  .def("delaunay_around_point", &CMeshDynTri2D::DelaunayAroundPoint);
  
  m.def("meshdyntri3d_initialize",&PyMeshDynTri3D_Initialize);
  m.def("meshdyntri2d_initialize",&PyMeshDynTri2D_Initialize);
  
  m.def("meshtri3d_read_ply",     &PyMeshTri3D_ReadPly,     py::return_value_policy::move);
  m.def("meshtri3d_read_obj",     &PyMeshTri3D_ReadObj,     py::return_value_policy::move);
  m.def("meshtri3d_read_nastran", &PyMeshTri3D_ReadNastran, py::return_value_policy::move);
  m.def("meshtri3d_write_obj",    &PyMeshTri3D_WriteObj);
  m.def("meshquad3d_subdiv",      &PyMeshQuad3D_Subviv,     py::return_value_policy::move);
  m.def("meshhex3d_subdiv",       &PyMeshHex3D_Subviv,      py::return_value_policy::move);
  m.def("meshquad2d_grid",        &PyMeshQuad2D_Grid,       py::return_value_policy::move);
  
  m.def("triangulation",&PyTriangulation,
        py::arg("aXY"),
        py::arg("edge_length")=0.03);
  
  m.def("jarray_mesh_psup",    &PyJArray_MeshPsup,    py::return_value_policy::move);
  m.def("jarray_add_diagonal", &PyJArray_AddDiagonal, py::return_value_policy::move);
  m.def("jarray_sort",         &PyJArray_Sort);
  m.def("elemQuad_dihedralTri",&GetElemQuad_DihedralTri);
  m.def("quality_meshTri2D",   &PyQuality_MeshTri2D);
  m.def("draw_mesh_facenorm",  &PyDrawMesh_FaceNorm);
  m.def("draw_mesh_edge",      &PyDrawMesh_Edge);
}
