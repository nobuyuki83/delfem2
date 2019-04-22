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


class CMeshDynTri{
public:
  void Initialize(const double* aPo, int nPo, int ndim,
                  const int* aTri, int nTri)
  {
    aEPo.resize(nPo);
    for(int ipo=0;ipo<nPo;ipo++){
      if( ndim == 3 ){
        aEPo[ipo].p.x = aPo[ipo*3+0];
        aEPo[ipo].p.y = aPo[ipo*3+1];
        aEPo[ipo].p.z = aPo[ipo*3+2];
      }
      else if( ndim == 2 ){
        aEPo[ipo].p.x = aPo[ipo*2+0];
        aEPo[ipo].p.y = aPo[ipo*2+1];
        aEPo[ipo].p.z = 0.0;
      }
    }
    for(int itri=0;itri<nTri;itri++){
      unsigned int i1 = aTri[itri*3+0];
      unsigned int i2 = aTri[itri*3+1];
      unsigned int i3 = aTri[itri*3+2];
      aEPo[i1].e = itri; aEPo[i1].d = 0;
      aEPo[i2].e = itri; aEPo[i2].d = 1;
      aEPo[i3].e = itri; aEPo[i3].d = 2;
    }
    /////////////////
    aETri.resize(nTri);
    for(int itri=0;itri<nTri;itri++){
      aETri[itri].v[0] = aTri[itri*3+0];
      aETri[itri].v[1] = aTri[itri*3+1];
      aETri[itri].v[2] = aTri[itri*3+2];
    }
    {
      std::vector<int> elsup_ind, elsup;
      JArray_MakeElSuP(elsup_ind, elsup,
                      aETri, (int)aEPo.size());
      MakeInnerRelationTri(aETri, (int)aEPo.size(),
                           elsup_ind, elsup);
    }
  }
  void Check()
  {
    CheckTri(aETri);
    CheckTri(aEPo, aETri);
  }
  std::vector<double> MinMax_XYZ() const {
    double x_min,x_max, y_min,y_max, z_min,z_max;
    x_min=x_max=aEPo[0].p.x;
    y_min=y_max=aEPo[0].p.y;
    z_min=z_max=aEPo[0].p.z;
    for(unsigned int ipo=0;ipo<aEPo.size();ipo++){
      updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                      aEPo[ipo].p.x, aEPo[ipo].p.y, aEPo[ipo].p.z);
    }
    return {x_min,x_max, y_min,y_max, z_min,z_max};
  }
  int insertPointElem(int itri0, double r0, double r1){
    const int ipo0 = aEPo.size();
    CEPo2 p0;
    {
      int i0 = aETri[itri0].v[0];
      int i1 = aETri[itri0].v[1];
      int i2 = aETri[itri0].v[2];
      p0.p = r0*aEPo[i0].p+r1*aEPo[i1].p+(1-r0-r1)*aEPo[i2].p;
    }
    aEPo.push_back(p0);
    InsertPoint_Elem(ipo0, itri0, aEPo, aETri);
    return ipo0;
  }
  void DelaunayAroundPoint(int ipo){
    ::DelaunayAroundPoint(ipo, aEPo, aETri);
  }
  void Draw_FaceNorm()const { DrawMeshDynTri_FaceNorm(aEPo,aETri); }
  void Draw_Edge() const { DrawMeshDynTri_Edge(aEPo,aETri); }
  void draw() const { this->Draw_Edge(); }
  int nTri() const { return aETri.size(); }
  void DeleteTriEdge(int itri, int iedge){ Collapse_ElemEdge(itri, iedge, aEPo, aETri); }
public:
  std::vector<CEPo2> aEPo;
  std::vector<ETri> aETri;
};

void PyMeshDynTri3D_Initialize
(CMeshDynTri& mesh,
 const py::array_t<double>& po,
 const py::array_t<int>& tri)
{
  mesh.Initialize(po.data(), po.shape()[0], po.shape()[1],
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
  return std::forward_as_tuple(aPos,aElm, aPtrVtxInd,aVtxInd);
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
  
  py::class_<CMeshDynTri>(m, "CppMeshDynTri")
  .def(py::init<>())
  .def("draw", &CMeshDynTri::draw)
  .def("draw_face", &CMeshDynTri::Draw_FaceNorm)
  .def("draw_edge", &CMeshDynTri::Draw_Edge)
  .def("check", &CMeshDynTri::Check)
  .def("ntri",  &CMeshDynTri::nTri)
  .def("delete_tri_edge", &CMeshDynTri::DeleteTriEdge)
  .def("minmax_xyz",      &CMeshDynTri::MinMax_XYZ)
  .def("insert_point_elem", &CMeshDynTri::insertPointElem)
  .def("delaunay_around_point", &CMeshDynTri::DelaunayAroundPoint);
  
  m.def("meshdyntri3d_initialize",&PyMeshDynTri3D_Initialize);
  
  m.def("meshtri3d_read_ply",     &PyMeshTri3D_ReadPly,     py::return_value_policy::move);
  m.def("meshtri3d_read_obj",     &PyMeshTri3D_ReadObj,     py::return_value_policy::move);
  m.def("meshtri3d_read_nastran", &PyMeshTri3D_ReadNastran, py::return_value_policy::move);
  m.def("meshtri3d_write_obj",    &PyMeshTri3D_WriteObj);
  m.def("meshquad3d_subdiv",      &PyMeshQuad3D_Subviv,     py::return_value_policy::move);
  m.def("meshhex3d_subdiv",       &PyMeshHex3D_Subviv,      py::return_value_policy::move);
  m.def("meshquad2d_grid",        &PyMeshQuad2D_Grid,       py::return_value_policy::move);
  
  m.def("triangulation",&Triangulation);
  m.def("triangulation",&Triangulation,
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
