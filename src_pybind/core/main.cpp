/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <map>
#include <deque>

#include "../py_funcs.h"

#include "delfem2/mat3.h"
#include "delfem2/primitive.h"
#include "delfem2/mshtopoio.h"
#include "delfem2/gridvoxel.h"
#include "delfem2/gridcube.h"
#include "delfem2/bv.h"
#include "delfem2/iss.h"
#include "delfem2/slice.h"
#include "delfem2/evalmathexp.h"

#include "delfem2/cad2_dtri2.h"
#include "delfem2/rig_geo3.h"

#include "tinygltf/tiny_gltf.h"
#include "io_gltf.h"

#include "stb_image.h" // stb is already compiled in io_gltf.cpp

namespace py = pybind11;
namespace dfm2 = delfem2;

// ----------------------------------------------------------

//void init_rigidbody(py::module &m);
void init_polyline(py::module &m);
void init_mshtopoio(py::module &m);
void init_field(py::module &m);
void init_fem(py::module &m);
void init_sdf(py::module &m);

// ---------------------------------
// img related

py::array_t<unsigned char> PyImRead(const std::string& d)
{
  int width, height, channels;
  unsigned char *img = stbi_load(d.c_str(),
                                 &width, &height, &channels, 0);
  py::array_t<unsigned char> npR({height,width,channels}, img);
  delete[] img;
  return npR;
}

// ---------------------------------
// voxel related

std::tuple<std::vector<double>,std::vector<unsigned int>> PyMeshQuad3D_VoxelGrid
(const dfm2::CGrid3<int>& vg)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aQuad;
  dfm2::MeshQuad3D_VoxelGrid(aXYZ,aQuad,
                             vg.ndivx, vg.ndivy, vg.ndivz,
                             vg.aVal);
  return std::make_tuple(aXYZ,aQuad);
}

std::tuple<std::vector<double>,std::vector<int>> PyMeshHex3D_VoxelGrid
(const dfm2::CGrid3<int>& vg)
{
  std::vector<double> aXYZ;
  std::vector<int> aHex;
  dfm2::MeshHex3D_VoxelGrid(aXYZ,aHex,
                            vg.ndivx, vg.ndivy, vg.ndivz,
                            vg.aVal);
  return std::make_tuple(aXYZ,aHex);
}

// -------------------------------------------

std::tuple<py::array_t<double>, py::array_t<unsigned int>>
NumpyXYTri_MeshDynTri2D
(dfm2::CMeshDynTri2D& dmesh)
{
  std::vector<double> aXY;
  std::vector<unsigned int> aTri;
  dmesh.Export_StlVectors(aXY,aTri);
  py::array_t<double> npXY({(int)aXY.size()/2,2}, aXY.data());
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  return std::make_tuple(npXY,npTri);
}

py::array_t<int> PyCad2D_GetPointsEdge
(const dfm2::CCad2D& cad,
 const std::vector<int>& aIE,
 const py::array_t<double>& aXY,
 double torelance)
{
  std::vector<int> aIdP;
  cad.GetPointsEdge(aIdP,
                    aXY.data(), aXY.shape()[0],
                    aIE,torelance);
  std::set<int> setIdP(aIdP.begin(),aIdP.end());
  aIdP.assign(setIdP.begin(),setIdP.end());
  return py::array_t<int>((int)aIdP.size(), aIdP.data());
}

py::array_t<double> PyMVC
(const py::array_t<double>& XY,
 const py::array_t<double>& XY_bound)
{
  assert( AssertNumpyArray2D(XY, -1, 2) );
  assert( AssertNumpyArray2D(XY_bound, -1, 2) );
  const int np = XY.shape()[0];
  const int npb = XY_bound.shape()[0];
  py::array_t<double> aW({np,npb});
  auto buff_w = aW.request();
  for(int ip=0;ip<np;++ip){
    dfm2::MeanValueCoordinate2D((double*)buff_w.ptr+ip*npb,
                                XY.at(ip,0), XY.at(ip,1),
                                XY_bound.data(), npb);
  }
  return aW;
}

py::array_t<double> PyRotMat3_Cartesian(const std::vector<double>& d)
{
  dfm2::CMat3d m;
  m.SetRotMatrix_Cartesian(d[0], d[1], d[2]);
  py::array_t<double> npR({3,3});
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      npR.mutable_at(i,j) = m.Get(i, j);
    }
  }
  return npR;
}

// -----------------------------------------------
// Rigging related from here

std::tuple<py::array_t<double>, py::array_t<unsigned int>, py::array_t<double>, py::array_t<unsigned int>>
PyGLTF_GetMeshInfo
(const dfm2::CGLTF& gltf,
 int imesh, int iprimitive)
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<double> aRigWeight;
  std::vector<unsigned int> aRigJoint;
  gltf.GetMeshInfo(aXYZ0,aTri,aRigWeight,aRigJoint,
                   imesh, iprimitive);
  const int np = aXYZ0.size()/3;
  assert( (int)aRigWeight.size() == np*4 );
  assert( (int)aRigJoint.size() == np*4 );
  py::array_t<double> npXYZ0({np,3}, aXYZ0.data());
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npRW({np,4}, aRigWeight.data());
  py::array_t<unsigned int> npRJ({np,4}, aRigJoint.data());
  return std::make_tuple(npXYZ0,npTri,npRW,npRJ);
}

class CBoneArray{
public:
  void SetTranslation(int ib, const std::vector<double>& aT){
    assert(aT.size()==3);
    aRigBone[ib].SetTranslation(aT[0], aT[1], aT[2]);
    UpdateBoneRotTrans(aRigBone);
  }
  void SetRotationBryant(int ib, const std::vector<double>& aRB){
    assert(aRB.size()==3);
    aRigBone[ib].SetRotationBryant(aRB[0], aRB[1], aRB[2]);
    UpdateBoneRotTrans(aRigBone);
  }
public:
  std::vector<dfm2::CRigBone> aRigBone;
};

CBoneArray
PyGLTF_GetBones
(const dfm2::CGLTF& gltf,
 int iskin)
{
  CBoneArray BA;
  gltf.GetBone(BA.aRigBone,
               iskin);
  return BA;
}

void PyUpdateRigSkin
(py::array_t<double>& npXYZ,
 const py::array_t<double>& npXYZ0,
 const py::array_t<unsigned int>& npTri,
 const CBoneArray& BA,
 const py::array_t<double>& npRigWeight,
 const py::array_t<unsigned int>& npRigJoint)
{
  assert( AssertNumpyArray2D(npXYZ, -1, 3) );
  assert( AssertNumpyArray2D(npXYZ0, -1, 3) );
  assert( AssertNumpyArray2D(npTri, -1, 3) );
  assert( AssertNumpyArray2D(npRigWeight, -1, 4) );
  assert( AssertNumpyArray2D(npRigJoint, -1, 4) );
  assert( npXYZ.shape()[0] == npXYZ0.shape()[0] );
  assert( npXYZ.shape()[0] == npRigWeight.shape()[0] );
  assert( npXYZ.shape()[0] == npRigJoint.shape()[0] );
  double* aXYZ = (double*)(npXYZ.request().ptr);
  dfm2::Skinning_LBS_LocalWeight(aXYZ,
                                 npXYZ0.data(), npXYZ0.shape()[0],
                                 npTri.data(), npTri.shape()[0],
                                 BA.aRigBone,
                                 npRigWeight.data(),
                                 npRigJoint.data());
}

// Rigging related ends here
// -----------------------------------------

void PyCad2D_ImportSVG
 (dfm2::CCad2D& cad,
  const std::string& path_svg,
  double scale_x,
  double scale_y)
{
  std::vector< std::vector<dfm2::CCad2D_EdgeGeo> > aaEdge;
  ReadSVG_LoopEdgeCCad2D(aaEdge,
                         path_svg);
  cad.Clear();
  for(unsigned int ie=0;ie<aaEdge.size();++ie){
    std::vector<dfm2::CCad2D_EdgeGeo> aEdge = aaEdge[ie];
    Transform_LoopEdgeCad2D(aEdge,false,true,scale_x,scale_y);
    if( AreaLoop(aEdge) < 0 ){ aEdge = InvertLoop(aEdge); }
    aEdge = RemoveEdgeWithZeroLength(aEdge);
    for(unsigned int ie=0;ie<aEdge.size();++ie){ aEdge[ie].GenMeshLength(-1); }
    cad.AddFace(aEdge);
  }
}


void PyIsoSurfaceToSVG
 (const py::array_t<double>& npXY,
  const py::array_t<unsigned int>& npTri,
  const py::array_t<double>& npVal,
  double scale)
{
  assert( AssertNumpyArray2D(npXY, -1, 2) );
  assert( AssertNumpyArray2D(npTri, -1, 3) );
  assert( npVal.ndim() == 1 );
  assert( npVal.size() == npXY.shape()[0] );
  assert( npVal.strides()[0] == sizeof(double) );
  std::vector<dfm2::CSegInfo> aSeg;
  dfm2::AddContour(aSeg,
                   0.0,
                   npTri.data(), npTri.shape()[0],
                   npVal.data());
  std::vector<double> aXY_Line(aSeg.size()*4);
  for(unsigned int iseg=0;iseg<aSeg.size();++iseg){
    double pA[2], pB[2];
    aSeg[iseg].Pos2D(pA, pB,
                     npXY.data(), npTri.data());
    aXY_Line[iseg*4+0] = pA[0]*scale;
    aXY_Line[iseg*4+1] = pA[1]*scale;
    aXY_Line[iseg*4+2] = pB[0]*scale;
    aXY_Line[iseg*4+3] = pB[1]*scale;
  }
  std::ofstream fout("hoge.svg");
  fout << "<?xml version=\"1.0\"?>" << std::endl;;
  fout << "<svg xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;;
  for(unsigned int il=0;il<aXY_Line.size()/4;++il){
    fout << "<line";
    fout << " x1=\"" << aXY_Line[il*4+0] << "\" y1=\"" << -aXY_Line[il*4+1] << "\"";
    fout << " x2=\"" << aXY_Line[il*4+2] << "\" y2=\"" << -aXY_Line[il*4+3] << "\"";
    fout << " stroke=\"black\" stroke-width=\"2\" />" << std::endl;
  }
  fout << "</svg>" << std::endl;
}


PYBIND11_MODULE(c_core, m) {
  m.doc() = "pybind11 delfem2 binding";
  // ------------------------
  
  init_mshtopoio(m);
  init_polyline(m);
  init_field(m);
  init_fem(m);
  init_sdf(m);
//  init_rigidbody(m);
  
  // ----------------------
  m.def("imread", &PyImRead);
  
  // -------------------------
  // axis arrigned boudning box
  py::class_<dfm2::CBV3d_AABB>(m,"AABB3", "3D axis aligned bounding box class")
  .def(py::init<>())
  .def(py::init<const std::vector<double>&, const std::vector<double>&>())
  .def(py::init<const std::vector<double>&>())
  .def("__str__",            &dfm2::CBV3d_AABB::str, "print x_min,x_max,y_min,y_max,z_min,z_max")
  .def("minmax_xyz",         &dfm2::CBV3d_AABB::AABBVec3)
  // 
  .def("set_minmax_xyz",     &dfm2::CBV3d_AABB::Set_AABBVec3)
  .def("add_minmax_xyz",     &dfm2::CBV3d_AABB::Add_AABBVec3)
  .def("list_xyz",           &dfm2::CBV3d_AABB::Point3D_Vox, "corner xyz coords in voxel point order")
  .def("diagonal_length",    &dfm2::CBV3d_AABB::DiagonalLength, "diagonal length of the bounding box")
  .def("max_length",         &dfm2::CBV3d_AABB::MaxLength, "diagonal length of the bounding box")
  .def("center",             &dfm2::CBV3d_AABB::Center, "center position")
  .def("is_active",          &dfm2::CBV3d_AABB::IsActive);
  
  // --------
  // voxel
  py::class_<dfm2::CGrid3<int>>(m, "CppVoxelGrid", "voxel grid class")
  .def(py::init<>())
  .def("add",&dfm2::CGrid3<int>::Set,"add voxel at the integer coordinate")
  .def("initialize",&dfm2::CGrid3<int>::Initialize,"add voxel at the integer coordinate");
  
  m.def("meshquad3d_voxelgrid",&PyMeshQuad3D_VoxelGrid);
  m.def("meshhex3d_voxelgrid", &PyMeshHex3D_VoxelGrid);
  
  // -------------
  // cad
  py::class_<dfm2::CCad2D>(m, "CppCad2D", "2D CAD class")
  .def(py::init<>())
  .def("minmax_xyz",     &dfm2::CCad2D::MinMaxXYZ)
  //
  .def("clear",          &dfm2::CCad2D::Clear)
  .def("pick",           &dfm2::CCad2D::Pick)
  .def("drag_picked",    &dfm2::CCad2D::DragPicked)
  .def("add_polygon",    &dfm2::CCad2D::AddPolygon)
  .def("add_vtx_edge",   &dfm2::CCad2D::AddVtxEdge)
  .def("add_vtx_face",   &dfm2::CCad2D::AddVtxFace)
  .def("xy_vtxctrl_face",&dfm2::CCad2D::XY_VtxCtrl_Face)
  .def("xy_vtx",         &dfm2::CCad2D::XY_Vtx)
  .def("ind_vtx_face",   &dfm2::CCad2D::Ind_Vtx_Face)
  .def("ind_edge_face",  &dfm2::CCad2D::Ind_Edge_Face)
  .def("ind_vtx_edge",   &dfm2::CCad2D::Ind_Vtx_Edge)
  .def("set_edge_type",  &dfm2::CCad2D::SetEdgeType)
  .def("edge_type",      &dfm2::CCad2D::GetEdgeType)
  .def("check",          &dfm2::CCad2D::Check)
  .def("nface",          &dfm2::CCad2D::nFace)
  .def("nvtx",           &dfm2::CCad2D::nVtx)
  .def("nedge",          &dfm2::CCad2D::nEdge)
  .def_readwrite("is_draw_face",  &dfm2::CCad2D::is_draw_face)
  .def_readwrite("ivtx_picked",   &dfm2::CCad2D::ivtx_picked)
  .def_readwrite("iedge_picked",  &dfm2::CCad2D::iedge_picked)
  .def_readwrite("iface_picked",  &dfm2::CCad2D::iface_picked);
  
  py::enum_<dfm2::CCad2D_EdgeGeo::EDGE_TYPE>(m, "CAD_EDGE_GEOM_TYPE")
  .value("CAD_EDGE_GEOM_LINE",             dfm2::CCad2D_EdgeGeo::LINE)
  .value("CAD_EDGE_GEOM_BEZIER_CUBIC",     dfm2::CCad2D_EdgeGeo::BEZIER_CUBIC)
  .value("CAD_EDGE_GEOM_BEZIER_QUADRATIC", dfm2::CCad2D_EdgeGeo::BEZIER_QUADRATIC)
  .export_values();
  
  m.def("cppCad2D_ImportSVG",
        &PyCad2D_ImportSVG);
  m.def("cppSVG_Polyline",
        &dfm2::Str_SVGPolygon);  
  
  py::class_<dfm2::CMesher_Cad2D>(m,"CppMesher_Cad2D")
  .def(py::init<>())
  .def("meshing",               &dfm2::CMesher_Cad2D::Meshing)
  .def("points_on_one_edge",    &dfm2::CMesher_Cad2D::IndPoint_IndEdge)
  .def("points_on_edges",       &dfm2::CMesher_Cad2D::IndPoint_IndEdgeArray)
  .def("points_on_faces",       &dfm2::CMesher_Cad2D::IndPoint_IndFaceArray)
  .def_readwrite("edge_length", &dfm2::CMesher_Cad2D::edge_length);

  m.def("cad_getPointsEdge",
        &PyCad2D_GetPointsEdge,
        py::arg("cad"),
        py::arg("list_edge_index"),
        py::arg("np_xy"),
        py::arg("tolerance") = 0.001,
        py::return_value_policy::move);
  
  m.def("numpyXYTri_MeshDynTri2D",&NumpyXYTri_MeshDynTri2D);
  
  py::class_<dfm2::CGLTF>(m,"CppGLTF")
  .def(py::init<>())
  .def("read", &dfm2::CGLTF::Read)
  .def("print", &dfm2::CGLTF::Print);
  
  py::class_<CBoneArray>(m,"CppBoneArray")
  .def("set_translation", &CBoneArray::SetTranslation)
  .def("set_rotation_bryant", &CBoneArray::SetRotationBryant)
  .def(py::init<>());
  
  m.def("CppGLTF_GetMeshInfo",   &PyGLTF_GetMeshInfo);
  m.def("CppGLTF_GetBones",      &PyGLTF_GetBones);
  m.def("update_rig_skin",       &PyUpdateRigSkin);
  m.def("update_bone_transform", &dfm2::UpdateBoneRotTrans);

  // ------------------------------------------

  py::class_<dfm2::CMathExpressionEvaluator>(m,"MathExpressionEvaluator")
  .def(py::init<>())
  .def("set_expression",&dfm2::CMathExpressionEvaluator::SetExp)
  .def("set_key",       &dfm2::CMathExpressionEvaluator::SetKey)
  .def("eval",          &dfm2::CMathExpressionEvaluator::Eval);
  
  m.def("mvc",              &PyMVC);
  m.def("rotmat3_cartesian", &PyRotMat3_Cartesian);
  
  m.def("isoline_svg", &PyIsoSurfaceToSVG);
}


