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
#include "delfem2/sdf.h"
#include "delfem2/isosurface_stuffing.h"
#include "delfem2/mathexpeval.h"

namespace py = pybind11;

//////////////////////////////////////////////////////////////////////////////////////////


void init_polyline(py::module &m);
void init_mshtopoio_gl(py::module &m);
void init_sampler(py::module &m);
void init_fbx(py::module &m);
void init_texture(py::module &m);
void init_rigidbody(py::module &m);
void init_field(py::module &m);
void init_fem(py::module &m);

std::tuple<std::vector<double>,std::vector<int>> PyMeshQuad3D_VoxelGrid
(const CGrid3D& vg)
{
  std::vector<double> aXYZ;
  std::vector<int> aQuad;
  vg.GetQuad(aXYZ, aQuad);
  return std::forward_as_tuple(aXYZ,aQuad);
}

std::tuple<std::vector<double>,std::vector<int>> PyMeshHex3D_VoxelGrid
(const CGrid3D& vg)
{
  std::vector<double> aXYZ;
  std::vector<int> aHex;
  vg.GetHex(aXYZ, aHex);
  return std::forward_as_tuple(aXYZ,aHex);
}

std::tuple<CMeshDynTri2D, py::array_t<int>, py::array_t<int>> MeshDynTri2D_Cad2D
(const CCad2D& cad, double len)
{
  CMeshDynTri2D dmesh;
  std::vector<int> aFlgPnt;
  std::vector<int> aFlgTri;
  cad.Meshing(dmesh, aFlgPnt,aFlgTri,
              len);
  assert( aFlgPnt.size() == dmesh.aVec2.size() );
  assert( aFlgTri.size() == dmesh.aETri.size() );
  ////
  py::array_t<int> npFlgPnt(aFlgPnt.size(), aFlgPnt.data());
  py::array_t<int> npFlgTri(aFlgTri.size(), aFlgTri.data());
  return std::forward_as_tuple(dmesh,npFlgPnt,npFlgTri);
}


std::tuple<py::array_t<double>, py::array_t<int>>
NumpyXYTri_MeshDynTri2D
(CMeshDynTri2D& dmesh)
{
  std::vector<double> aXY;
  std::vector<int> aTri;
  dmesh.Export_StlVectors(aXY,aTri);
  py::array_t<double> npXY({(int)aXY.size()/2,2}, aXY.data());
  py::array_t<int> npTri({(int)aTri.size()/3,3}, aTri.data());
  return std::tie(npXY,npTri);
}

py::array_t<int> PyCad2D_GetPointsEdge
(const CCad2D& cad,
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



std::tuple<py::array_t<double>, py::array_t<int>> PyIsoSurface
(const std::vector<const CSDF3*>& apSDF)
{
  class CMyInput : public CInputIsosurfaceStuffing
  {
  public:
    CMyInput(const std::vector<const CSDF3*>& apSDF){
      this->apSDF = apSDF;
    }
    virtual double SignedDistance(double px, double py, double pz) const{
      double n[3];
      double max_dist = apSDF[0]->Projection(px,py,pz, n);
      for(unsigned int ipct=1;ipct<apSDF.size();ipct++){
        double dist0,n0[3];
        dist0 = apSDF[ipct]->Projection(px,py,pz, n0);
        if( dist0 < max_dist ) continue;
        max_dist = dist0;
        n[0] = n0[0];
        n[1] = n0[1];
        n[2] = n0[2];
      }
      return max_dist;
    }
    virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                       double px, double py, double pz) const
    {
      sdf = this->SignedDistance(px,py,pz);
      /*
      const double rad0 = sp.radius_;
      const double rad1 = sqrt(px*px+py*py+pz*pz);
      if( rad1 > rad0*0.5 ){ ilevel_vol = 0; }
      else{ ilevel_vol = 1; }
       */
      ilevel_srf = 1;
      nlayer = 1;
      ilevel_vol = -1;
    }
  public:
    std::vector<const CSDF3*> apSDF;
  } input(apSDF);
  
  std::vector<double> aXYZ;
  std::vector<int> aTet;
  std::vector<int> aIsOnSurfXYZ;
  double rad = 1.5;
  double cent[3] = {0,0,0};
  IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                     input, 0.2, rad*4.0, cent);
  
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<int> npTet({(int)aTet.size()/4,4}, aTet.data());
  return std::tie(npXYZ,npTet);
}

py::array_t<double> PyMVC
(const py::array_t<double>& XY,
 const py::array_t<double>& XY_bound)
{
  assert(XY.ndim()==2);
  assert(XY.shape()[1]==2);
  assert(XY_bound.ndim()==2);
  assert(XY_bound.shape()[1]==2);
  const int np = XY.shape()[0];
  const int npb = XY_bound.shape()[0];
  py::array_t<double> aW({np,npb});
  auto buff_w = aW.request();
  for(int ip=0;ip<np;++ip){
    MeanValueCoordinate2D((double*)buff_w.ptr+ip*npb,
                          XY.at(ip,0), XY.at(ip,1),
                          XY_bound.data(), npb);
  }
  return aW;
}



PYBIND11_MODULE(libdelfem2, m) {
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
  init_fem(m);
  
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
  py::class_<CGrid3D>(m, "CppVoxelGrid", "voxel grid class")
  .def(py::init<>())
  .def("add",&CGrid3D::Add,"add voxel at the integer coordinate");
  
  m.def("meshquad3d_voxelgrid",&PyMeshQuad3D_VoxelGrid);
  m.def("meshhex3d_voxelgrid",&PyMeshHex3D_VoxelGrid);
  
  ///////////////////////////////////
  // SDF
  py::class_<CSDF3>(m, "SDF");
  
  py::class_<CSignedDistanceField3D_Sphere, CSDF3>(m, "SDF_Sphere")
  .def(py::init<>())
  .def(py::init<double,const std::vector<double>&,bool>())
  .def("draw",  &CSignedDistanceField3D_Sphere::Draw);
  
  m.def("isosurface", &PyIsoSurface);
  
  ///////////////////////////////////
  // cad
  py::class_<CCad2D>(m, "CppCad2D", "2D CAD class")
  .def(py::init<>())
  .def("draw",        &CCad2D::Draw)
  .def("pick",        &CCad2D::Pick)
  .def("drag_picked", &CCad2D::DragPicked)
  .def("minmax_xyz",  &CCad2D::MinMaxXYZ)
  .def("add_polygon", &CCad2D::AddPolygon)
  .def("meshing",     &CCad2D::Meshing)
  .def("xy_vtx_face", &CCad2D::XY_Vtx_Face)
  .def("ind_vtx_face", &CCad2D::Ind_Vtx_Face)
  .def("ind_edge_face",&CCad2D::Ind_Edge_Face)
  .def("add_vtx_edge", &CCad2D::AddVtxEdge)
  .def("check",       &CCad2D::Check)
  .def("nface",       &CCad2D::nFace)
  .def("nvtx",        &CCad2D::nVtx)
  .def("nedge",       &CCad2D::nEdge)
  .def_readwrite("is_draw_face", &CCad2D::is_draw_face)
  .def_readwrite("ivtx_picked",  &CCad2D::ivtx_picked)
  .def_readwrite("iedge_picked",  &CCad2D::iedge_picked)
  .def_readwrite("iface_picked",  &CCad2D::iface_picked);
  

  m.def("cad_getPointsEdge",
        &PyCad2D_GetPointsEdge,
        py::arg("cad"),
        py::arg("list_edge_index"),
        py::arg("np_xy"),
        py::arg("tolerance") = 0.001,
        py::return_value_policy::move);
  
  m.def("meshDynTri2D_CppCad2D",&MeshDynTri2D_Cad2D);
  m.def("numpyXYTri_MeshDynTri2D",&NumpyXYTri_MeshDynTri2D);
  
  ////////////////////////////////////

  py::class_<CColorMap>(m,"ColorMap")
  .def(py::init<>())
  .def(py::init<double, double, const std::string&>());
  
  py::class_<CMathExpressionEvaluator>(m,"MathExpressionEvaluator")
  .def(py::init<>())
  .def("set_expression",&CMathExpressionEvaluator::SetExp)
  .def("set_key",       &CMathExpressionEvaluator::SetKey)
  .def("eval",          &CMathExpressionEvaluator::Eval);
  
   ///////////////////////////////////
  // gl misc
  m.def("setSomeLighting",  &setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl",       &setUpGLSL, "compile shader program");
  m.def("glew_init",        &glewInit);
  m.def("draw_sphere",      &DrawSphereAt );
  m.def("mvc",              &PyMVC);
}


