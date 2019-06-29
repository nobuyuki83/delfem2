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

void init_sampler(py::module &m);
void init_texture(py::module &m);

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


std::tuple<py::array_t<double>, py::array_t<unsigned int>>
NumpyXYTri_MeshDynTri2D
(CMeshDynTri2D& dmesh)
{
  std::vector<double> aXY;
  std::vector<unsigned int> aTri;
  dmesh.Export_StlVectors(aXY,aTri);
  py::array_t<double> npXY({(int)aXY.size()/2,2}, aXY.data());
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
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



std::tuple<py::array_t<double>, py::array_t<unsigned int>> PyIsoSurface
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
  std::vector<unsigned int> aTet;
  std::vector<int> aIsOnSurfXYZ;
  double rad = 1.5;
  double cent[3] = {0,0,0};
  IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                     input, 0.2, rad*4.0, cent);
  
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<unsigned int> npTet({(int)aTet.size()/4,4}, aTet.data());
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



PYBIND11_MODULE(c_gl, m) {
  m.doc() = "pybind11 delfem2 binding";
  ///////////////////////////////////
  init_sampler(m);
  init_texture(m);
   
 ////////////////////////////////////

  py::class_<CColorMap>(m,"ColorMap")
  .def(py::init<>())
  .def(py::init<double, double, const std::string&>());
  
  ///////////////////////////////////
  // gl misc
  m.def("setSomeLighting",  &setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl",       &setUpGLSL, "compile shader program");
  m.def("glew_init",        &glewInit);
  m.def("draw_sphere",      &DrawSphereAt );
}


