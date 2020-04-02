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

#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/primitive.h"
#include "delfem2/iss.h"

#include "../py_funcs.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

// ------------------------------------------------------------------------

std::tuple<py::array_t<double>, py::array_t<unsigned int>> PyIsoSurface
(const std::vector<const dfm2::CSDF3*>& apSDF)
{
  class CMyInput : public delfem2::CInput_IsosurfaceStuffing
  {
  public:
    CMyInput(const std::vector<const dfm2::CSDF3*>& apSDF){
      this->apSDF = apSDF;
    }
    virtual double SignedDistance(double px, double py, double pz) const{
      double n[3];
      double max_dist = apSDF[0]->Projection(n,
                                             px,py,pz);
      for(std::size_t ipct=1;ipct<apSDF.size();ipct++){
        double dist0,n0[3];
        dist0 = apSDF[ipct]->Projection(n0,
                                        px,py,pz);
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
    std::vector<const dfm2::CSDF3*> apSDF;
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


void PyProjectPointOutsideSDF
(py::array_t<double>& npXYZt,
 const std::vector<const dfm2::CSDF3*>& apSDF)
{
  assert( AssertNumpyArray2D(npXYZt, -1, 3) );
  const int np = npXYZt.shape()[0];
  double* pXYZt = (double*)(npXYZt.request().ptr);
  for(int ip=0;ip<np;++ip){
    const double px = npXYZt.at(ip,0);
    const double py = npXYZt.at(ip,1);
    const double pz = npXYZt.at(ip,2);
    double n[3];
    double max_dist = apSDF[0]->Projection(n,
                                           px,py,pz);
    for(unsigned int ipct=1;ipct<apSDF.size();ipct++){
      double dist0,n0[3];
      dist0 = apSDF[ipct]->Projection(n0,
                                      px,py,pz);
      if( dist0 < max_dist ){ continue; }
      max_dist = dist0;
      n[0] = n0[0];
      n[1] = n0[1];
      n[2] = n0[2];
    }
    if( max_dist <= 0 ){ continue; }
    pXYZt[ip*3+0] += n[0]*max_dist;
    pXYZt[ip*3+1] += n[1]*max_dist;
    pXYZt[ip*3+2] += n[2]*max_dist;
  }
}

class CPyCollision_Points_MeshTri3D
{
public:
  void SetMesh(const py::array_t<double>& xyz,
               const py::array_t<unsigned int>& tri,
               double margin)
  {
    contact_clearance = margin;
    bvh.Init(xyz.data(), xyz.shape()[0],
             tri.data(), tri.shape()[0],
             margin);
  }
  void Project(py::array_t<double>& npXYZt,
               const py::array_t<double>& npXYZ,
               const py::array_t<unsigned int>& npTri,
               const py::array_t<double>& npNorm,
               double rad_explore)
  {
    assert( AssertNumpyArray2D(npXYZ, -1, 3) );
    assert( AssertNumpyArray2D(npTri, -1, 3) );
    assert( AssertNumpyArray2D(npNorm, -1, 3) );
    assert( AssertNumpyArray2D(npXYZt, -1, 3) );
    assert( npXYZ.shape()[0] == npNorm.shape()[0] );
    assert( aInfoNearest.empty() || (int)aInfoNearest.size() == npXYZt.shape()[0] );
    double* pXYZt = (double*)(npXYZt.request().ptr);
    Project_PointsIncludedInBVH_Outside_Cache(pXYZt, aInfoNearest,
                                              npXYZt.shape()[0],
                                              contact_clearance,bvh,
                                              npXYZ.data(), npXYZ.shape()[0],
                                              npTri.data(), npTri.shape()[0],
                                              npNorm.data(), rad_explore);
  }
public:
  double contact_clearance;
  std::vector<dfm2::CInfoNearest<double>> aInfoNearest;
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere,double> bvh;
};

void init_sdf(py::module &m){
  
  // --------------------
  // SDF
  py::class_<dfm2::CSDF3>(m, "CppSDF3");
  
  py::class_<delfem2::CSphere<double>, dfm2::CSDF3>(m, "CppSDF3_Sphere")
  .def(py::init<>())
  .def(py::init<double,const std::vector<double>&,bool>())
  .def_readwrite("cent", &delfem2::CSphere<double>::cent_)
  .def_readwrite("rad",  &delfem2::CSphere<double>::radius_);
  
  py::class_<CPyCollision_Points_MeshTri3D>(m, "CppClliderPointsMeshTri3D")
  .def("set_mesh", &CPyCollision_Points_MeshTri3D::SetMesh)
  .def("project",  &CPyCollision_Points_MeshTri3D::Project)
  .def(py::init<>());
  
  m.def("project_points_outside_sdf", &PyProjectPointOutsideSDF);
  
  m.def("isosurface", &PyIsoSurface);  
}


