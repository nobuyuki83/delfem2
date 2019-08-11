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
#include "delfem2/primitive.h"
#include "delfem2/isosurface_stuffing.h"

namespace py = pybind11;

//////////////////////////////////////////////////////////////////////////////////////////


std::tuple<py::array_t<double>, py::array_t<unsigned int>> PyIsoSurface
(const std::vector<const CSDF3*>& apSDF)
{
  class CMyInput : public CInput_IsosurfaceStuffing
  {
  public:
    CMyInput(const std::vector<const CSDF3*>& apSDF){
      this->apSDF = apSDF;
    }
    virtual double SignedDistance(double px, double py, double pz) const{
      double n[3];
      double max_dist = apSDF[0]->Projection(n,
                                             px,py,pz);
      for(unsigned int ipct=1;ipct<apSDF.size();ipct++){
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



void init_sdf(py::module &m){
  
  ///////////////////////////////////
  // SDF
  py::class_<CSDF3>(m, "CppSDF3");
  py::class_<CSphere, CSDF3>(m, "CppSDF3_Sphere")
  .def(py::init<>())
  .def(py::init<double,const std::vector<double>&,bool>())
  .def_readwrite("cent", &CSphere::cent_)
  .def_readwrite("rad",  &CSphere::radius_);
  
  m.def("isosurface", &PyIsoSurface);  
}


