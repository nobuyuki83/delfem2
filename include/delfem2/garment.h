/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GARMENT_H
#define DFM2_GARMENT_H

#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/cad2d_v2dtri.h"

namespace delfem2 {

class CRigidTrans_2DTo3D
{
public:
  CRigidTrans_2DTo3D(){
    radinv_x = 0.0;
    R.SetIdentity();
  }
  delfem2::CVec3d Transform(const delfem2::CVec2d& pi) const {
    delfem2::CVec3d p0(pi.x()-org2.x(), pi.y()-org2.y(),0.0);
    delfem2::CVec3d p2 = p0;
    if( radinv_x < 1.0e-5 ) {
      double x0 = p0.x();
      p2.p[0] = x0;
      p2.p[2] = -0.5*radinv_x*x0*x0;
    }
    else{
      double x0 = p0.x();
      p2.p[0] = (1.0/radinv_x)*sin(radinv_x*x0);
      p2.p[2] = -(1.0/radinv_x)*(1-cos(radinv_x*x0));
    }
    delfem2::CVec3d p3;
    delfem2::MatVec3(p3.p, R.mat,p2.p);
    delfem2::CVec3d p4 = org3 + p3;
    return p4;
  }
public:
  delfem2::CVec2d org2;
  delfem2::CVec3d org3;
  delfem2::CMat3d R;
  double radinv_x;
};


void MeshingPattern
 (std::vector<CDynTri>& aETri,
  std::vector<CVec2d>& aVec2,
  std::vector<double>& aXYZ, // deformed vertex positions
  std::vector<unsigned int>& aLine,
  const std::vector<CRigidTrans_2DTo3D>& aRT23,
  const delfem2::CCad2D& cad,
  const std::vector<unsigned int>& aIESeam,
  const double mesher_edge_length)
{
  delfem2::CMesher_Cad2D mesher;
  mesher.edge_length = mesher_edge_length;
  const double el = mesher.edge_length;
  for(int ie=0;ie<aIESeam.size()/2;++ie){
    const unsigned int ie0 = aIESeam[ie*2+0];
    const unsigned int ie1 = aIESeam[ie*2+1];
    const double l0 = cad.aEdge[ie0].LengthMesh();
    const double l1 = cad.aEdge[ie1].LengthMesh();
    const unsigned int ndiv = (int)(0.5*(l0+l1)/el+1);
    mesher.mapIdEd_NDiv[ie0] = ndiv;
    mesher.mapIdEd_NDiv[ie1] = ndiv;
  }
  delfem2::CMeshDynTri2D dmesh;
  mesher.Meshing(dmesh,
                 cad);
  dmesh.Check();
  aETri = dmesh.aETri;
  aVec2 = dmesh.aVec2;
  // ---
  const int np = aVec2.size();
  aXYZ.resize(np*3);
  for(int ifc=0;ifc<cad.aFace.size();++ifc){
    const CRigidTrans_2DTo3D& rt23 = aRT23[ifc];
    std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,ifc), cad);
    for(int ip : aIP){
      rt23.Transform(aVec2[ip]).CopyTo(aXYZ.data()+ip*3);
    }
  }
  for(int ie=0;ie<aIESeam.size()/2;++ie){
    unsigned int ie0 = aIESeam[ie*2+0];
    unsigned int ie1 = aIESeam[ie*2+1];
    std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(ie0, true, cad);
    std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(ie1, true, cad);
    const int npe = aIP0.size();
    assert( aIP1.size() == npe );
    for(int iip=0;iip<npe;++iip){
      int ip0 = aIP0[iip];
      int ip1 = aIP1[npe-iip-1];
      aLine.push_back(ip0);
      aLine.push_back(ip1);
    }
  }
}

}



#endif /* garment_h */
