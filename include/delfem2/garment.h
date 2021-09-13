/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GARMENT_H
#define DFM2_GARMENT_H

#include "delfem2/cad2_dtri2.h"
#include "delfem2/pbd_geo3dtri23.h"
#include "delfem2/pbd_geo3.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/srchbv3aabb.h"
#include "delfem2/srchbv3sphere.h"

namespace delfem2 {

/**
 * @brief Rigid Transformation from 2D to 3D
 */
class CRigidTrans_2DTo3D
{
public:
  CRigidTrans_2DTo3D(){
    radinv_x = 0.0;
    R.SetIdentity();
  }
  delfem2::CVec3d Transform(const delfem2::CVec2d& pi) const {
    delfem2::CVec3d p0(pi.x-org2.x, pi.y-org2.y,0.0);
    delfem2::CVec3d p2 = p0;
    if( radinv_x < 1.0e-5 ) {
      double x0 = p0.x;
      p2.p[0] = x0;
      p2.p[2] = -0.5*radinv_x*x0*x0;
    }
    else{
      double x0 = p0.x;
      p2.p[0] = (1.0/radinv_x)*sin(radinv_x*x0);
      p2.p[2] = -(1.0/radinv_x)*(1-cos(radinv_x*x0));
    }
    delfem2::CVec3d p3;
    delfem2::MatVec3(p3.data(), R.data(),p2.data());
    delfem2::CVec3d p4 = org3 + p3;
    return p4;
  }
public:
  delfem2::CVec2d org2;
  delfem2::CVec3d org3;
  delfem2::CMat3d R;
  double radinv_x;
};


/**
 * @brief meshing clothing pattern such that the opposing sides of the cloth panel have the same number of discretization
 */
void MeshingPattern
 (std::vector<CDynTri>& aETri,
  std::vector<CVec2d>& aVec2,
  std::vector<double>& aXYZ, // deformed vertex positions
  std::vector<unsigned int>& aLine,
  CMesher_Cad2D& mesher,
  const std::vector<CRigidTrans_2DTo3D>& aRT23,
  const delfem2::CCad2D& cad,
  const std::vector<unsigned int>& aIESeam,
  const double mesher_edge_length)
{
  mesher.edge_length = mesher_edge_length;
  const double el = mesher.edge_length;
  for(unsigned int ie=0;ie<aIESeam.size()/2;++ie){
    const unsigned int ie0 = aIESeam[ie*2+0];
    const unsigned int ie1 = aIESeam[ie*2+1];
    if( ie0 >= cad.aEdge.size() || ie1 >= cad.aEdge.size() ) continue;
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
  const unsigned int np = aVec2.size();
  aXYZ.resize(np*3);
  for(unsigned int ifc=0;ifc<cad.aFace.size();++ifc){
    const CRigidTrans_2DTo3D& rt23 = aRT23[ifc];
    std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,ifc), cad);
    for(int ip : aIP){
      rt23.Transform(aVec2[ip]).CopyTo(aXYZ.data()+ip*3);
    }
  }
  for(unsigned int ie=0;ie<aIESeam.size()/2;++ie){
    unsigned int ie0 = aIESeam[ie*2+0];
    unsigned int ie1 = aIESeam[ie*2+1];
    std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(ie0, true, cad);
    std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(ie1, true, cad);
    const unsigned int npe = aIP0.size();
    assert( aIP1.size() == npe );
    for(unsigned int iip=0;iip<npe;++iip){
      int ip0 = aIP0[iip];
      int ip1 = aIP1[npe-iip-1];
      aLine.push_back(ip0);
      aLine.push_back(ip1);
    }
  }
}

class CProjectorMesh{
public:
  void Init(){
    aNorm_Body.resize(aXYZ_Body.size());
    delfem2::Normal_MeshTri3D(
        aNorm_Body.data(),
        aXYZ_Body.data(), aXYZ_Body.size()/3,
        aTri_Body.data(), aTri_Body.size()/3);
    bvh_Body.Init(
        aXYZ_Body.data(), aXYZ_Body.size()/3,
        aTri_Body.data(), aTri_Body.size()/3,
        0.01);
  }
  void Project(double *aXYZt, unsigned int nXYZt) {
    Project_PointsIncludedInBVH_Outside_Cache(
        aXYZt, aInfoNearest_Cloth,
        nXYZt,
        contact_clearance, bvh_Body,
        aXYZ_Body.data(), aXYZ_Body.size() / 3,
        aTri_Body.data(), aTri_Body.size() / 3,
        aNorm_Body.data(), rad_explore);
  }
public:
  const double contact_clearance = 0.001;
  const double rad_explore = 0.1;
  std::vector<CInfoNearest<double>> aInfoNearest_Cloth;
  std::vector<double> aXYZ_Body, aNorm_Body;
  std::vector<unsigned int> aTri_Body;
  CBVH_MeshTri3D<CBV3d_Sphere, double> bvh_Body;
};

/**
 *
 */
class CProjector_RigMesh{
public:
  void UpdatePose(bool isUpdateTopo){
    UpdateBoneRotTrans(aBone);
    SkinningSparse_LBS(aXYZ1_Body,
        aXYZ0_Body, aBone, aSkinningSparseWeight, aSkinningSparseIdBone);
    if( isUpdateTopo ){
      bvh_Body.Init(
          aXYZ1_Body.data(), aXYZ1_Body.size()/3,
          aTri_Body.data(), aTri_Body.size()/3,
          contact_clearance);
    }
    bvh_Body.UpdateGeometry(
        aXYZ1_Body.data(), aXYZ1_Body.size()/3,
        aTri_Body.data(), aTri_Body.size()/3,
        contact_clearance);
    aNorm_Body.resize(aXYZ1_Body.size());
    Normal_MeshTri3D(
        aNorm_Body.data(),
        aXYZ1_Body.data(), aXYZ1_Body.size()/3,
        aTri_Body.data(), aTri_Body.size()/3);
  }
  void Project(double* aXYZt, unsigned int nXYZt){
    Project_PointsIncludedInBVH_Outside_Cache(
        aXYZt, aInfoNearest_Cloth,
        nXYZt,
        contact_clearance, bvh_Body,
        aXYZ1_Body.data(), aXYZ1_Body.size() / 3,
        aTri_Body.data(), aTri_Body.size() / 3,
        aNorm_Body.data(), rad_explore);
  }
public:
  const double contact_clearance = 0.001;
  const double rad_explore = 0.1;
  //
  std::vector<CInfoNearest<double>> aInfoNearest_Cloth;
  //
  std::vector<double> aXYZ0_Body, aXYZ1_Body, aNorm_Body;
  std::vector<unsigned int> aTri_Body;
  std::vector<CRigBone> aBone;
  std::vector<double> aSkinningSparseWeight;
  std::vector<unsigned int> aSkinningSparseIdBone;
  CBVH_MeshTri3D<CBV3d_Sphere,double> bvh_Body;
};

/**
 * @param[in] aXYZ_Contact the array of 3D coordinate of the contact target
 * @param[in] dt time step
 * @param[in] bend_stiff_ratio bending stiffness ratio of the clothing minimium:0 maximum:1
 */
template <typename PROJECTOR>
void StepTime_PbdClothSim(
    std::vector<double>& aXYZ, // deformed vertex positions
    std::vector<double>& aXYZt,
    std::vector<double>& aUVW, // deformed vertex velocity
    const std::vector<int>& aBCFlag,  // boundary condition flag (0:free 1:fixed)
    const std::vector<CDynTri>& aETri,
    const std::vector<CVec2d>& aVec2,
    const std::vector<unsigned int>& aLine,
    PROJECTOR& projector,
    const double dt,
    const double gravity[3],
    const double bend_stiff_ratio)
{
  PBD_Pre3D(
      aXYZt,
      dt, gravity, aXYZ, aUVW, aBCFlag);
  PBD_TriStrain(
      aXYZt.data(),
      aXYZt.size()/3, aETri, aVec2);
  PBD_Bend(
      aXYZt.data(),
      aXYZt.size()/3, aETri, aVec2, bend_stiff_ratio);
  PBD_Seam(
      aXYZt.data(),
      aXYZt.size()/3, aLine.data(), aLine.size()/2);
  projector.Project(aXYZt.data(), aXYZt.size()/3);
  PBD_Post(
      aXYZ, aUVW,
      dt, aXYZt, aBCFlag);
}

}



#endif /* garment_h */
