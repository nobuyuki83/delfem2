/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_SRCH_V3BVHMSHTOPO_H
#define DFM2_SRCH_V3BVHMSHTOPO_H

#include "delfem2/srchuni_v3.h" // CPointElemSurf
#include "delfem2/bvh.h"
#include "delfem2/mshtopo.h" // sourrounding relationship
#include "delfem2/vec3.h"
#include <math.h>
#include <vector>

namespace delfem2 {
  
// potential maximum distance of the nearest point
template <typename T, typename REAL>
void BVH_NearestPoint_MeshTri3D (
    double& dist_min,
    CPointElemSurf<REAL>& pes,
    //
    double px, double py, double pz,
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<T>& aBB);
  
// potential maximum distance of the nearest point
template <typename T, typename REAL>
void BVH_NearestPoint_IncludedInBVH_MeshTri3D(
    double& dist_tri, // minimum distance to triangle
    double& dist_bv, // minimum distance to leaf bounding volume
    CPointElemSurf<REAL>& pes,
    //
    double px, double py, double pz,
    double rad_exp, // exploring distance
    const double* aXYZ, unsigned int nXYZ,
    const unsigned int* aTri, unsigned int nTri,
    int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<T>& aBB);

template <typename BV, typename REAL>
class CBVH_MeshTri3D
{
public:
  CBVH_MeshTri3D(){}
  void Init(const double* pXYZ, unsigned int nXYZ,
            const unsigned int* pTri, unsigned int nTri,
            double margin)
  {
    assert( margin >= 0 );
    { // make BVH topology
      std::vector<double> aElemCenter(nTri*3);
      for(unsigned int itri=0;itri<nTri;++itri){
        int i0 = pTri[itri*3+0];
        int i1 = pTri[itri*3+1];
        int i2 = pTri[itri*3+2];
        double x0 = (pXYZ[i0*3+0]+pXYZ[i1*3+0]+pXYZ[i2*3+0])/3.0;
        double y0 = (pXYZ[i0*3+1]+pXYZ[i1*3+1]+pXYZ[i2*3+1])/3.0;
        double z0 = (pXYZ[i0*3+2]+pXYZ[i1*3+2]+pXYZ[i2*3+2])/3.0;
        aElemCenter[itri*3+0] = x0;
        aElemCenter[itri*3+1] = y0;
        aElemCenter[itri*3+2] = z0;
      }
      {
        std::vector<unsigned int> aTriSuTri;
        ElSuEl_MeshElem(aTriSuTri,
                        pTri, nTri,
                        delfem2::MESHELEM_TRI, nXYZ);
        iroot_bvh = BVHTopology_TopDown_MeshElem(aNodeBVH,
                                                 3,aTriSuTri,
                                                 aElemCenter);
      }
    }
    this->UpdateGeometry(pXYZ, nXYZ,
                         pTri, nTri,
                         margin);
    assert( aBB_BVH.size() == aNodeBVH.size() );
  }
  void UpdateGeometry(const double* pXYZ, unsigned int nXYZ,
                      const unsigned int* pTri, unsigned int nTri,
                      double margin)
  {
    assert( margin >= 0 );
    CLeafVolumeMaker_Mesh<BV,REAL> lvm(
        margin,
        pXYZ,nXYZ,
        pTri,nTri,3);
    BVH_BuildBVHGeometry(
        aBB_BVH,
        //
        iroot_bvh, aNodeBVH,
        lvm);
    assert( aBB_BVH.size() == aNodeBVH.size() );
  }
  double Nearest_Point_IncludedInBVH(
      CPointElemSurf<REAL>& pes,
      const CVec3<REAL>& p0,
      double rad_exp, // look leaf inside this radius
      const double* aXYZ, unsigned int nXYZ,
      const unsigned int* aTri, unsigned int nTri) const{
    assert( aBB_BVH.size() == aNodeBVH.size() );
    double dist = -1, dist_min = rad_exp;
    pes.itri = UINT_MAX;
    delfem2::BVH_NearestPoint_IncludedInBVH_MeshTri3D(
        dist,dist_min, pes,
        p0.x(), p0.y(), p0.z(), rad_exp,
        aXYZ, nXYZ, aTri, nTri,
        iroot_bvh, aNodeBVH, aBB_BVH);
    if( pes.itri == UINT_MAX ){ return dist_min; }
    return dist;
  }
  CPointElemSurf<REAL> NearestPoint_Global(
      const CVec3<REAL>& p0,
      const std::vector<double>& aXYZ,
      const std::vector<unsigned int>& aTri) const {
    assert( aBB_BVH.size() == aNodeBVH.size() );
    CPointElemSurf<REAL> pes;
    double dist_min = -1;
    BVH_NearestPoint_MeshTri3D(dist_min, pes,
                               p0.x(), p0.y(), p0.z(),
                               aXYZ, aTri, iroot_bvh, aNodeBVH, aBB_BVH);
    return pes;
  }
    // inside positive
  double SignedDistanceFunction(
      CVec3<REAL>& n0,
      //
      const CVec3<REAL>& p0,
      const std::vector<double>& aXYZ,
      const std::vector<unsigned int>& aTri,
      const std::vector<double>& aNorm) const
  {
    assert( aBB_BVH.size() == aNodeBVH.size() );
    CPointElemSurf<REAL> pes;
    {
      double dist_min = -1;
      delfem2::BVH_NearestPoint_MeshTri3D(dist_min, pes,
                                          p0.x(), p0.y(), p0.z(),
                                          aXYZ, aTri,
                                          iroot_bvh, aNodeBVH, aBB_BVH);
    }
    const CVec3<REAL> q0 = pes.Pos_Tri(aXYZ, aTri);
    double dist = (q0-p0).Length();
    if( !aBB_BVH[iroot_bvh].isInclude_Point(p0.x(),p0.y(),p0.z()) ){ // outside
      n0 = (p0-q0).Normalize();
      return -dist;
    }
    const CVec3<REAL> n1 = pes.UNorm_Tri(aXYZ, aTri, aNorm);
    if( dist < 1.0e-6 ){
      n0 = n1;
      if( (q0-p0)*n1 > 0 ){ return dist; } //inside
      return -dist; // outside
    }
    CVec3<REAL> dir = (CG_Tri3(pes.itri, aTri, aXYZ)-p0).Normalize();
    if( (q0-p0)*n1 < 0 ){ dir = -dir; } // probaby outside so shoot ray outside
    std::vector<unsigned int> aIndElem;
    BVH_GetIndElem_Predicate(aIndElem,
        CIsBV_IntersectRay<BV>(p0.p, dir.p),
        iroot_bvh, aNodeBVH, aBB_BVH);
    std::map<double,CPointElemSurf<REAL>> mapDepthPES1;
    IntersectionRay_MeshTri3DPart(mapDepthPES1,
                                  p0, dir,
                                  aTri, aXYZ, aIndElem,
                                  0.0);
    if( mapDepthPES1.size() %2 == 0 ){ // outside
      n0 = (p0-q0).Normalize();
      return -dist;
    }
    n0 = (q0-p0).Normalize();
    return +dist;
  }
public:
  int iroot_bvh;
  std::vector<delfem2::CNodeBVH2> aNodeBVH; // array of BVH node
  std::vector<BV> aBB_BVH;
};
  
template <typename T, typename REAL>
void Project_PointsIncludedInBVH_Outside(
    std::vector<double>& aXYZt,
    double cc,
    const CBVH_MeshTri3D<T,REAL>& bvh,
    const std::vector<double>& aXYZ0,
    const std::vector<unsigned int>& aTri0,
    const std::vector<double>& aNorm0);
  
template <typename REAL>
class CInfoNearest
{
  public:
    CInfoNearest(){
      is_active = false;
    }
  public:
    CPointElemSurf<REAL> pes;
    CVec3<REAL> pos;
    double sdf;
    bool is_active;
};

template <typename T, typename REAL>
void Project_PointsIncludedInBVH_Outside_Cache(double* aXYZt,
                                               std::vector<CInfoNearest<REAL>>& aInfoNearest,
                                               unsigned int nXYZt,
                                               double cc,
                                               const CBVH_MeshTri3D<T,REAL>& bvh,
                                               const double* pXYZ0, unsigned int nXYZ0,
                                               const unsigned int* pTri0, unsigned int nTri0,
                                               const double* pNorm0,
                                               double rad_explore);
  
} // namespace delfem2

// ---------------------------------------------------
// implemnatation for the template functions

// potential maximum distance of the nearest point
template <typename T, typename REAL>
void delfem2::BVH_NearestPoint_MeshTri3D(
    double& dist_min,
    CPointElemSurf<REAL>& pes,
    //
    double px, double py, double pz,
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<T>& aBB)
{
  double min0, max0;
  aBB[ibvh].Range_DistToPoint(min0,max0, px,py,pz);
  assert( min0 >= 0 && max0 >= min0 );
  //
  if( dist_min>=0 && min0>dist_min ){ return; }
  const unsigned int ichild0 = aBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == UINT_MAX ){ // leaf
    const unsigned int itri0 = ichild0;
    CPointElemSurf<REAL> pes_tmp;
    double dist = DistanceToTri(
        pes_tmp,
        CVec3<REAL>(px,py,pz),
        itri0, aXYZ,aTri);
    if( dist_min<0 || dist < dist_min ){
      dist_min = dist;
      pes = pes_tmp;
    }
    return;
  }
  //
  BVH_NearestPoint_MeshTri3D(dist_min,pes, px,py,pz,aXYZ,aTri, ichild0,aBVH,aBB);
  BVH_NearestPoint_MeshTri3D(dist_min,pes, px,py,pz,aXYZ,aTri, ichild1,aBVH,aBB);
}

// potential maximum distance of the nearest point
template <typename BV, typename REAL>
void delfem2::BVH_NearestPoint_IncludedInBVH_MeshTri3D(
    double& dist_tri, // minimum distance to triangle
    double& dist_bv, // minimum distance to leaf bounding volume
    CPointElemSurf<REAL>& pes,
    //
    double px, double py, double pz,
    double rad_exp, // exploring distance
    const double* aXYZ, unsigned int nXYZ,
    const unsigned int* aTri, unsigned int nTri,
    int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<BV>& aBB)
{
  assert(ibvh>=0&&ibvh<aBB.size());
  double min0,max0;
  aBB[ibvh].Range_DistToPoint(min0,max0,
                              px,py,pz);
  if( min0 > rad_exp ){ return; }
  //
  const unsigned int ichild0 = aBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == UINT_MAX ){ // leaf
    if( min0 < dist_bv ){ dist_bv = min0; }
    if( min0 == 0.0 ){
      dist_bv = 0.0;
      const unsigned int itri0 = ichild0;
      CPointElemSurf<REAL> pes_tmp;
      const double dist0 = DistanceToTri(
          pes_tmp,
          //
          CVec3<REAL>(px,py,pz),
          itri0, aXYZ,nXYZ, aTri,nTri);
      if( dist_tri<0 || dist0 < dist_tri ){
        dist_tri = dist0;
        pes = pes_tmp;
      }
    }
    return;
  }
  BVH_NearestPoint_IncludedInBVH_MeshTri3D(
      dist_tri,dist_bv, pes,
      //
      px,py,pz,rad_exp,
      aXYZ,nXYZ,aTri,nTri,
      ichild0,aBVH,aBB);
  BVH_NearestPoint_IncludedInBVH_MeshTri3D(
      dist_tri,dist_bv, pes,
      //
      px,py,pz,rad_exp,
      aXYZ,nXYZ,aTri,nTri,
      ichild1,aBVH,aBB);
}

// ----------------------

template <typename BV, typename REAL>
void delfem2::Project_PointsIncludedInBVH_Outside(
    std::vector<double>& aXYZt,
    double cc,
    const delfem2::CBVH_MeshTri3D<BV, REAL>& bvh,
    const std::vector<double>& aXYZ0,
    const std::vector<unsigned int>& aTri0,
    const std::vector<double>& aNorm0)
{
  for(unsigned int ip=0;ip<aXYZt.size()/3;++ip){
    CVec3<REAL> p0(aXYZt[ip*3+0], aXYZt[ip*3+1], aXYZt[ip*3+2] );
    CPointElemSurf<REAL> pes;
    bvh.Nearest_Point_IncludedInBVH(pes,
                                    p0, 0.0,
                                    aXYZ0.data(),aXYZ0.size()/3,
                                    aTri0.data(),aTri0.size()/3);
    if( pes.itri == UINT_MAX ){ continue; }
    CVec3<REAL> n0;
    double sdf = SDFNormal_NearestPoint(n0,
                                        p0,pes,aXYZ0,aTri0,aNorm0);
    if( sdf+cc < 0 ) continue;
    aXYZt[ip*3+0] += (sdf+cc)*n0.x();
    aXYZt[ip*3+1] += (sdf+cc)*n0.y();
    aXYZt[ip*3+2] += (sdf+cc)*n0.z();
  }
}

template <typename REAL>
class CInfoNearest
{
public:
  CInfoNearest(){
    is_active = false;
  }
public:
  delfem2::CPointElemSurf<REAL> pes;
  delfem2::CVec3<REAL> pos;
  double sdf;
  bool is_active;
};

/**
 *
 * @tparam BV
 * @tparam REAL
 * @param aXYZt
 * @param aInfoNearest
 * @param nXYZt
 * @param cc (in) contact clearance
 * @param bvh
 * @param pXYZ0
 * @param nXYZ0
 * @param pTri0
 * @param nTri0
 * @param pNorm0
 * @param rad_explore
 */
template <typename BV, typename REAL>
void delfem2::Project_PointsIncludedInBVH_Outside_Cache(
    double* aXYZt,
    std::vector<delfem2::CInfoNearest<REAL>>& aInfoNearest,
    unsigned int nXYZt,
    double cc,
    const delfem2::CBVH_MeshTri3D<BV,REAL>& bvh,
    const double* pXYZ0, unsigned int nXYZ0,
    const unsigned int* pTri0, unsigned int nTri0,
    const double* pNorm0,
    double rad_explore)
{
  const unsigned int np = nXYZt;
  aInfoNearest.resize(np);
  for(unsigned int ip=0;ip<np;++ip){
    CVec3<REAL> p0(aXYZt[ip*3+0], aXYZt[ip*3+1], aXYZt[ip*3+2] );
    if( aInfoNearest[ip].is_active ){
      double dp = Distance(p0,aInfoNearest[ip].pos);
      if( aInfoNearest[ip].sdf + dp + cc < 0 ){
        continue;
      }
    }
    aInfoNearest[ip].pos = p0;
    double dist0 = bvh.Nearest_Point_IncludedInBVH(aInfoNearest[ip].pes,
                                                   aInfoNearest[ip].pos, rad_explore,
                                                   pXYZ0, nXYZ0,
                                                   pTri0, nTri0);
    if( aInfoNearest[ip].pes.itri == UINT_MAX ){
      if( aInfoNearest[ip].is_active ){
        if( aInfoNearest[ip].sdf < 0 ){ aInfoNearest[ip].sdf = -dist0; }
        else{                           aInfoNearest[ip].sdf = +dist0; }
      }
      else{
        aInfoNearest[ip].sdf = -dist0;
        aInfoNearest[ip].is_active = true;
      }
      continue;
    }
    CVec3<REAL> n0;
    double sdf = SDFNormal_NearestPoint(
        n0,
        aInfoNearest[ip].pos,
        aInfoNearest[ip].pes,
        pXYZ0, nXYZ0,
        pTri0, nTri0,
        pNorm0);
    aInfoNearest[ip].sdf = sdf;
    aInfoNearest[ip].is_active = true;
    if( sdf+cc < 0 ) continue;
    aXYZt[ip*3+0] += (sdf+cc)*n0.x();
    aXYZt[ip*3+1] += (sdf+cc)*n0.y();
    aXYZt[ip*3+2] += (sdf+cc)*n0.z();
  }
}


#endif
