/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef SRCH_V3BVHMSHTOPO_H
#define SRCH_V3BVHMSHTOPO_H

#include <math.h>
#include <vector>

#include "delfem2/bvh.h"
#include "delfem2/mshtopo.h" // sourrounding relationship
#include "delfem2/vec3.h"

#include "delfem2/srchuni_v3.h" // CPointElemSurf

// potential maximum distance of the nearest point
template <typename T>
void BVH_NearestPoint_MeshTri3D
(double& dist_min,
 CPointElemSurf& pes,
 /////
 double px, double py, double pz,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<T>& aBB)
{
  double min0, max0;
  aBB[ibvh].Range_DistToPoint(min0,max0, px,py,pz);
  assert( min0 >= 0 && max0 >= min0 );
  ////
  if( dist_min>=0 && min0>dist_min ){ return; }
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf
    const unsigned int itri0 = ichild0;
    CPointElemSurf pes_tmp;
    double dist = DistanceToTri(pes_tmp,
                                CVector3(px,py,pz),
                                itri0, aXYZ,aTri);
    if( dist_min<0 || dist < dist_min ){
      dist_min = dist;
      pes = pes_tmp;
    }
    return;
  }
  /////
  BVH_NearestPoint_MeshTri3D(dist_min,pes, px,py,pz,aXYZ,aTri, ichild0,aBVH,aBB);
  BVH_NearestPoint_MeshTri3D(dist_min,pes, px,py,pz,aXYZ,aTri, ichild1,aBVH,aBB);
}

// potential maximum distance of the nearest point
template <typename T>
void BVH_NearestPoint_IncludedInBVH_MeshTri3D
(double& dist_tri, // minimum distance to triangle
 double& dist_bv, // minimum distance to leaf bounding volume
 CPointElemSurf& pes,
 /////
 double px, double py, double pz,
 double rad_exp, // exploring distance
 const double* aXYZ, unsigned int nXYZ,
 const unsigned int* aTri, unsigned int nTri,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<T>& aBB)
{
  double min0,max0;
  aBB[ibvh].Range_DistToPoint(min0,max0,
                              px,py,pz);
  if( min0 > rad_exp ){ return; }
  ////
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf
    if( min0 < dist_bv ){ dist_bv = min0; }
    if( min0 == 0.0 ){
      dist_bv = 0.0;
      const unsigned int itri0 = ichild0;
      CPointElemSurf pes_tmp;
      const double dist0 = DistanceToTri(pes_tmp,
                                         CVector3(px,py,pz),
                                         itri0, aXYZ,nXYZ, aTri,nTri);
      if( dist_tri<0 || dist0 < dist_tri ){
        dist_tri = dist0;
        pes = pes_tmp;
      }
    }
    return;
  }
  /////
  BVH_NearestPoint_IncludedInBVH_MeshTri3D(dist_tri,dist_bv, pes,
                                           px,py,pz,rad_exp,
                                           aXYZ,nXYZ,aTri,nTri,
                                           ichild0,aBVH,aBB);
  BVH_NearestPoint_IncludedInBVH_MeshTri3D(dist_tri,dist_bv, pes,
                                           px,py,pz,rad_exp,
                                           aXYZ,nXYZ,aTri,nTri,
                                           ichild1,aBVH,aBB);
}

// ----------------------

template <typename T>
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
        std::vector<int> aTriSurRel;
        makeSurroundingRelationship(aTriSurRel,
                                    pTri, nTri,
                                    delfem2::MESHELEM_TRI, nXYZ);
        iroot_bvh = BVH_MakeTreeTopology(aNodeBVH,
                                        3,aTriSurRel,
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
    BVH_BuildBVHGeometry(iroot_bvh,
                         margin,
                         pXYZ,nXYZ,
                         pTri,3,nTri,
                         aNodeBVH,aBB_BVH);
    assert( aBB_BVH.size() == aNodeBVH.size() );
  }
  double Nearest_Point_IncludedInBVH(CPointElemSurf& pes,
                                     const CVector3& p0,
                                     double rad_exp, // look leaf inside this radius
                                     const double* aXYZ, unsigned int nXYZ,
                                     const unsigned int* aTri, unsigned int nTri) const{
    assert( aBB_BVH.size() == aNodeBVH.size() );
    double dist = -1, dist_min = rad_exp;
    pes.itri = -1;
    BVH_NearestPoint_IncludedInBVH_MeshTri3D(dist,dist_min, pes,
                                             p0.x, p0.y, p0.z, rad_exp,
                                             aXYZ, nXYZ, aTri, nTri,
                                             iroot_bvh, aNodeBVH, aBB_BVH);
    if( pes.itri == -1 ){ return dist_min; }
    return dist;
    /*
    std::vector<int> aIndTri_Cand;
    BVH_GetIndElem_IncludePoint(aIndTri_Cand,
                                p0.x, p0.y, p0.z,
                                iroot_bvh,
                                aNodeBVH,aBB_BVH);
    return Nearest_Point_MeshTri3DPart(p0,aXYZ,aTri,
                                       aIndTri_Cand);
     */
  }
  CPointElemSurf NearestPoint_Global(const CVector3& p0,
                                     const std::vector<double>& aXYZ,
                                     const std::vector<unsigned int>& aTri) const {
    assert( aBB_BVH.size() == aNodeBVH.size() );
    CPointElemSurf pes;
    double dist_min = -1;
    BVH_NearestPoint_MeshTri3D(dist_min, pes,
                               p0.x, p0.y, p0.z,
                               aXYZ, aTri, iroot_bvh, aNodeBVH, aBB_BVH);
    return pes;
  }
  // inside positive
  double SignedDistanceFunction(CVector3& n0,
                                const CVector3& p0,
                                const std::vector<double>& aXYZ,
                                const std::vector<unsigned int>& aTri,
                                const std::vector<double>& aNorm) const
  {
    assert( aBB_BVH.size() == aNodeBVH.size() );
    CPointElemSurf pes;
    {
      double dist_min = -1;
      BVH_NearestPoint_MeshTri3D(dist_min, pes,
                                 p0.x, p0.y, p0.z,
                                 aXYZ, aTri,
                                 iroot_bvh, aNodeBVH, aBB_BVH);
    }
    const CVector3 q0 = pes.Pos_Tri(aXYZ, aTri);
    double dist = (q0-p0).Length();
    if( !aBB_BVH[iroot_bvh].isInclude_Point(p0.x,p0.y,p0.z) ){ // outside
      n0 = (p0-q0).Normalize();
      return -dist;
    }
    const CVector3 n1 = pes.UNorm_Tri(aXYZ, aTri, aNorm);
    if( dist < 1.0e-6 ){
      n0 = n1;
      if( (q0-p0)*n1 > 0 ){ return dist; } //inside
      return -dist; // outside
    }
    CVector3 dir = (cg_Tri(pes.itri, aTri, aXYZ)-p0).Normalize();
    if( (q0-p0)*n1 < 0 ){ dir = -dir; } // probaby outside so shoot ray outside
    std::vector<int> aIndElem;
    double ps0[3]; p0.CopyValueTo(ps0);
    double pd0[3]; dir.CopyValueTo(pd0);
    BVH_GetIndElem_IntersectRay(aIndElem, ps0, pd0,
                                iroot_bvh, aNodeBVH, aBB_BVH);
    std::map<double,CPointElemSurf> mapDepthPES1;
    IntersectionRay_MeshTri3DPart(mapDepthPES1,
                                  p0, dir,
                                  aTri, aXYZ, aIndElem);
    if( mapDepthPES1.size() %2 == 0 ){ // outside
      n0 = (p0-q0).Normalize();
      return -dist;
    }
    n0 = (q0-p0).Normalize();
    return +dist;
  }
public:
  int iroot_bvh;
  std::vector<CNodeBVH> aNodeBVH; // array of BVH node
  std::vector<T> aBB_BVH;
};


template <typename T>
void Project_PointsIncludedInBVH_Outside
(std::vector<double>& aXYZt,
 double cc,
 const CBVH_MeshTri3D<T>& bvh,
 const std::vector<double>& aXYZ0,
 const std::vector<unsigned int>& aTri0,
 const std::vector<double>& aNorm0)
{
  for(unsigned int ip=0;ip<aXYZt.size()/3;++ip){
    CVector3 p0(aXYZt[ip*3+0], aXYZt[ip*3+1], aXYZt[ip*3+2] );
    CPointElemSurf pes;
    bvh.Nearest_Point_IncludedInBVH(pes,
                                    p0, 0.0,
                                    aXYZ0.data(),aXYZ0.size()/3,
                                    aTri0.data(),aTri0.size()/3);
    /////
    if( pes.itri == -1 ){ continue; }
    CVector3 n0;
    double sdf = SDFNormal_NearestPoint(n0,
                                        p0,pes,aXYZ0,aTri0,aNorm0);
//    std::cout << sdf+cc << std::endl;
    if( sdf+cc < 0 ) continue;
    aXYZt[ip*3+0] += (sdf+cc)*n0.x;
    aXYZt[ip*3+1] += (sdf+cc)*n0.y;
    aXYZt[ip*3+2] += (sdf+cc)*n0.z;
  }
}

class CInfoNearest
{
public:
  CInfoNearest(){
    is_active = false;
  }
public:
  CPointElemSurf pes;
  CVector3 pos;
  double sdf;
  bool is_active;
};

template <typename T>
void Project_PointsIncludedInBVH_Outside_Cache
(double* aXYZt,
 std::vector<CInfoNearest>& aInfoNearest,
 unsigned int nXYZt,
 double cc,
 const CBVH_MeshTri3D<T>& bvh,
 const double* pXYZ0, unsigned int nXYZ0,
 const unsigned int* pTri0, unsigned int nTri0,
 const double* pNorm0,
 double rad_explore)
{
  const unsigned int np = nXYZt;
  aInfoNearest.resize(np);
  for(unsigned int ip=0;ip<np;++ip){
    CVector3 p0(aXYZt[ip*3+0], aXYZt[ip*3+1], aXYZt[ip*3+2] );
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
    if( aInfoNearest[ip].pes.itri == -1 ){
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
    CVector3 n0;
    double sdf = SDFNormal_NearestPoint(n0,
                                        aInfoNearest[ip].pos, aInfoNearest[ip].pes,
                                        pXYZ0, nXYZ0,
                                        pTri0, nTri0,
                                        pNorm0);
    aInfoNearest[ip].sdf = sdf;
    aInfoNearest[ip].is_active = true;
    if( sdf+cc < 0 ) continue;
    aXYZt[ip*3+0] += (sdf+cc)*n0.x;
    aXYZt[ip*3+1] += (sdf+cc)*n0.y;
    aXYZt[ip*3+2] += (sdf+cc)*n0.z;
  }
}


#endif
