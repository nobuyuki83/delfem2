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

#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/vec3.h"

#include "delfem2/srchuni_v3.h"

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
    assert( aBB[ibvh].is_active );
    const unsigned int itri0 = ichild0;
    const unsigned int i0 = aTri[itri0*3+0];
    const unsigned int i1 = aTri[itri0*3+1];
    const unsigned int i2 = aTri[itri0*3+2];
    const CVector3 p0(aXYZ[i0*3+0]-px, aXYZ[i0*3+1]-py, aXYZ[i0*3+2]-pz);
    const CVector3 p1(aXYZ[i1*3+0]-px, aXYZ[i1*3+1]-py, aXYZ[i1*3+2]-pz);
    const CVector3 p2(aXYZ[i2*3+0]-px, aXYZ[i2*3+1]-py, aXYZ[i2*3+2]-pz);
    double r0,r1;
    CVector3 p_min = Nearest_Origin_Tri(r0,r1, p0,p1,p2);
    assert( r0 > -1.0e-10 && r1 > -1.0e-10 && (1-r0-r1) > -1.0e-10 );
    double dist = p_min.Length();
    if( dist_min<0 || dist < dist_min ){
      dist_min = dist;
      pes = CPointElemSurf(itri0,r0,r1);
    }
    return;
  }
  /////
  BVH_NearestPoint_MeshTri3D(dist_min,pes, px,py,pz,aXYZ,aTri, ichild0,aBVH,aBB);
  BVH_NearestPoint_MeshTri3D(dist_min,pes, px,py,pz,aXYZ,aTri, ichild1,aBVH,aBB);
}

template <typename T>
class CBVH_MeshTri3D
{
public:
  CBVH_MeshTri3D(){}
  void Init(const std::vector<double>& aXYZ,
            const std::vector<unsigned int>& aTri,
            double margin)
  {
    assert( margin >= 0 );
    { // make BVH topology
      const int ntri = aTri.size()/3;
      std::vector<double> aElemCenter(ntri*3);
      for(int itri=0;itri<ntri;++itri){
        int i0 = aTri[itri*3+0];
        int i1 = aTri[itri*3+1];
        int i2 = aTri[itri*3+2];
        double x0 = (aXYZ[i0*3+0]+aXYZ[i1*3+0]+aXYZ[i2*3+0])/3.0;
        double y0 = (aXYZ[i0*3+1]+aXYZ[i1*3+1]+aXYZ[i2*3+1])/3.0;
        double z0 = (aXYZ[i0*3+2]+aXYZ[i1*3+2]+aXYZ[i2*3+2])/3.0;
        aElemCenter[itri*3+0] = x0;
        aElemCenter[itri*3+1] = y0;
        aElemCenter[itri*3+2] = z0;
      }
      std::vector<int> aTriSurRel;
      makeSurroundingRelationship(aTriSurRel,
                                  aTri.data(), aTri.size()/3,
                                  MESHELEM_TRI, aXYZ.size()/3);
      iroot_bvh = BVH_MakeTreeTopology(aNodeBVH,
                                       3,aTriSurRel,
                                       aElemCenter);
    }
    this->UpdateGeometry(aXYZ,aTri,margin);
    assert( aBB_BVH.size() == aNodeBVH.size() );
  }
  void UpdateGeometry(const std::vector<double>& aXYZ,
                      const std::vector<unsigned int>& aTri,
                      double margin)
  {
    assert( margin >= 0 );
    BVH_BuildBVHGeometry(iroot_bvh,
                         margin,
                         aXYZ,aTri,3,aNodeBVH,aBB_BVH);
    assert( aBB_BVH.size() == aNodeBVH.size() );
  }
  CPointElemSurf NearestPoint_IncludedInBVH(const CVector3& p0,
                                            const std::vector<double>& aXYZ,
                                            const std::vector<unsigned int>& aTri){
    assert( aBB_BVH.size() == aNodeBVH.size() );
    std::vector<int> aIndTri_Cand;
    BVH_GetIndElem_IncludePoint(aIndTri_Cand,
                                p0.x, p0.y, p0.z,
                                iroot_bvh,
                                aNodeBVH,aBB_BVH);
    return Nearest_Point_MeshTri3DPart(p0,aXYZ,aTri,
                                       aIndTri_Cand);
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


#endif
