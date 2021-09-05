/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/vec3.h"
#include "delfem2/geoproximity3_v3.h"
#include "delfem2/gridcube.h"
#include <cstdio>
#include <climits>

void delfem2::Pick_CubeGrid(
    unsigned int& icube_pic,
    int& iface_pic,
    const double src_pic_[3],
    const double dir_pic_[3],
    double elen,
    const double org_[3],
    const std::vector<CCubeGrid>& aCube)
{
  const int noelElemFace_Vox[8][4] = { // this numbering is corresponds to VTK_VOX
    { 0, 4, 6, 2 }, // -x
    { 1, 3, 7, 5 }, // +x
    { 0, 1, 5, 4 }, // -y
    { 2, 6, 7, 3 }, // +y
    { 0, 2, 3, 1 }, // -z
    { 4, 5, 7, 6 }  // +z
  };
  const CVec3d normalHexFace[6] = {
    CVec3d(-1, 0, 0),
    CVec3d(+1, 0, 0),
    CVec3d( 0,-1, 0),
    CVec3d( 0,+1, 0),
    CVec3d( 0, 0,-1),
    CVec3d( 0, 0,+1)
  };
  
  CVec3d dir_pic(dir_pic_), org(org_), src_pic(src_pic_);
  //
  icube_pic = UINT_MAX;
  double depth_min = 0;
  for(unsigned int ivox=0;ivox<aCube.size();++ivox){
    if( !aCube[ivox].is_active ) continue;
    int ih = aCube[ivox].ivx;
    int jh = aCube[ivox].ivy;
    int kh = aCube[ivox].ivz;
    CVec3d cnt =  org + elen*CVec3d(ih+0.5,jh+0.5,kh+0.5);
    {
      CVec3d q = nearest_Line_Point(cnt, src_pic, dir_pic);
      if( (q-cnt).norm() > elen  ) continue;
    }
    CVec3d aP[8] = {
      org + elen*CVec3d(ih+0,jh+0,kh+0),
      org + elen*CVec3d(ih+1,jh+0,kh+0),
      org + elen*CVec3d(ih+0,jh+1,kh+0),
      org + elen*CVec3d(ih+1,jh+1,kh+0),
      org + elen*CVec3d(ih+0,jh+0,kh+1),
      org + elen*CVec3d(ih+1,jh+0,kh+1),
      org + elen*CVec3d(ih+0,jh+1,kh+1),
      org + elen*CVec3d(ih+1,jh+1,kh+1) };
    for(int iface=0;iface<6;++iface){
      const CVec3d& n = normalHexFace[iface];
      const CVec3d& p0 = aP[noelElemFace_Vox[iface][0]];
      const CVec3d& p1 = aP[noelElemFace_Vox[iface][1]];
      //      const CVector3& p2 = aP[noelHexFace[iface][2]];
      const CVec3d& p3 = aP[noelElemFace_Vox[iface][3]];
      const CVec3d pi = intersection_Plane_Line(p0,n, src_pic,dir_pic);
      const double r0 = (pi-p0).dot(p1-p0)/(elen*elen);
      const double r1 = (pi-p0).dot(p3-p0)/(elen*elen);
      if( r0>0 && r0<1 && r1>0 && r1< 1 ){
        double depth = (pi-src_pic).dot(dir_pic)/dir_pic.squaredNorm();
        if( icube_pic == UINT_MAX || depth < depth_min ){
          depth_min = depth;
          icube_pic = ivox;
          iface_pic = iface;
        }
      }
    }
  }
}

void delfem2::Adj_CubeGrid(
    int& ivx, int& ivy, int& ivz,
    unsigned int icube,
    int iface,
    std::vector<CCubeGrid>& aCube)
{
  ivx = aCube[icube].ivx;
  ivy = aCube[icube].ivy;
  ivz = aCube[icube].ivz;
  if( iface == 0 ){ ivx -= 1; }
  if( iface == 1 ){ ivx += 1; }
  if( iface == 2 ){ ivy -= 1; }
  if( iface == 3 ){ ivy += 1; }
  if( iface == 4 ){ ivz -= 1; }
  if( iface == 5 ){ ivz += 1; }
}

void delfem2::Add_CubeGrid
 (std::vector<CCubeGrid>& aCube,
  int ivx1, int ivy1, int ivz1)
{
  for(auto & ic : aCube){
    if( ic.ivx == ivx1 && ic.ivy == ivy1 && ic.ivz == ivz1 ){
      if( ic.is_active ){
        return;
      }
      else{
        ic.is_active = true;
        return;
      }
    }
  }
  {
    CCubeGrid v(ivx1,ivy1,ivz1);
    aCube.push_back(v);
  }
}

void delfem2::Del_CubeGrid
 (std::vector<CCubeGrid>& aCube,
  int i1, int j1, int k1)
{
  for(auto & ic : aCube){
    if( ic.ivx == i1 && ic.ivy == j1 && ic.ivz == k1 ){
      if( ic.is_active ){
        ic.is_active = false;
        return;
      }
      else{
        return;
      }
    }
  }
}


void delfem2::AABB_CubeGrid
 (int aabb[6],
  const std::vector<CCubeGrid>& aCube)
{
  if( aCube.empty() ){
    aabb[0] = +1; aabb[1] = -1;
    aabb[2] = +1; aabb[3] = -1;
    aabb[4] = +1; aabb[5] = -1;
    return;
  }
  aabb[0] = aCube[0].ivx;  aabb[1] = aabb[0]+1;
  aabb[2] = aCube[0].ivy;  aabb[3] = aabb[0]+1;
  aabb[4] = aCube[0].ivz;  aabb[5] = aabb[0]+1;
  
  for(std::size_t ic=1;ic<aCube.size();++ic){
    if( aCube[ic].ivx+0 < aabb[0] ){ aabb[0] = aCube[ic].ivx+0; }
    if( aCube[ic].ivx+1 > aabb[1] ){ aabb[1] = aCube[ic].ivx+1; }
    if( aCube[ic].ivy+0 < aabb[2] ){ aabb[2] = aCube[ic].ivy+0; }
    if( aCube[ic].ivy+1 > aabb[3] ){ aabb[3] = aCube[ic].ivy+1; }
    if( aCube[ic].ivz+0 < aabb[4] ){ aabb[4] = aCube[ic].ivz+0; }
    if( aCube[ic].ivz+1 > aabb[5] ){ aabb[5] = aCube[ic].ivz+1; }
  }
}

