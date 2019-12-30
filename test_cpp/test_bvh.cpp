/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec3.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/sdf.h"
#include "delfem2/primitive.h"
#include "delfem2/mshmisc.h"

#include "delfem2/srchuni_v3.h"
#include "delfem2/objfunc_v23.h"
#include "delfem2/srch_v3bvhmshtopo.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ---------------------------------------------

TEST(bvh,inclusion_sphere)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    delfem2::Rotate(aXYZ, 0.2, 0.3, 0.4);
  }
  //  std::cout << "ntri: " << aTri.size()/3 << std::endl;
  dfm2::CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.03);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> udist(-1.0, 1.0);
  for(int itr=0;itr<10000;++itr){
    CVector3 p0(udist(rng), udist(rng), udist(rng));
    for(int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3D_Sphere& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.isInclude_Point(p0.x, p0.y, p0.z);
      if( !is_intersect && node.ichild[1] != -1 ){ // branch
        const int ichild0 = node.ichild[0];
        const int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].isInclude_Point(p0.x, p0.y, p0.z) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].isInclude_Point(p0.x, p0.y, p0.z) );
      }
    }
  }
}

TEST(bvh,inclusion_aabb)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    delfem2::Rotate(aXYZ,
                    0.2, 0.3, 0.4);
  }
  //  std::cout << "ntri: " << aTri.size()/3 << std::endl;
  dfm2::CBVH_MeshTri3D<dfm2::CBV3D_AABB> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.03);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> udist(-1.0, 1.0);
  for(int itr=0;itr<10000;++itr){
    CVector3 p0(udist(rng), udist(rng), udist(rng));
    for(int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3D_AABB& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.isInclude_Point(p0.x, p0.y, p0.z);
      if( !is_intersect && node.ichild[1] != -1 ){ // branch
        const int ichild0 = node.ichild[0];
        const int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].isInclude_Point(p0.x, p0.y, p0.z) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].isInclude_Point(p0.x, p0.y, p0.z) );
      }
    }
  }
}

TEST(bvh,nearestinc_sphere)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    delfem2::Rotate(aXYZ, 0.2, 0.3, 0.4);
  }
//  std::cout << "ntri: " << aTri.size()/3 << std::endl;
  std::vector<double> aNorm(aXYZ.size());
  Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  dfm2::CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.03);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> udist(-5.0, 5.0);
  for(int itr=0;itr<1000;++itr){
    CVector3 p0(udist(rng), udist(rng), udist(rng));
    {
      p0.SetNormalizedVector();
      if( itr % 2 == 0 ){ p0 *= 1.02; } // outside included in bvh
      else{               p0 *= 0.98; } // inside in included in bvh
    }
    CPointElemSurf pes1;
    double dist1 = bvh.Nearest_Point_IncludedInBVH(pes1,p0,0.1,
                                                   aXYZ.data(), aXYZ.size()/3,
                                                   aTri.data(), aTri.size()/3);
    EXPECT_LE( dist1, 0.1 );
    EXPECT_GE( dist1, 0.0 );
    EXPECT_TRUE( pes1.Check(aXYZ, aTri,1.0e-10) );
    CVector3 q1 = pes1.Pos_Tri(aXYZ, aTri);
    {
      CPointElemSurf pes0 = Nearest_Point_MeshTri3D(p0, aXYZ, aTri);
      CVector3 q0 = pes0.Pos_Tri(aXYZ, aTri);
      EXPECT_LT(Distance(q0,q1),1.0e-10);
    }
    CVector3 n0 = pes1.UNorm_Tri(aXYZ, aTri, aNorm);
    EXPECT_EQ( n0*(p0-q1)>0, itr%2==0);
    //////
    {
      CPointElemSurf pes2;
      double dist_tri = -1, dist_bv = 0.1;
      dfm2::BVH_NearestPoint_IncludedInBVH_MeshTri3D(dist_tri, dist_bv, pes2,
                                                     p0.x, p0.y, p0.z, 0.1,
                                                     aXYZ.data(), aXYZ.size()/3,
                                                     aTri.data(), aXYZ.size()/3,
                                                     bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
      CVector3 q2 = pes2.Pos_Tri(aXYZ, aTri);
      EXPECT_LT(Distance(q2,q1),1.0e-10);
    }
  }
}


TEST(bvh,nearest_range) // find global nearest from range
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    delfem2::Rotate(aXYZ,
                    0.2, 0.3, 0.4);
  }
  dfm2::CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.0);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> udist(-5.0, 5.0);
  for(int itr=0;itr<1000;++itr){
    CVector3 p0(udist(rng), udist(rng), udist(rng));
    {
      double dist_min=-1, dist_max = -1;
      dfm2::BVH_Range_DistToNearestPoint(dist_min, dist_max,
                                         p0.x, p0.y, p0.z,
                                         bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
      bool is_max = false;
      for(int it=0;it<aTri.size()/3;++it){
        dfm2::CBV3D_Sphere bb;
        for(int inoel=0;inoel<3;++inoel){
          const int ino0 = aTri[it*3+inoel];
          bb.AddPoint(aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2], 0.0);
        }
        double min0, max0;
        bb.Range_DistToPoint(min0, max0, p0.x, p0.y, p0.z);
        EXPECT_GE( max0, dist_max );
        EXPECT_GE( min0, dist_min );
        if( max0 < dist_max+1.0e-10 ){ is_max = true; }
      }
      EXPECT_TRUE( is_max );
      std::vector<int> aIndElem;
      BVH_GetIndElem_InsideRange(aIndElem,
                                 dist_min,dist_max,
                                 p0.x, p0.y, p0.z,
                                 bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
      EXPECT_GT(aIndElem.size(), 0);
      std::vector<int> aFlg(aTri.size()/3,0);
      for(int iit=0;iit<aIndElem.size();++iit){
        int itri0 = aIndElem[iit];
        aFlg[itri0] = 1;
      }
      for(int itri=0;itri<aTri.size()/3;++itri){
        dfm2::CBV3D_Sphere bb;
        for(int inoel=0;inoel<3;++inoel){
          const int ino0 = aTri[itri*3+inoel];
          bb.AddPoint(aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2], 0.0);
        }
        double min0, max0;
        bb.Range_DistToPoint(min0, max0, p0.x, p0.y, p0.z);
        if( aFlg[itri] == 1 ){ // inside range
          EXPECT_LE(min0,dist_max);
          EXPECT_GE(max0,dist_min);
        }
        else{ // outside range
          EXPECT_TRUE((min0>dist_max)||(max0<dist_min));
        }
      }
    }
  }
}

TEST(bvh,nearest_point) // find global nearest directry
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    delfem2::Rotate(aXYZ, 0.2, 0.3, 0.4);
  }
  dfm2::CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.0);
  for(int itr=0;itr<1000;++itr){
    CVector3 p0;
    {
      p0.x = 10.0*(rand()/(RAND_MAX+1.0)-0.5);
      p0.y = 10.0*(rand()/(RAND_MAX+1.0)-0.5);
      p0.z = 10.0*(rand()/(RAND_MAX+1.0)-0.5);
    }
    CPointElemSurf pes1 = bvh.NearestPoint_Global(p0,aXYZ,aTri);
    EXPECT_TRUE( pes1.Check(aXYZ, aTri,1.0e-10) );
    CVector3 q1 = pes1.Pos_Tri(aXYZ, aTri);
    {
      CPointElemSurf pes0 = Nearest_Point_MeshTri3D(p0, aXYZ, aTri);
      CVector3 q0 = pes0.Pos_Tri(aXYZ, aTri);
      EXPECT_LT(Distance(q0,q1),1.0e-10);
    }
  }
}


TEST(bvh,sdf) // find global nearest directry
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    delfem2::Rotate(aXYZ,
                    0.2, 0.3, 0.4);
  }
  std::vector<double> aNorm(aXYZ.size());
  Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  dfm2::CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.0);
  for(int itr=0;itr<1000;++itr){
    CVector3 p0;
    {
      p0.x = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      p0.y = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      p0.z = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
    }
    if( (p0.Length()-1.0)<1.0e-3 ){ continue; }
    CVector3 n0;
    double sdf = bvh.SignedDistanceFunction(n0,
                                            p0, aXYZ, aTri, aNorm);
    EXPECT_NEAR(1-p0.Length(), sdf, 1.0e-2);
    EXPECT_NEAR(n0*p0.Normalize(), 1.0, 1.0e-2 );
  }
}


TEST(bvh,lineintersection)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    delfem2::Rotate(aXYZ,
                    0.2, 0.3, 0.4);
  }
  std::vector<double> aNorm(aXYZ.size());
  Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  delfem2::CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           1.0e-5);
  for(int itr=0;itr<100;++itr){
    CVector3 s0, d0;
    {
      s0.x = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      s0.y = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      s0.z = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.x = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.y = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.z = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.SetNormalizedVector();
    }
    double ps0[3]; s0.CopyValueTo(ps0);
    double pd0[3]; d0.CopyValueTo(pd0);
    for(int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3D_Sphere& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.IsIntersectLine(ps0,pd0);
      if( !is_intersect && node.ichild[1] != -1 ){ // branch
        const int ichild0 = node.ichild[0];
        const int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].IsIntersectLine(ps0,pd0) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].IsIntersectLine(ps0,pd0) );
      }
    }
    std::vector<int> aIndElem;
    BVH_GetIndElem_IntersectLine(aIndElem, ps0, pd0,
                                 bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
    std::vector<int> aFlg(aTri.size()/3,0);
    for(unsigned int iit=0;iit<aIndElem.size();++iit){
      int itri0 = aIndElem[iit];
      aFlg[itri0] = 1;
    }
    for(unsigned int itri=0;itri<aTri.size()/3;++itri){
      dfm2::CBV3D_Sphere bb;
      for(int inoel=0;inoel<3;++inoel){
        const int ino0 = aTri[itri*3+inoel];
        bb.AddPoint(aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2], 1.0e-5);
      }
      bool res = bb.IsIntersectLine(ps0, pd0);
      EXPECT_EQ( res , aFlg[itri] == 1 );
    }
  }
}


TEST(bvh,rayintersection)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    delfem2::MeshTri3D_Sphere(aXYZ, aTri,
                              1.0, 64, 32);
    delfem2::Rotate(aXYZ,
                    0.2, 0.3, 0.4);
  }
  std::vector<double> aNorm(aXYZ.size());
  Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  delfem2::CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           1.0e-5);
  for(int itr=0;itr<100;++itr){
    CVector3 s0, d0;
    {
      s0.x = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      s0.y = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      s0.z = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.x = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.y = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.z = 3.0*(rand()/(RAND_MAX+1.0)-0.5);
      d0.SetNormalizedVector();
    }
    double ps0[3]; s0.CopyValueTo(ps0);
    double pd0[3]; d0.CopyValueTo(pd0);
    for(int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3D_Sphere& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.IsIntersectRay(ps0,pd0);
      if( !is_intersect && node.ichild[1] != -1 ){ // branch
        const int ichild0 = node.ichild[0];
        const int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].IsIntersectRay(ps0,pd0) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].IsIntersectRay(ps0,pd0) );
      }
    }
    std::vector<int> aIndElem;
    BVH_GetIndElem_IntersectRay(aIndElem, ps0, pd0,
                                bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
    std::vector<int> aFlg(aTri.size()/3,0);
    for(int itri0 : aIndElem){
      aFlg[itri0] = 1;
    }
    for(unsigned int itri=0;itri<aTri.size()/3;++itri){
      dfm2::CBV3D_Sphere bb;
      for(int inoel=0;inoel<3;++inoel){
        const int ino0 = aTri[itri*3+inoel];
        bb.AddPoint(aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2], 1.0e-5);
      }
      bool res = bb.IsIntersectRay(ps0, pd0);
      EXPECT_EQ( res , aFlg[itri] == 1 );
    }
    {
      std::map<double,CPointElemSurf> mapDepthPES0;
      IntersectionRay_MeshTri3D(mapDepthPES0,
                                s0, d0, aTri, aXYZ);
      std::map<double,CPointElemSurf> mapDepthPES1;
      IntersectionRay_MeshTri3DPart(mapDepthPES1,
                                    s0, d0,
                                    aTri, aXYZ, aIndElem);
      EXPECT_EQ(mapDepthPES0.size(),mapDepthPES1.size());
      int N = mapDepthPES0.size();
      auto itr0 = mapDepthPES0.begin();
      auto itr1 = mapDepthPES0.begin();
      for(int i=0;i<N;++i){
        EXPECT_FLOAT_EQ(itr0->first,itr1->first);
        CPointElemSurf pes0 = itr0->second;
        CPointElemSurf pes1 = itr1->second;
        CVector3 q0 = pes0.Pos_Tri(aXYZ, aTri);
        CVector3 q1 = pes1.Pos_Tri(aXYZ, aTri);
        EXPECT_NEAR(Distance(q0,q1), 0.0, 1.0e-10);
      }
    }
  }
}

void mark_child(std::vector<int>& aFlg,
                unsigned int inode0,
                const std::vector<dfm2::CNodeBVH2>& aNode)
{
  assert( inode0 < aNode.size() );
  if( aNode[inode0].ichild[1] == -1 ){ // leaf
    const unsigned int in0 = aNode[inode0].ichild[0];
    assert( in0 < aFlg.size() );
    aFlg[in0] += 1;
    return;
  }
  const unsigned int in0 = aNode[inode0].ichild[0];
  const unsigned int in1 = aNode[inode0].ichild[1];
  mark_child(aFlg, in0, aNode);
  mark_child(aFlg, in1, aNode);
}

TEST(bvh,morton_code)
{
  std::vector<double> aXYZ; // 3d points
  const unsigned int N = 10000;
  aXYZ.resize(N*3);
  const double minmax_xyz[6] = {-1,+1, -1,+1, -1,+1};
  dfm2::CBV3D_AABB bb(minmax_xyz);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> udist(0.0, 1.0);
  for(int i=0;i<N;++i){
    aXYZ[i*3+0] = (bb.x_max -  bb.x_min) * udist(rng) + bb.x_min;
    aXYZ[i*3+1] = (bb.y_max -  bb.y_min) * udist(rng) + bb.y_min;
    aXYZ[i*3+2] = (bb.z_max -  bb.z_min) * udist(rng) + bb.z_min;
  }
  std::vector<unsigned int> aSortedId;
  std::vector<unsigned int> aSortedMc;
  dfm2::GetSortedMortenCode(aSortedId,aSortedMc,
                            aXYZ,minmax_xyz);
  for(int ini=0;ini<aSortedMc.size()-1;++ini){
    const std::pair<int,int> range = dfm2::determineRange(aSortedMc.data(), aSortedMc.size()-1, ini);
    int isplit = dfm2::findSplit(aSortedMc.data(), range.first, range.second);
    const std::pair<int,int> rangeA = dfm2::determineRange(aSortedMc.data(), aSortedMc.size()-1, isplit);
    const std::pair<int,int> rangeB = dfm2::determineRange(aSortedMc.data(), aSortedMc.size()-1, isplit+1);
    assert( range.first == rangeA.first );
    assert( range.second == rangeB.second );
    {
      const int last1 = ( isplit == range.first ) ? isplit : rangeA.second;
      const int first1 = ( isplit+1 == range.second ) ? isplit+1 : rangeB.first;
      assert( last1+1 == first1 );
    }
  }
  // ---------------
  std::vector<dfm2::CNodeBVH2> aNodeBVH;
  dfm2::BVH_TreeTopology_Morton(aNodeBVH,
                                aSortedId,aSortedMc);
  std::vector<int> aFlg(N,0);
  mark_child(aFlg, 0, aNodeBVH);
  for(int i=0;i<N;++i){
    EXPECT_EQ(aFlg[i],1);
  }
}
