/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h" // need to be defined in the beginning

#include "delfem2/srchuni_v3.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/vec3.h"
#include "delfem2/srhbv3sphere.h"
#include "delfem2/srhbv3aabb.h"
#include "delfem2/bvh.h"
#include "delfem2/primitive.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include <random>

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ---------------------------------------------

TEST(bvh,inclusion_sphere)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
                            0.2, 0.3, 0.4);
  }
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.03);
  {
    EXPECT_EQ(bvh.aNodeBVH.size(), bvh.aBB_BVH.size());
    for (unsigned int ibb = 0; ibb < bvh.aBB_BVH.size(); ++ibb) {
      unsigned int iroot = bvh.aNodeBVH[ibb].iparent;
      if (iroot == UINT_MAX) { continue; }
      EXPECT_TRUE(bvh.aNodeBVH[iroot].ichild[0] == ibb || bvh.aNodeBVH[iroot].ichild[1] == ibb);
      const dfm2::CBV3d_Sphere &aabbp = bvh.aBB_BVH[iroot];
      const dfm2::CBV3d_Sphere &aabbc = bvh.aBB_BVH[ibb];
      EXPECT_TRUE(aabbp.IsActive());
      EXPECT_TRUE(aabbc.IsActive());
      EXPECT_TRUE(aabbp.IsInclude(aabbc,1.0e-7)); // parent includes children
    }
  }
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> dist_m1p1(-1.0, 1.0);
  for(int itr=0;itr<1000;++itr){
    dfm2::CVec3d p0(dist_m1p1(rng), dist_m1p1(rng), dist_m1p1(rng));
    for(unsigned int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3d_Sphere& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.isInclude_Point(p0.x(), p0.y(), p0.z());
      if( !is_intersect && node.ichild[1] != -1 ){ // branch
        const unsigned int ichild0 = node.ichild[0];
        const unsigned int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].isInclude_Point(p0.x(), p0.y(), p0.z()) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].isInclude_Point(p0.x(), p0.y(), p0.z()) );
      }
    }
  }
}

TEST(bvh,inclusion_aabb)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
                            0.2, 0.3, 0.4);
  }
  //  std::cout << "ntri: " << aTri.size()/3 << std::endl;
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_AABB, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.03);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> udist(-1.0, 1.0);
  for(int itr=0;itr<1000;++itr){
    dfm2::CVec3d p0(udist(rng), udist(rng), udist(rng));
    for(int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3d_AABB& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.isInclude_Point(p0.x(), p0.y(), p0.z());
      if( !is_intersect && node.ichild[1] != -1 ){ // branch
        const unsigned int ichild0 = node.ichild[0];
        const unsigned int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].isInclude_Point(p0.x(), p0.y(), p0.z()) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].isInclude_Point(p0.x(), p0.y(), p0.z()) );
      }
    }
  }
}

TEST(bvh,nearestinc_sphere)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
                            0.2, 0.3, 0.4);
  }
//  std::cout << "ntri: " << aTri.size()/3 << std::endl;
  std::vector<double> aNorm(aXYZ.size());
  dfm2::Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.03);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> udist(-5.0, 5.0);
  for(int itr=0;itr<1000;++itr){
    dfm2::CVec3d p0(udist(rng), udist(rng), udist(rng));
    {
      p0.SetNormalizedVector();
      if( itr % 2 == 0 ){ p0 *= 1.02; } // outside included in bvh
      else{               p0 *= 0.98; } // inside in included in bvh
    }
    dfm2::CPtElm2<double> pes1;
    double dist1 = bvh.Nearest_Point_IncludedInBVH(pes1,p0,0.1,
                                                   aXYZ.data(), aXYZ.size()/3,
                                                   aTri.data(), aTri.size()/3);
    EXPECT_LE( dist1, 0.1 );
    EXPECT_GE( dist1, 0.0 );
    EXPECT_TRUE( pes1.Check(aXYZ, aTri,1.0e-10) );
    dfm2::CVec3d q1 = pes1.Pos_Tri(aXYZ, aTri);
    {
      dfm2::CPtElm2d pes0 = Nearest_Point_MeshTri3D(p0, aXYZ, aTri);
      dfm2::CVec3d q0 = pes0.Pos_Tri(aXYZ, aTri);
      EXPECT_LT(Distance(q0,q1),1.0e-10);
    }
    dfm2::CVec3d n0 = pes1.UNorm_Tri(aXYZ, aTri, aNorm);
    EXPECT_EQ( n0*(p0-q1)>0, itr%2==0);
    // ---------------------
    {
      dfm2::CPtElm2d pes2;
      double dist_tri = -1, dist_bv = 0.1;
      dfm2::BVH_NearestPoint_IncludedInBVH_MeshTri3D(
          dist_tri, dist_bv, pes2,
          p0.x(), p0.y(), p0.z(), 0.1,
          aXYZ.data(), aXYZ.size()/3,
          aTri.data(), aXYZ.size()/3,
          bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
      dfm2::CVec3d q2 = pes2.Pos_Tri(aXYZ, aTri);
      EXPECT_LT(Distance(q2,q1),1.0e-10);
    }
  }
}


TEST(bvh,nearest_range) // find global nearest from range
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri, 1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
                            0.2, 0.3, 0.4);
  }
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.0);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> udist(-5.0, 5.0);
  for(int itr=0;itr<1000;++itr){
    const dfm2::CVec3d p0(udist(rng), udist(rng), udist(rng)); // random points
    double dist_min=+1, dist_max = -1;
    dfm2::BVH_Range_DistToNearestPoint(
        dist_min, dist_max,
        p0.data(),
        bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
    {
      bool is_max = false;
      for(int it=0;it<aTri.size()/3;++it){
        dfm2::CBV3d_Sphere bb_tri;
        for(int inoel=0;inoel<3;++inoel){
          const int ino0 = aTri[it*3+inoel];
          bb_tri.AddPoint(aXYZ.data()+ino0*3, 0.0);
        }
        double min0, max0;
        bb_tri.Range_DistToPoint(
            min0, max0,
            p0.x(), p0.y(), p0.z());
        EXPECT_GE( max0, dist_max );
        EXPECT_GE( min0, dist_min );
        if( max0 < dist_max+1.0e-10 ){ is_max = true; }
      }
      EXPECT_TRUE( is_max );
    }
    std::vector<unsigned int> aIndElem;
    dfm2::BVH_GetIndElem_Predicate(
        aIndElem,
        dfm2::CIsBV_InsideRange<dfm2::CBV3d_Sphere>(p0.p,dist_min,dist_max),
        bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
    EXPECT_GT(aIndElem.size(), 0);
    std::vector<int> aFlg(aTri.size()/3,0);
    for(auto itri0 : aIndElem){ aFlg[itri0] = 1; }
    for(unsigned int itri=0;itri<aTri.size()/3;++itri){
      dfm2::CBV3d_Sphere bb_tri;
      for(int inoel=0;inoel<3;++inoel){
        const int ino0 = aTri[itri*3+inoel];
        bb_tri.AddPoint(aXYZ.data()+ino0*3, 0.0);
      }
      double min0, max0;
      bb_tri.Range_DistToPoint(min0, max0, p0.x(), p0.y(), p0.z());
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

TEST(bvh,nearest_point) // find global nearest directry
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> randomDist(-5,5);
  ///
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri,
        1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
        0.2, 0.3, 0.4);
  }
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.0);
  for(int itr=0;itr<1000;++itr){
    const dfm2::CVec3d p0(
      randomDist(randomEng),
      randomDist(randomEng),
      randomDist(randomEng) );
    dfm2::CPtElm2<double> pes1 = bvh.NearestPoint_Global(p0,aXYZ,aTri);
    EXPECT_TRUE( pes1.Check(aXYZ, aTri,1.0e-10) );
    dfm2::CVec3d q1 = pes1.Pos_Tri(aXYZ, aTri);
    {
      dfm2::CPtElm2<double> pes0 = Nearest_Point_MeshTri3D(p0, aXYZ, aTri);
      dfm2::CVec3d q0 = pes0.Pos_Tri(aXYZ, aTri);
      EXPECT_LT(Distance(q0,q1),1.0e-10);
    }
  }
}


TEST(bvh,sdf) // find global nearest directry
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> randomDist(-1.5,1.5);
  // -----
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri, 
        1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
        0.2, 0.3, 0.4);
  }
  std::vector<double> aNorm(aXYZ.size());
  dfm2::Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           0.0);
  for(int itr=0;itr<1000;++itr){
    const dfm2::CVec3d p0(
      randomDist(randomEng),
      randomDist(randomEng),
      randomDist(randomEng) );
    if( (p0.Length()-1.0)<1.0e-3 ){ continue; }
    dfm2::CVec3d n0;
    double sdf = bvh.SignedDistanceFunction(n0,
                                            p0, aXYZ, aTri, aNorm);
    EXPECT_NEAR(1-p0.Length(), sdf, 1.0e-2);
    EXPECT_NEAR(n0*p0.Normalize(), 1.0, 1.0e-2 );
  }
}


TEST(bvh,lineintersection)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_m2p2(-2.,2.);
  // -----
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri,
        1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
        0.2, 0.3, 0.4);
  }
  std::vector<double> aNorm(aXYZ.size());
  dfm2::Normal_MeshTri3D(aNorm.data(),
      aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
           aTri.data(), aTri.size()/3,
           1.0e-5);
  for(int itr=0;itr<100;++itr){
    const dfm2::CVec3d s0 = dfm2::CVec3d::Random(dist_m2p2,randomDevice);
    const dfm2::CVec3d d0 = dfm2::CVec3d::Random(dist_m2p2,randomDevice).Normalize();
    for(int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3d_Sphere& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.IsIntersectLine(s0.p,d0.p);
      if( !is_intersect && node.ichild[1] != UINT_MAX ){ // branch
        const unsigned int ichild0 = node.ichild[0];
        const unsigned int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].IsIntersectLine(s0.p,d0.p) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].IsIntersectLine(s0.p,d0.p) );
      }
    }
    std::vector<unsigned int> aIndElem;
    dfm2::BVH_GetIndElem_Predicate(aIndElem,
        dfm2::CIsBV_IntersectLine<dfm2::CBV3d_Sphere>(s0.p,d0.p),
        bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
    std::vector<int> aFlg(aTri.size()/3,0);
    for(int itri0 : aIndElem){
      aFlg[itri0] = 1;
    }
    for(unsigned int itri=0;itri<aTri.size()/3;++itri){
      dfm2::CBV3d_Sphere bb;
      for(int inoel=0;inoel<3;++inoel){
        const int ino0 = aTri[itri*3+inoel];
        bb.AddPoint(aXYZ.data()+ino0*3, 1.0e-5);
      }
      bool res = bb.IsIntersectLine(s0.p, d0.p);
      EXPECT_EQ( res , aFlg[itri] == 1 );
    }
  }
}


TEST(bvh,rayintersection)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_m2p2(-2.,2.);
  // ----
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  { // make a unit sphere
    dfm2::MeshTri3D_Sphere(aXYZ, aTri,
                              1.0, 64, 32);
    dfm2::Rotate_Points3(aXYZ,
                            0.2, 0.3, 0.4);
  }
  std::vector<double> aNorm(aXYZ.size());
  dfm2::Normal_MeshTri3D(aNorm.data(),
      aXYZ.data(), aXYZ.size()/3,
      aTri.data(), aTri.size()/3);
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;
  bvh.Init(aXYZ.data(), aXYZ.size()/3,
      aTri.data(), aTri.size()/3,
      1.0e-5);
  for(int itr=0;itr<100;++itr){
    const dfm2::CVec3d s0 = dfm2::CVec3d::Random(dist_m2p2,randomDevice);
    const dfm2::CVec3d d0 = dfm2::CVec3d::Random(dist_m2p2,randomDevice).Normalize();
    for(unsigned int ibvh=0;ibvh<bvh.aNodeBVH.size();++ibvh){
      const dfm2::CBV3d_Sphere& bv = bvh.aBB_BVH[ibvh];
      const dfm2::CNodeBVH2& node = bvh.aNodeBVH[ibvh];
      bool is_intersect = bv.IsIntersectRay(s0.p,d0.p);
      if( !is_intersect && node.ichild[1] != UINT_MAX ){ // branch
        const unsigned int ichild0 = node.ichild[0];
        const unsigned int ichild1 = node.ichild[1];
        EXPECT_FALSE( bvh.aBB_BVH[ichild0].IsIntersectRay(s0.p,d0.p) );
        EXPECT_FALSE( bvh.aBB_BVH[ichild1].IsIntersectRay(s0.p,d0.p) );
      }
    }
    std::vector<unsigned int> aIndElem;
    dfm2::BVH_GetIndElem_Predicate(aIndElem,
        dfm2::CIsBV_IntersectRay<dfm2::CBV3d_Sphere>(s0.p, d0.p),
        bvh.iroot_bvh, bvh.aNodeBVH, bvh.aBB_BVH);
    std::vector<int> aFlg(aTri.size()/3,0);
    for(int itri0 : aIndElem){
      aFlg[itri0] = 1;
    }
    for(unsigned int itri=0;itri<aTri.size()/3;++itri){
      dfm2::CBV3d_Sphere bb;
      for(int inoel=0;inoel<3;++inoel){
        const int ino0 = aTri[itri*3+inoel];
        bb.AddPoint(aXYZ.data()+ino0*3, 1.0e-5);
      }
      bool res = bb.IsIntersectRay(s0.p, d0.p);
      EXPECT_EQ( res , aFlg[itri] == 1 );
    }
    {
      std::map<double,dfm2::CPtElm2d > mapDepthPES0;
      IntersectionRay_MeshTri3(mapDepthPES0,
          s0, d0, aTri, aXYZ,
          0.0);
      std::map<double,dfm2::CPtElm2d > mapDepthPES1;
      IntersectionRay_MeshTri3DPart(mapDepthPES1,
          s0, d0,
          aTri, aXYZ, aIndElem,
          0.0);
      EXPECT_EQ(mapDepthPES0.size(),mapDepthPES1.size());
      int N = mapDepthPES0.size();
      auto itr0 = mapDepthPES0.begin();
      auto itr1 = mapDepthPES0.begin();
      for(int i=0;i<N;++i){
        EXPECT_FLOAT_EQ(itr0->first,itr1->first);
        dfm2::CPtElm2d pes0 = itr0->second;
        dfm2::CPtElm2d pes1 = itr1->second;
        const dfm2::CVec3d q0 = pes0.Pos_Tri(aXYZ, aTri);
        const dfm2::CVec3d q1 = pes1.Pos_Tri(aXYZ, aTri);
        EXPECT_NEAR(Distance(q0,q1), 0.0, 1.0e-10);
      }
    }
  }
}

void mark_child(std::vector<int>& aFlgBranch,
                std::vector<int>& aFlgLeaf,
                std::vector<int>& aFlgID,
                unsigned int nID,
                unsigned int inode0,
                const std::vector<dfm2::CNodeBVH2>& aNode)
{
  EXPECT_TRUE( inode0 < aNode.size() );
  if( aNode[inode0].ichild[1] == -1 ){ // leaf
    EXPECT_TRUE(inode0>=nID-1 && inode0<nID*2-1);
    aFlgLeaf[inode0-(nID-1)] += 1;
    const unsigned int in0 = aNode[inode0].ichild[0];
    EXPECT_TRUE( in0 < aFlgID.size() );
    aFlgID[in0] += 1;
    return;
  }
  EXPECT_TRUE(inode0<nID-1);
  aFlgBranch[inode0] += 1;
  const unsigned int in0 = aNode[inode0].ichild[0];
  const unsigned int in1 = aNode[inode0].ichild[1];
  mark_child(aFlgBranch, aFlgLeaf, aFlgID, nID, in0, aNode);
  mark_child(aFlgBranch, aFlgLeaf, aFlgID, nID, in1, aNode);
}

TEST(bvh,clz) {
  const unsigned int n0 = dfm2::nbits_leading_zero(0);
  EXPECT_EQ(n0,32);
}

TEST(bvh,morton_code)
{
  std::vector<double> aXYZ; // 3d points
  const double min_xyz[3] = {-1,-1,-1};
  const double max_xyz[3] = {+1,+1,+1};
  dfm2::CBV3d_AABB bb(min_xyz,max_xyz);
  {
    const unsigned int N = 10000;
    aXYZ.resize(N*3);
    std::random_device randomDevice;
    std::mt19937 randomEng(randomDevice());
    std::uniform_real_distribution<> dist_01(0.0, 1.0);
    for(int i=0;i<N;++i){
      aXYZ[i*3+0] = (bb.bbmax[0] -  bb.bbmin[0]) * dist_01(randomEng) + bb.bbmin[0];
      aXYZ[i*3+1] = (bb.bbmax[1] -  bb.bbmin[1]) * dist_01(randomEng) + bb.bbmin[1];
      aXYZ[i*3+2] = (bb.bbmax[2] -  bb.bbmin[2]) * dist_01(randomEng) + bb.bbmin[2];
    }
    for(int iip=0;iip<3;++iip){ // hash collision
      const auto ip = static_cast<unsigned int>(N * dist_01(randomEng));
      assert( ip < N );
      const double x0 = aXYZ[ip*3+0];
      const double y0 = aXYZ[ip*3+1];
      const double z0 = aXYZ[ip*3+2];
      for(int itr=0;itr<2;itr++){
        aXYZ.insert(aXYZ.begin(), z0);
        aXYZ.insert(aXYZ.begin(), y0);
        aXYZ.insert(aXYZ.begin(), x0);
      }
    }
  }
  std::vector<unsigned int> aSortedId;
  std::vector<std::uint32_t> aSortedMc;
  dfm2::SortedMortenCode_Points3(aSortedId,aSortedMc,
                                 aXYZ,
                                 min_xyz,max_xyz);
  for(unsigned int ini=0;ini<aSortedMc.size()-1;++ini){
    const std::pair<int,int> range = dfm2::MortonCode_DeterminRange(aSortedMc.data(), aSortedMc.size(), ini);
    unsigned int isplit = dfm2::MortonCode_FindSplit(aSortedMc.data(), range.first, range.second);
    const std::pair<int,int> rangeA = dfm2::MortonCode_DeterminRange(aSortedMc.data(), aSortedMc.size(), isplit);
    const std::pair<int,int> rangeB = dfm2::MortonCode_DeterminRange(aSortedMc.data(), aSortedMc.size(), isplit+1);
    EXPECT_EQ( range.first, rangeA.first );
    EXPECT_EQ( range.second, rangeB.second );
    {
      const int last1 = ( isplit == range.first ) ? isplit : rangeA.second;
      const int first1 = ( isplit+1 == range.second ) ? isplit+1 : rangeB.first;
      EXPECT_EQ( last1+1, first1 );
    }
  }
  // ---------------
  std::vector<dfm2::CNodeBVH2> aNodeBVH;
  dfm2::BVHTopology_Morton(aNodeBVH,
                           aSortedId,aSortedMc);
  {
    const unsigned int N = aXYZ.size()/3;
    std::vector<int> aFlgBranch(N-1,0);
    std::vector<int> aFlgLeaf(N,0);
    std::vector<int> aFlgID(N,0);
    mark_child(aFlgBranch,aFlgLeaf,aFlgID, N,
               0,aNodeBVH);
    for(unsigned int i=0;i<N;++i){
      EXPECT_EQ(aFlgLeaf[i],1);
      EXPECT_EQ(aFlgID[i],1);
    }
    for(int i=0;i<N-1;++i){
      EXPECT_EQ(aFlgBranch[i],1);
    }
  }
  // ----------------
  std::vector<dfm2::CBV3_Sphere<double>> aAABB;
  dfm2::CLeafVolumeMaker_Point<dfm2::CBV3_Sphere<double>,double> lvm(
      aXYZ.data(), aXYZ.size()/3);
  dfm2::BVH_BuildBVHGeometry(aAABB,
      0, aNodeBVH, lvm);
  for(int itr=0;itr<100;++itr){
    const double cur_time = itr*0.07 + 0.02;
    const  double p0[3] = {
      1.5*(bb.bbmax[0]-bb.bbmin[0])*sin(cur_time*1)-(bb.bbmax[0]+bb.bbmin[0])*0.5,
      1.5*(bb.bbmax[1]-bb.bbmin[1])*sin(cur_time*2)-(bb.bbmax[1]+bb.bbmin[1])*0.5,
      1.5*(bb.bbmax[2]-bb.bbmin[2])*sin(cur_time*3)-(bb.bbmax[2]+bb.bbmin[2])*0.5 };
    double dist = -1;
    unsigned int ip_nearest = 0;
    dfm2::BVH_IndPoint_NearestPoint(
        ip_nearest, dist, p0, 0,
        aNodeBVH,aAABB);
    for(unsigned int ip=0;ip<aXYZ.size()/3;++ip){
      double dist0 = dfm2::Distance3(p0, aXYZ.data()+ip*3);
      EXPECT_GE(dist0,dist);
    }
  }
}
