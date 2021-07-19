/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <vector>
#include <random>
#include <bitset>
#include "gtest/gtest.h"
#include "delfem2/points.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/vec3.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"

#include "delfem2/cuda/cu_matvec.h"
#include "delfem2/cuda/cu_bvh.h"

namespace dfm2 = delfem2;

// ----------------------------------------------


TEST(matvec,vecscale)
{
  std::uniform_int_distribution<unsigned int> distUI(1,2000);
  std::uniform_real_distribution <float> distF(-1.0, +1.0);
  std::random_device randomDevice;
  std::mt19937 engin(randomDevice());

  for(int itr=0;itr<1000;++itr) {
    const unsigned int n = distUI(engin);
    std::vector<float> in(n), out(n);
    for (int i = 0; i < n; i++) in[i] = distF(engin);

    const float scale = 2.0;
    dfm2::cuda::cuda_VecScale(out.data(), in.data(), scale, n);

    for (int i = 0; i < n; i++) {
      EXPECT_FLOAT_EQ(in[i]*scale,out[i]);
    }
  }
}


TEST(matvec,dot)
{
  std::random_device randomDevice;
  std::mt19937 engine(randomDevice());
  std::uniform_int_distribution<unsigned int> dist0(1, 20000);
  std::uniform_real_distribution<float> dist1(-1.0, 1.0);

  for(int itr=0;itr<500;++itr) {
    const unsigned int n = dist0(engine);
    std::vector<float> A(n), B(n);
    for (int i = 0; i < n; ++i) {
      A[i] = dist1(engine);
      B[i] = dist1(engine);
    }

    float dot0 = 0.0;
    for (int i = 0; i < n; ++i) { dot0 += A[i] * B[i]; }

    float dot1 = dfm2::cuda::cuda_Dot(A.data(), B.data(), n);

    EXPECT_NEAR(dot0,dot1,1.0e-3);
  }
}

TEST(matvec,matmat) {
  std::random_device randomDevice;
  std::mt19937 engin(randomDevice());
  std::uniform_int_distribution<unsigned int> dist0(1, 400);
  std::uniform_real_distribution<float> dist1(-1.0, 1.0);
  // ------------------------------
  for (int itr = 0; itr < 10; ++itr) {
    const unsigned int n = dist0(engin);
    std::vector<float> A(n*n), B(n*n);
    for(int i=0;i<n*n;++i){
      A[i] = dist1(engin);
      B[i] = dist1(engin);
    }
    // ---------------------------
    std::vector<float> C0(n*n);
    dfm2::cuda::cuda_MatMat(C0.data(),
                            A.data(), B.data(), n);
    // ---------------------------
    std::vector<float> C1(n*n);
    for(int i=0;i<n;++i) {
      for(int j=0;j<n;++j){
        float tmp = 0.0;
        for(int k=0;k<n;++k){ tmp += A[n*i+k] * B[n*k+j]; }
        C1[n*i+j] = tmp;
      }
    }
    // ---------------------------
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        EXPECT_NEAR(C0[i * n + j], C1[i * n + j], 1.0e-4);
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------






TEST(bvh,minmax_po3d)
{
  std::random_device randomDevice;
  std::mt19937 engine(randomDevice());
  std::uniform_int_distribution<unsigned int> dist0(1, 20000);
  std::uniform_real_distribution<float> dist1(-1.0, 1.0);
  std::uniform_real_distribution<float> dist2(-2.0, 2.0);

  for(int itr=0;itr<300;++itr) {
    std::vector<float> aXYZ;
    {
      const unsigned int np = dist0(engine);
      aXYZ.resize(np * 3);
      for (int ip = 0; ip < np; ++ip) {
        aXYZ[ip * 3 + 0] = dist1(engine) * 1.0f + dist2(engine);
        aXYZ[ip * 3 + 1] = dist1(engine) * 0.5f + dist2(engine);
        aXYZ[ip * 3 + 2] = dist1(engine) * 0.3f + dist2(engine);
      }
    }

    float min3A[3], max3A[3];
    dfm2::BoundingBox3_Points3(min3A,max3A,
        aXYZ.data(), aXYZ.size()/3);

    float min3B[3], max3B[3];
    dfm2::cuda::cuda_Min3Max3_Points3F(min3B,max3B,
        aXYZ.data(), aXYZ.size() / 3);

    EXPECT_FLOAT_EQ(min3A[0], min3B[0]);
    EXPECT_FLOAT_EQ(min3A[1], min3B[1]);
    EXPECT_FLOAT_EQ(min3A[2], min3B[2]);
    EXPECT_FLOAT_EQ(max3A[0], max3B[0]);
    EXPECT_FLOAT_EQ(max3A[1], max3B[1]);
    EXPECT_FLOAT_EQ(max3A[2], max3B[2]);
  }
}

TEST(bvh,meshtri3d_centrad)
{
  std::mt19937 engin(0);
  std::uniform_int_distribution<unsigned int> dist0(3, 100);
  // ----------------------
  std::vector<float> aXYZ;
  std::vector<unsigned int> aTri;
  for(int itr=0;itr<100;++itr) {
    unsigned int nr = dist0(engin);
    unsigned int nl = dist0(engin);
    dfm2::MeshTri3_Torus(aXYZ, aTri,
                         0.5f, 0.2f, nr, nl );
    const unsigned int nTri = aTri.size() / 3;
    // -----------------------------------------------------
    std::vector<float> aXYZ_c0;
    float max_rad0 = -1;
    max_rad0 = dfm2::CentsMaxRad_MeshTri3(aXYZ_c0,
                                          aXYZ, aTri);
    // ------------------
    std::vector<float> aXYZ_c1(nTri * 3);
    float max_rad1;
    dfm2::cuda::cuda_CentsMaxRad_MeshTri3F(
        aXYZ_c1.data(), &max_rad1,
        aXYZ.data(), aXYZ.size() / 3,
        aTri.data(), nTri);
    // ------------------
    for (unsigned int itri = 0; itri < nTri; ++itri) {
      EXPECT_FLOAT_EQ(aXYZ_c0[itri * 3 + 0], aXYZ_c1[itri * 3 + 0]);
      EXPECT_FLOAT_EQ(aXYZ_c0[itri * 3 + 1], aXYZ_c1[itri * 3 + 1]);
      EXPECT_FLOAT_EQ(aXYZ_c0[itri * 3 + 2], aXYZ_c1[itri * 3 + 2]);
    }
    EXPECT_FLOAT_EQ(max_rad0, max_rad1);
  }
}

void FlagBVHMortonCode(
    std::vector<int>& aFlgBranch,
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
  aFlgBranch[inode0] += 1;
  const unsigned int in0 = aNode[inode0].ichild[0];
  const unsigned int in1 = aNode[inode0].ichild[1];
  FlagBVHMortonCode(aFlgBranch, aFlgLeaf, aFlgID, nID, in0, aNode);
  FlagBVHMortonCode(aFlgBranch, aFlgLeaf, aFlgID, nID, in1, aNode);
}

TEST(bvh,morton_code) {
  std::uniform_real_distribution<float> udist0(0.0, 1.0);
  std::uniform_int_distribution<unsigned int> udist1(0, 100000);
  std::random_device randomDevice;
  std::mt19937 rng(randomDevice());
  // -----------------------------------
  for(int itr=0;itr<10;++itr) {
    const float bbmin[3] = {0.f, 0.f, 0.f};
    const float bbmax[3] = {1.f, 1.f, 1.f};
    std::vector<float> aXYZ; // 3d points
    { // making points
      const unsigned int N = udist1(rng);
      aXYZ.resize(N * 3);
      for (int i = 0; i < N; ++i) {
        aXYZ[i * 3 + 0] = udist0(rng);
        aXYZ[i * 3 + 1] = udist0(rng);
        aXYZ[i * 3 + 2] = udist0(rng);
      }
      std::uniform_int_distribution<unsigned int> udist2(0, N-1);
      for(int iip=0;iip<3;++iip){ // hash collision
        const unsigned int ip = udist2(rng);
        assert( N >= 0 && ip < N);
        const float x0 = aXYZ[ip*3+0];
        const float y0 = aXYZ[ip*3+1];
        const float z0 = aXYZ[ip*3+2];
        for(int iduplicate=0;iduplicate<2;iduplicate++){
          aXYZ.insert(aXYZ.begin(), z0);
          aXYZ.insert(aXYZ.begin(), y0);
          aXYZ.insert(aXYZ.begin(), x0);
        }
      }
    }
    // ---------------------------------------------
    const unsigned int N = aXYZ.size()/3;
    std::vector<unsigned int> aSortedId(N);
    std::vector<std::uint32_t> aSortedMc(N);
    dfm2::cuda::cuda_MortonCode_Points3FSorted(aSortedId.data(), aSortedMc.data(),
                                               aXYZ.data(), aXYZ.size() / 3,
                                               bbmin, bbmax);
    { // check sorted morton code
      for(unsigned int imc=1;imc<aSortedMc.size();++imc){
        std::uint32_t mc0 = aSortedMc[imc-1];
        std::uint32_t mc1 = aSortedMc[imc+0];
        EXPECT_LE( mc0, mc1 );
      }
      for(unsigned int imc=0;imc<aSortedMc.size();++imc){
        std::uint32_t mc0 = aSortedMc[imc];
        unsigned int ip = aSortedId[imc];
        float x0 = aXYZ[ip*3+0];
        float y0 = aXYZ[ip*3+1];
        float z0 = aXYZ[ip*3+2];
        float x1 = (x0-bbmin[0])/(bbmax[0]-bbmin[0]);
        float y1 = (y0-bbmin[1])/(bbmax[1]-bbmin[1]);
        float z1 = (z0-bbmin[2])/(bbmax[2]-bbmin[2]);
        std::uint32_t mc1 = dfm2::MortonCode(x1,y1,z1);
        EXPECT_EQ( mc0, mc1 );
      }
    }

    std::vector<dfm2::CNodeBVH2> aNodeBVH(N*2-1);
    dfm2::cuda::cuda_MortonCode_BVHTopology(aNodeBVH.data(),
                                            aSortedId.data(), aSortedMc.data(), N);
    { // check topology
      std::vector<dfm2::CNodeBVH2> aNodeBVH1(N*2-1);
      dfm2::BVHTopology_Morton(aNodeBVH1,
          aSortedId, aSortedMc);
      EXPECT_EQ(aNodeBVH.size(), aNodeBVH1.size());
      for(unsigned int ibb=0;ibb<aNodeBVH.size();++ibb) {
        EXPECT_EQ( aNodeBVH[ibb].iparent, aNodeBVH1[ibb].iparent );
        EXPECT_EQ( aNodeBVH[ibb].ichild[0], aNodeBVH1[ibb].ichild[0] );
        EXPECT_EQ( aNodeBVH[ibb].ichild[1], aNodeBVH1[ibb].ichild[1] );
      }
      for(unsigned int ibb=0;ibb<aNodeBVH.size();++ibb) {
        const unsigned int iroot = aNodeBVH[ibb].iparent;
        if( iroot != UINT_MAX ) {
          EXPECT_TRUE(aNodeBVH[iroot].ichild[0] == ibb || aNodeBVH[iroot].ichild[1] == ibb);
        }
        if( aNodeBVH[ibb].ichild[1] == UINT_MAX ){ // leaf
          const unsigned int itri = aNodeBVH[ibb].ichild[0];
          EXPECT_GE(itri,0);
          EXPECT_LT(itri, N);
        }
      }
      std::vector<int> aFlgBranch(aXYZ.size()/3-1,0);
      std::vector<int> aFlgLeaf(aXYZ.size()/3,0);
      std::vector<int> aFlgID(aXYZ.size()/3,0);
      FlagBVHMortonCode(aFlgBranch, aFlgLeaf, aFlgID, N, 0, aNodeBVH);
      for(unsigned int i=0;i<N;++i){
        EXPECT_EQ(aFlgID[i],1);
        EXPECT_EQ(aFlgLeaf[i],1);
      }
      for(unsigned int i=0;i<N-1;++i) {
        EXPECT_EQ(aFlgBranch[i],1);
      }
    } // end checking bvh topology
    // ------------------------------------------------------
    std::vector<dfm2::CBV3_Sphere<float>> aAABB1;
    dfm2::BVH_BuildBVHGeometry(aAABB1, 0, aNodeBVH,
                               dfm2::CLeafVolumeMaker_Point<dfm2::CBV3_Sphere<float>,float>( aXYZ.data(), aXYZ.size()/3) );

    { // compare nearest points
      const unsigned int npt = 100;
      std::vector<float> aXYZ_Test(npt*3);
      for(unsigned int ipt=0;ipt<npt;++ipt) {
        const float cur_time = (float)ipt * 0.07f + 0.02f;
        aXYZ_Test[ipt * 3 + 0] = 1.5f * (bbmax[0] - bbmin[0]) * sin(cur_time * 1) - (bbmax[0] + bbmin[0]) * 0.5f;
        aXYZ_Test[ipt * 3 + 1] = 1.5f * (bbmax[1] - bbmin[1]) * sin(cur_time * 2) - (bbmax[1] + bbmin[1]) * 0.5f;
        aXYZ_Test[ipt * 3 + 2] = 1.5f * (bbmax[2] - bbmin[2]) * sin(cur_time * 3) - (bbmax[2] + bbmin[2]) * 0.5f;
      }
      std::vector<unsigned int> aId(npt);
      dfm2::cuda::cuda_BVH_NearestPoint(aId.data(),
          aXYZ_Test.data(), aXYZ_Test.size()/3,
          aNodeBVH.data(), aNodeBVH.size(),
          aAABB1.data());
      for(unsigned int ipt=0;ipt<npt;++ipt){
        const float* p0 = aXYZ_Test.data() + ipt*3;
        unsigned int ip_nearest = aId[ipt];
        float dist_min = dfm2::Distance3(p0, aXYZ.data()+ip_nearest*3);
        for(unsigned int ip=0;ip<aXYZ.size()/3;++ip) {
          double dist0 = dfm2::Distance3(p0, aXYZ.data() + ip * 3);
          EXPECT_GE(dist0, dist_min);
        }
      }
    } // compare nearest point

  }
}


TEST(bvh,aabb_tri)
{
  std::random_device randomDevice;
  std::mt19937 engin(randomDevice());
  std::uniform_int_distribution<unsigned int> dist0(3, 30);
  std::uniform_real_distribution<float> dist1(0, 1);
  std::vector<float> aXYZ;
  std::vector<unsigned int> aTri;
  //
  for(int itr=0;itr<1;++itr) {
    unsigned int nr = dist0(engin);
    unsigned int nl = dist0(engin);
    dfm2::MeshTri3_Torus(aXYZ, aTri,
                         0.5f, 0.2f, nr, nl);
    for (int ip = 0; ip < aXYZ.size() / 3; ++ip) {
      aXYZ[ip * 3 + 0] += 0.05f * dist1(engin);
      aXYZ[ip * 3 + 1] += 0.05f * dist1(engin);
      aXYZ[ip * 3 + 2] += 0.05f * dist1(engin);
    }
    const unsigned int nTri = aTri.size() / 3;
    // -----------------------------------------------------
    std::vector<float> aXYZ_c(nTri * 3);
    float max_rad1;
    dfm2::cuda::cuda_CentsMaxRad_MeshTri3F(
        aXYZ_c.data(), &max_rad1,
        aXYZ.data(), aXYZ.size() / 3,
        aTri.data(), nTri);
    float bbmin3[3], bbmax3[3];
    dfm2::cuda::cuda_Min3Max3_Points3F(bbmin3, bbmax3,
                                       aXYZ_c.data(), aXYZ_c.size() / 3);
    std::vector<dfm2::CNodeBVH2> aNodeBVH(nTri * 2 - 1);
    {
      std::vector<unsigned int> aSortedId(nTri);
      std::vector<std::uint32_t> aSortedMc(nTri);
      dfm2::cuda::cuda_MortonCode_Points3FSorted(aSortedId.data(), aSortedMc.data(),
                                                 aXYZ_c.data(), aXYZ_c.size() / 3,
                                                 bbmin3, bbmax3);
      dfm2::cuda::cuda_MortonCode_BVHTopology(aNodeBVH.data(),
                                              aSortedId.data(), aSortedMc.data(), nTri);
      // -------------------------------------
      { // check topology against CPU
        std::vector<dfm2::CNodeBVH2> aNodeBVH1(nTri * 2 - 1);
        dfm2::BVHTopology_Morton(aNodeBVH1,
                                      aSortedId, aSortedMc);
        EXPECT_EQ(aNodeBVH.size(), aNodeBVH1.size());
        for (int ibb = 0; ibb < aNodeBVH.size(); ++ibb) {
          EXPECT_EQ(aNodeBVH[ibb].iparent, aNodeBVH1[ibb].iparent);
          EXPECT_EQ(aNodeBVH[ibb].ichild[0], aNodeBVH1[ibb].ichild[0]);
          EXPECT_EQ(aNodeBVH[ibb].ichild[1], aNodeBVH1[ibb].ichild[1]);
        }
      }
    }
    { // check if visiting only once
      const unsigned int N = nTri;
      std::vector<int> aFlgBranch(N - 1, 0);
      std::vector<int> aFlgLeaf(N, 0);
      std::vector<int> aFlgID(N, 0);
      FlagBVHMortonCode(aFlgBranch, aFlgLeaf, aFlgID, N, 0, aNodeBVH);
      for (int i = 0; i < N; ++i) {
        EXPECT_EQ(aFlgID[i], 1);
        EXPECT_EQ(aFlgLeaf[i], 1);
      }
      for (int i = 0; i < N - 1; ++i) {
        EXPECT_EQ(aFlgBranch[i], 1);
      }
    }
    { // AABB as bounding volume
      std::vector<dfm2::CBV3f_AABB> aAABB(nTri * 2 - 1);
      dfm2::cuda::cuda_BVHGeometry_AABB3f(aAABB.data(),
                                          aNodeBVH.data(),
                                          aXYZ.data(), aXYZ.size() / 3,
                                          aTri.data(), nTri);
      { // check gpu computed AABB is same as cpu computed AABB
        std::vector<dfm2::CBV3f_AABB> aAABB1(nTri * 2 - 1);
        dfm2::BVH_BuildBVHGeometry(
            aAABB1,
            0, aNodeBVH,
            dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3f_AABB,float>(
                0.f,
                aXYZ.data(), aXYZ.size() / 3,
                aTri.data(), aTri.size() / 3, 3) );
        for (int ibb = 0; ibb < aAABB.size(); ++ibb) {
          EXPECT_FLOAT_EQ(aAABB[ibb].bbmin[0], aAABB1[ibb].bbmin[0]);
          EXPECT_FLOAT_EQ(aAABB[ibb].bbmin[1], aAABB1[ibb].bbmin[1]);
          EXPECT_FLOAT_EQ(aAABB[ibb].bbmin[2], aAABB1[ibb].bbmin[2]);
          EXPECT_FLOAT_EQ(aAABB[ibb].bbmax[0], aAABB1[ibb].bbmax[0]);
          EXPECT_FLOAT_EQ(aAABB[ibb].bbmax[1], aAABB1[ibb].bbmax[1]);
          EXPECT_FLOAT_EQ(aAABB[ibb].bbmax[2], aAABB1[ibb].bbmax[2]);
        }
      }
      // ---------------------------------------------
      EXPECT_EQ(aNodeBVH.size(), aAABB.size());
      for (int ibb = 0; ibb < aAABB.size(); ++ibb) {
        const unsigned int iroot = aNodeBVH[ibb].iparent;
        if (iroot == UINT_MAX) { continue; }
        EXPECT_TRUE(aNodeBVH[iroot].ichild[0] == ibb || aNodeBVH[iroot].ichild[1] == ibb);
        const dfm2::CBV3f_AABB &aabbp = aAABB[iroot];
        const dfm2::CBV3f_AABB &aabbc = aAABB[ibb];
        EXPECT_TRUE(aabbp.IsActive());
        EXPECT_TRUE(aabbc.IsActive());
        EXPECT_TRUE(aabbp.IsInclude_AABB3(aabbc)); // parent includes children
      }
    } // end: AABB as bounding volume
    { // begin: sphere as bounding volume
      std::vector<dfm2::CBV3f_Sphere> aSphere(nTri * 2 - 1);
      dfm2::cuda::cuda_BVHGeometry_Sphere(aSphere.data(),
                                          aNodeBVH.data(),
                                          aXYZ.data(), aXYZ.size() / 3,
                                          aTri.data(), nTri);
      {
        std::vector<dfm2::CBV3f_Sphere> aSphere1(nTri * 2 - 1);
        dfm2::BVH_BuildBVHGeometry(
            aSphere1,
            0, aNodeBVH,
            dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3f_Sphere,float>(
            0.f,
            aXYZ.data(), aXYZ.size() / 3,
            aTri.data(), aTri.size() / 3, 3) );
        for (unsigned int ibb = 0; ibb < aSphere.size(); ++ibb) {
          EXPECT_NEAR(aSphere[ibb].r, aSphere1[ibb].r, 1.0e-6);
          EXPECT_NEAR(aSphere[ibb].c[0], aSphere1[ibb].c[0], 1.0e-6);
          EXPECT_NEAR(aSphere[ibb].c[1], aSphere1[ibb].c[1], 1.0e-6);
          EXPECT_NEAR(aSphere[ibb].c[2], aSphere1[ibb].c[2], 1.0e-6);
        }
        // ---------------------------------------------
        EXPECT_EQ(aNodeBVH.size(), aSphere.size());
        for (unsigned int ibb = 0; ibb < aSphere.size(); ++ibb) {
          const unsigned int iroot = aNodeBVH[ibb].iparent;
          if (iroot == UINT_MAX) { continue; }
          EXPECT_TRUE(aNodeBVH[iroot].ichild[0] == ibb || aNodeBVH[iroot].ichild[1] == ibb);
          const dfm2::CBV3f_Sphere &aabbp = aSphere1[iroot];
          const dfm2::CBV3f_Sphere &aabbc = aSphere1[ibb];
          EXPECT_TRUE(aabbp.IsActive());
          EXPECT_TRUE(aabbc.IsActive());
          EXPECT_TRUE(aabbp.IsInclude(aabbc, 1.0e-6)); // parent includes children
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------------------------------





