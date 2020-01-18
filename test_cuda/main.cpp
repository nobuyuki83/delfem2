/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include "gtest/gtest.h"

#include "delfem2/mshmisc.h"
#include "delfem2/primitive.h"
#include "delfem2/vec3.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"

#include "delfem2/cuda/cu_matvec.h"

namespace dfm2 = delfem2;

// ----------------------------------------------


TEST(matvec,vecscale)
{
  std::uniform_int_distribution<unsigned int> distUI(1,2000);
  std::uniform_real_distribution <float> distF(-1.0, +1.0);
  std::mt19937 engin(0);

  for(int itr=0;itr<1000;++itr) {
    const int n = distUI(engin);
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
  std::mt19937 engine(0);
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
  std::mt19937 engin(0);
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



TEST(matvec,minmax_po3d)
{
  std::mt19937 engine(0);
  std::uniform_int_distribution<unsigned int> dist0(1, 20000);
  std::uniform_real_distribution<float> dist1(-1.0, 1.0);
  std::uniform_real_distribution<float> dist2(-2.0, 2.0);

  for(int itr=0;itr<300;++itr) {
    std::vector<float> aXYZ;
    {
      const unsigned int np = dist0(engine);
      aXYZ.resize(np * 3);
      for (int ip = 0; ip < np; ++ip) {
        aXYZ[ip * 3 + 0] = dist1(engine) * 1.0 + dist2(engine);
        aXYZ[ip * 3 + 1] = dist1(engine) * 0.5 + dist2(engine);
        aXYZ[ip * 3 + 2] = dist1(engine) * 0.3 + dist2(engine);
      }
    }

    float min3A[3], max3A[3];
    dfm2::Min3Max3_Points3(min3A,max3A,
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

// ---------------------------------------------------------------------------------------------------------------------

TEST(matvec,meshtri3d_centrad)
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
                         0.5, 0.20, nr, nl );
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


TEST(bvh,morton_code) {
  std::vector<float> aXYZ; // 3d points
  std::uniform_real_distribution<> udist(0.0, 1.0);
  std::mt19937 rng(0);

  const unsigned int N = 10000;
  aXYZ.resize(N * 3);
  for (int i = 0; i < N; ++i) {
    aXYZ[i * 3 + 0] = udist(rng);
    aXYZ[i * 3 + 1] = udist(rng);
    aXYZ[i * 3 + 2] = udist(rng);
  }
  std::vector<uint32_t> aMC0(N);
  for(int i=0;i<N;++i) {
    aMC0[i] = dfm2::MortonCode(aXYZ[i*3+0],aXYZ[i*3+1],aXYZ[i*3+2]);
  }

  std::vector<uint32_t> aMC1(N);
  dfm2::cuda::cuda_MortonCode_Points3F(aMC1.data(),aXYZ.data(),aXYZ.size()/3);

  for(int i=0;i<N;++i){
    EXPECT_EQ(aMC0[i],aMC1[i]);
  }



}









