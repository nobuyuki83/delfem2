#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include "gtest/gtest.h"

#include "delfem2/mshmisc.h"

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

    float bb3A[6];
    dfm2::BB3_Points3(bb3A,
                      aXYZ.data(), aXYZ.size()/3);

    float bb3B[6];
    dfm2::cuda::cuda_MinMax_Point3D(bb3B,
                                    aXYZ.data(), aXYZ.size() / 3);

    for(int i=0;i<6;++i) {
      EXPECT_FLOAT_EQ(bb3A[0], bb3B[0]);
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
  std::mt19937 engine(0);
  std::uniform_int_distribution<unsigned int> dist0(1, 400);
  std::uniform_real_distribution<float> dist1(-1.0, 1.0);

  for (int itr = 0; itr < 10; ++itr) {
    const unsigned int n = dist0(engine);
    std::vector<float> A(n*n), B(n*n);
    for(int i=0;i<n*n;++i){
      A[i] = dist1(engine);
      B[i] = dist1(engine);
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









