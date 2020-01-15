#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include "gtest/gtest.h"

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

