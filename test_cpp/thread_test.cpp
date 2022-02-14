/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstring>
#include <random>

#include "gtest/gtest.h"
#include "delfem2/thread.h"


#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(thread, parallel_for0) {
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_int_distribution<unsigned int> dist0(0, 100);
  std::uniform_int_distribution<unsigned int> dist1(0, 5);
  for (unsigned int itr = 0; itr < 100; ++itr) {
    const unsigned int N = dist0(rdeng);
    std::vector<int> aIn(N);
    for (unsigned int i = 0; i < N; ++i) { aIn[i] = i; }
    std::vector<int> aOut(aIn.size());
    auto func0 = [&aIn, &aOut](int i) { aOut[i] = aIn[i] * aIn[i]; }; // square
    const unsigned int nthread = dist1(rdeng);
    dfm2::parallel_for(
      static_cast<unsigned int>(aIn.size()),
      func0,
      nthread);
    std::vector<int> aTrg(N);
    for (unsigned int i = 0; i < N; ++i) { aTrg[i] = i * i; }
    EXPECT_TRUE(0 == std::memcmp(aOut.data(), aTrg.data(), aIn.size() * sizeof(int)));
  }
}

TEST(thread, parallel_for1) {
  std::mt19937 rdeng(std::random_device{}());
  std::uniform_int_distribution<unsigned int> dist0(0, 100);
  std::uniform_int_distribution<unsigned int> dist1(0, 5);
  for (unsigned int itr = 0; itr < 100; ++itr) {
    const unsigned int N = dist0(rdeng);
    const unsigned int nthread = dist1(rdeng);
    std::vector<unsigned int> out0(N*N,0);
    auto func0 = [&N, &out0](int i, int j) { out0[N*j+i] = N*j+i; };
    delfem2::parallel_for(N,N,func0);
    for(unsigned int i=0;i<N*N;++i){
      EXPECT_EQ(out0[i],i);
    }
  }
}
