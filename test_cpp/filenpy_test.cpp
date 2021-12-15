/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>

#include "gtest/gtest.h"
#include "delfem2/filenpy_str.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(funcs, numpy_load_2df) {
  std::string path = std::string(PATH_INPUT_DIR) + "/numpy_array4x4_float.npy";
  int ndim0, ndim1;
  std::vector<float> aData;
  bool res = dfm2::LoadNumpy_2Dim(
    ndim0, ndim1, aData,
    path);
  EXPECT_TRUE(res);
  EXPECT_EQ(ndim0, 4);
  EXPECT_EQ(ndim1, 4);
  EXPECT_EQ(aData.size(), ndim0 * ndim1);
}

TEST(funcs, numpy_load_2dd) {
  std::string path = std::string(PATH_INPUT_DIR) + "/numpy_array4x4_double.npy";
  int ndim0, ndim1;
  std::vector<double> aData;
  bool res = dfm2::LoadNumpy_2Dim(
    ndim0, ndim1, aData,
    path);
  EXPECT_TRUE(res);
  EXPECT_EQ(ndim0, 4);
  EXPECT_EQ(ndim1, 4);
  EXPECT_EQ(aData.size(), ndim0 * ndim1);
}

TEST(funcs, numpy_load_1df) {
  std::string path = std::string(PATH_INPUT_DIR) + "/numpy_array4_float.npy";
  int ndim0;
  std::vector<float> aData;
  bool res = dfm2::LoadNumpy_1DimF(
    ndim0, aData,
    path);
  EXPECT_TRUE(res);
  EXPECT_EQ(ndim0, 4);
  EXPECT_EQ(aData.size(), ndim0);
}