//
// Created by Nobuyuki Umetani on 2021-11-19.
//

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/fem_solidhyper.h"

void MakeRandomHex(
    double aP0[8][3],
    double aU0[8][3]) {
  aP0[0][0] = -1;
  aP0[0][1] = -1;
  aP0[0][2] = -1;
  aP0[1][0] = +1;
  aP0[1][1] = -1;
  aP0[1][2] = -1;
  aP0[2][0] = +1;
  aP0[2][1] = +1;
  aP0[2][2] = -1;
  aP0[3][0] = -1;
  aP0[3][1] = +1;
  aP0[3][2] = -1;
  aP0[4][0] = -1;
  aP0[4][1] = -1;
  aP0[4][2] = +1;
  aP0[5][0] = +1;
  aP0[5][1] = -1;
  aP0[5][2] = +1;
  aP0[6][0] = +1;
  aP0[6][1] = +1;
  aP0[6][2] = +1;
  aP0[7][0] = -1;
  aP0[7][1] = +1;
  aP0[7][2] = +1;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist01(0, 1);
  for (int ino = 0; ino < 8; ++ino) {
    for (int idim = 0; idim < 3; ++idim) {
      aP0[ino][idim] += 0.1 * dist01(rndeng);
      aU0[ino][idim] = 0.2 * dist01(rndeng);
    }
  }
}

TEST(fem_solidhyper, Check_WdWddW_Mooneyrivlin2reduced) {
  double aP0[8][3], aU0[8][3];
  MakeRandomHex(aP0, aU0);
//
  double eps = 1.0e-6;
  double c1 = 1.0;
  double c2 = 1.0;
  double W0, dW0[8][3], ddW0[8][8][3][3];
  {
    W0 = 0.0;
    for (int i = 0; i < 8 * 8 * 3 * 3; ++i) { (&ddW0[0][0][0][0])[i] = 0.0; }
    for (int i = 0; i < 8 * 3; ++i) { (&dW0[0][0])[i] = 0.0; }
    double vol0 = 0.0;
    delfem2::AddWdWddW_Solid3HyperMooneyrivlin2Reduced_Hex(
        W0, dW0, ddW0, vol0,
        c1, c2, aP0, aU0, 1);
  }
  for (int ino = 0; ino < 8; ++ino) {
    for (int idim = 0; idim < 3; ++idim) {
      double W1, dW1[8][3];
      {
        double ddW1[8][8][3][3];
        W1 = 0.0;
        for (int i = 0; i < 8 * 8 * 3 * 3; ++i) { (&ddW1[0][0][0][0])[i] = 0.0; }
        for (int i = 0; i < 8 * 3; ++i) { (&dW1[0][0])[i] = 0.0; }
        double vol1 = 0.0;
        double aU1[8][3];
        for (int jno = 0; jno < 8; ++jno) {
          for (int jdim = 0; jdim < 3; ++jdim) {
            aU1[jno][jdim] = aU0[jno][jdim];
          }
        }
        aU1[ino][idim] += eps;
        delfem2::AddWdWddW_Solid3HyperMooneyrivlin2Reduced_Hex(
            W1, dW1, ddW1, vol1,
            c1, c2, aP0, aU1, 1);
      }
      {
        EXPECT_NEAR((W1 - W0) / eps, dW0[ino][idim], 1.0e-4);
      }
      for (int jno = 0; jno < 8; ++jno) {
        for (int jdim = 0; jdim < 3; ++jdim) {
          double tmp0 = (dW1[jno][jdim] - dW0[jno][jdim]) / eps;
          EXPECT_NEAR(tmp0, ddW0[ino][jno][idim][jdim], 1.0e-4);
        }
      }
    }
  }
}

TEST(fem_solidhyper, Check_WdWddW_Compression) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist01(0, 1);
  //
  double aP0[8][3], aU0[8][3];
  MakeRandomHex(aP0, aU0);
  for (int ino = 0; ino < 8; ++ino) {
    for (int idim = 0; idim < 3; ++idim) {
      aP0[ino][idim] += 0.1 * dist01(rndeng);
      aU0[ino][idim] += 0.2 * dist01(rndeng);
    }
  }
  double eps = 1.0e-6;
  double c1 = 1.0;
  double c2 = 1.0;
  double W0, dW0[8][3], ddW0[8][8][3][3];
  {
    W0 = 0.0;
    for (int i = 0; i < 8 * 8 * 3 * 3; ++i) { (&ddW0[0][0][0][0])[i] = 0.0; }
    for (int i = 0; i < 8 * 3; ++i) { (&dW0[0][0])[i] = 0.0; }
    double vol0 = 0.0;
    delfem2::AddWdWddW_Solid3Compression_Hex(
        W0, dW0, ddW0, vol0,
        1., aP0, aU0, 1);
  }
  for (int ino = 0; ino < 8; ++ino) {
    for (int idim = 0; idim < 3; ++idim) {
      double W1, dW1[8][3];
      {
        double ddW1[8][8][3][3];
        W1 = 0.0;
        for (int i = 0; i < 8 * 8 * 3 * 3; ++i) { (&ddW1[0][0][0][0])[i] = 0.0; }
        for (int i = 0; i < 8 * 3; ++i) { (&dW1[0][0])[i] = 0.0; }
        double vol1 = 0.0;
        double aU1[8][3];
        for (int jno = 0; jno < 8; ++jno) {
          for (int jdim = 0; jdim < 3; ++jdim) {
            aU1[jno][jdim] = aU0[jno][jdim];
          }
        }
        aU1[ino][idim] += eps;
        delfem2::AddWdWddW_Solid3Compression_Hex(
            W1, dW1, ddW1, vol1,
            1., aP0, aU1, 1);
      }
      {
        EXPECT_NEAR((W1 - W0) / eps, dW0[ino][idim], 1.0e-4);
      }
      for (int jno = 0; jno < 8; ++jno) {
        for (int jdim = 0; jdim < 3; ++jdim) {
          double tmp0 = (dW1[jno][jdim] - dW0[jno][jdim]) / eps;
          EXPECT_NEAR(tmp0, ddW0[ino][jno][idim][jdim], 1.0e-4);
        }
      }
    }
  }
}
