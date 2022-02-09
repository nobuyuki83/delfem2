

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/mat3.h"
#include "delfem2/mat3_funcs.h"
#include "delfem2/sampling.h"
#include "delfem2/svd3.h"

namespace dfm2 = delfem2;

TEST(mat3, eigen3) {
  std::mt19937 rdeng(std::random_device{}());
  std::uniform_real_distribution<double> dist(-50.0, 50.0);
  for (int itr = 0; itr < 10000; itr++) {
    double sm[6];
    for (double &v: sm) {
      v = dist(rdeng);
    }
    double l[3];
    dfm2::CMat3d U;
    dfm2::EigenSym3(
        U.data(), l,
        sm, 20);
    {
      double diffU = (U.transpose() * U - dfm2::CMat3d::Identity()).squaredNorm();
      EXPECT_NEAR(diffU, 0.0, 1.0e-10);
    }
    {
      double L[9] = {l[0], 0, 0, 0, l[1], 0, 0, 0, l[2]};
      dfm2::CMat3d UL;
      dfm2::MatMat3(UL.data(), U.data(), L);
      dfm2::CMat3d ULUt;
      dfm2::MatMatT3(ULUt.data(), UL.data(), U.data());
      dfm2::CMat3d SM;
      SM.SetSymetric(sm);
      double diff = (ULUt - SM).squaredNorm();
      EXPECT_NEAR(diff, 0.0, 1.0e-6);
    }
  }
  // -----------------------------
  for (int itr = 0; itr < 100; itr++) {
    double sm[6];
    for (double &v: sm) {
      v = dist(rdeng);
    }
    sm[5] = -sm[4];
    double l[3];
    dfm2::CMat3d U;
    dfm2::EigenSym3(
        U.data(), l,
        sm, 20);
    {
      double diffU = (U.transpose() * U - dfm2::CMat3d::Identity()).squaredNorm();
      EXPECT_NEAR(diffU, 0.0, 1.0e-10);
    }
    {
      double L[9] = {l[0], 0, 0, 0, l[1], 0, 0, 0, l[2]};
      dfm2::CMat3d UL;
      dfm2::MatMat3(UL.data(), U.data(), L);
      dfm2::CMat3d ULUt;
      dfm2::MatMatT3(ULUt.data(), UL.data(), U.data());
      dfm2::CMat3d SM;
      SM.SetSymetric(sm);
      double diff = (ULUt - SM).squaredNorm();
      EXPECT_NEAR(diff, 0.0, 1.0e-6);
    }
  }
}

TEST(mat3, svd3) {
  for (int itr = 0; itr < 10000; itr++) {
    dfm2::CMat3d M;
    M.SetRandom();
    auto[U, UG, V] = dfm2::Svd3<dfm2::CMat3d>(M, 20);
    {
      double diffU = (U.transpose() * U - dfm2::CMat3d::Identity()).squaredNorm();
      EXPECT_NEAR(diffU, 0.0, 1.0e-6);
    }
    {
      double diffV = (V.transpose() * V - dfm2::CMat3d::Identity()).squaredNorm();
      EXPECT_NEAR(diffV, 0.0, 1.0e-10);
    }
    {
      dfm2::CMat3d UGVt = U * UG * V.transpose();
      double diff = (UGVt - M).squaredNorm();
      EXPECT_NEAR(diff, 0.0, 1.0e-10);
    }
  }
}

TEST(mat3, rot_comp) {
  for (int itr = 0; itr < 10000; itr++) {
    dfm2::CMat3d M;
    M.SetRandom();
    dfm2::CMat3d R;
    dfm2::GetRotPolarDecomp(R.data(), M.data(), 40);
    {
      double diff = (R.transpose() * R - dfm2::CMat3d::Identity()).squaredNorm();
      EXPECT_NEAR(diff, 0.0, 1.0e-5);
    }
    {
      dfm2::CMat3d MR = M.MatMat(R.transpose());
      double diff0 = (MR - MR.Sym()).squaredNorm();
      EXPECT_NEAR(diff0, 0.0, 1.0e-5);
    }
    {
      dfm2::CMat3d RM = (R.transpose()).MatMat(M);
      double diff1 = (RM - RM.Sym()).squaredNorm();
      EXPECT_NEAR(diff1, 0.0, 1.0e-5);
    }
  }
}

TEST(mat3, quat) {
  std::random_device randomDevice;
  std::uniform_real_distribution<double> dist(-50.0, +50.0);
  std::mt19937 mtd(randomDevice());
  for (int itr = 0; itr < 10000; itr++) {
    double quat0[4] = {dist(mtd), dist(mtd), dist(mtd), dist(mtd)};
    dfm2::Normalize_Quat(quat0);
    dfm2::CMat3d R0;
    R0.SetRotMatrix_Quaternion(quat0);
    {
      double diff = (R0.transpose() * R0 - dfm2::CMat3d::Identity()).squaredNorm();
      EXPECT_NEAR(diff, 0.0, 1.0e-14);
    }
    { // q0 -> R0 -> q1 -> R1
      const std::array<double, 4> quat1 = R0.GetQuaternion();
      dfm2::CMat3d R1;
      R1.SetRotMatrix_Quaternion(quat1.data());
      double diff = (R1 - R0).squaredNorm();
      EXPECT_NEAR(diff, 0.0, 1.0e-20);
    }
    {
      dfm2::CVec3d v0(dist(mtd), dist(mtd), dist(mtd));
      dfm2::CVec3d qv0 = dfm2::QuatVec3(quat0, v0);
      dfm2::CVec3d Rv0 = dfm2::MatVec(R0, v0);
      EXPECT_LT((qv0 - Rv0).norm(), 1.0e-20);
    }
  }

}

TEST(mat3, mat3_quat_eulerangle) {
  std::uniform_real_distribution<double> dist_m1p1(-1, 1);
  std::mt19937 mtd(std::random_device{}());
  for (int itr = 0; itr < 10000; itr++) {
    double ea0[3] = {
        dist_m1p1(mtd) * M_PI * 0.5,
        dist_m1p1(mtd) * M_PI * 0.5,
        dist_m1p1(mtd) * M_PI * 0.5};
    {
      double quat0[4];
      delfem2::Quaternion_EulerAngle(quat0, {ea0[0], ea0[1], ea0[2]}, {2, 1, 0});
      double mat0[9];
      delfem2::Mat3_Quat(mat0, quat0);
      double ea1[3];
      delfem2::EulerAngle_Mat3(ea1, mat0, {2, 1, 0});
      EXPECT_NEAR(ea0[0], ea1[0], 1.0e-9);
      EXPECT_NEAR(ea0[1], ea1[1], 1.0e-9);
      EXPECT_NEAR(ea0[2], ea1[2], 1.0e-9);
    }
    {
      double quat0[4];
      delfem2::Quaternion_EulerAngle(quat0, {ea0[0], ea0[1], ea0[2]}, {2, 0, 1});
      double mat0[9];
      delfem2::Mat3_Quat(mat0, quat0);
      double ea1[3];
      delfem2::EulerAngle_Mat3(ea1, mat0, {2, 0, 1});
      EXPECT_NEAR(ea0[0], ea1[0], 1.0e-9);
      EXPECT_NEAR(ea0[1], ea1[1], 1.0e-9);
      EXPECT_NEAR(ea0[2], ea1[2], 1.0e-9);
    }
  }
}

TEST(mat3, diff_svd) {

  std::mt19937 rndeng(0);
  std::uniform_real_distribution<double> dist_01(0, 1);
  double eps = 1.0e-5;

  for (unsigned int itr = 0; itr < 1; ++itr) {
    // make random rest shape of a tetrahedron
    const dfm2::CMat3d F0 = dfm2::RandomVec<9>(dist_01, rndeng);

    const auto[U0, S0, V0] = dfm2::Svd3<dfm2::CMat3d>(F0, 30);
    EXPECT_LT((U0 * (S0 * V0.transpose()) - F0).squaredNorm(), 1.0e-20);
    double diff[9][3][3];
    Svd3Differential(
        diff,
        U0, S0, V0);

    // check SVD3 differential
    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        dfm2::CMat3d F1 = F0;
        F1(i, j) += eps;
        const auto[U1, S1, V1] = dfm2::Svd3<dfm2::CMat3d>(F1, 30);
        {  // check diff. of U.  Asym(diff) = U^T U'
          dfm2::CMat3d dU = U0.transpose() * (U1 - U0) / eps;
          EXPECT_NEAR(dU(1, 2), diff[0][i][j], 2.0e-3 * (1.0 + fabs(diff[0][i][j])));
          EXPECT_NEAR(dU(2, 0), diff[1][i][j], 2.0e-3 * (1.0 + fabs(diff[1][i][j])));
          EXPECT_NEAR(dU(0, 1), diff[2][i][j], 2.0e-3 * (1.0 + fabs(diff[2][i][j])));
        }
        {  // check diff. of singular value
          double v0 = (S1(0, 0) - S0(0, 0)) / eps;
          double v1 = (S1(1, 1) - S0(1, 1)) / eps;
          double v2 = (S1(2, 2) - S0(2, 2)) / eps;
          EXPECT_NEAR(v0, diff[3][i][j], 5.0e-4 * (1.0 + fabs(diff[3][i][j])));
          EXPECT_NEAR(v1, diff[4][i][j], 5.0e-4 * (1.0 + fabs(diff[4][i][j])));
          EXPECT_NEAR(v2, diff[5][i][j], 5.0e-4 * (1.0 + fabs(diff[5][i][j])));
        }
        {   // check diff. of V. Asym(diff) = V^T V'
          dfm2::CMat3d dV = V0.transpose() * (V1 - V0) / eps;
          EXPECT_NEAR(dV(1, 2), diff[6][i][j], 2.0e-3 * (1.0 + fabs(diff[6][i][j])));
          EXPECT_NEAR(dV(2, 0), diff[7][i][j], 2.0e-3 * (1.0 + fabs(diff[7][i][j])));
          EXPECT_NEAR(dV(0, 1), diff[8][i][j], 2.0e-3 * (1.0 + fabs(diff[8][i][j])));
        }
      }
    }
  }
}
