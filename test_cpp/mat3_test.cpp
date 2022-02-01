

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/mat3.h"
#include "delfem2/geosolidelm_v3.h"
#include "delfem2/sampling.h"
#include "delfem2/svd3.h"

namespace dfm2 = delfem2;

TEST(mat3, eigen3) {
  std::random_device randomDevice;
  std::mt19937 rdeng(randomDevice());
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
      dfm2::CVec3d qv0 = dfm2::QuatVec(quat0, v0);
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
      EXPECT_NEAR(ea0[0], ea1[0], 2.0e-10);
      EXPECT_NEAR(ea0[1], ea1[1], 2.0e-10);
      EXPECT_NEAR(ea0[2], ea1[2], 2.0e-10);
    }
    {
      double quat0[4];
      delfem2::Quaternion_EulerAngle(quat0, {ea0[0], ea0[1], ea0[2]}, {2, 0, 1});
      double mat0[9];
      delfem2::Mat3_Quat(mat0, quat0);
      double ea1[3];
      delfem2::EulerAngle_Mat3(ea1, mat0, {2, 0, 1});
      EXPECT_NEAR(ea0[0], ea1[0], 2.0e-10);
      EXPECT_NEAR(ea0[1], ea1[1], 2.0e-10);
      EXPECT_NEAR(ea0[2], ea1[2], 2.0e-10);
    }
  }
}

TEST(mat3, diff_svd) {

  std::mt19937 rndeng(0);
  std::uniform_real_distribution<double> dist_01(0, 1);
  double eps = 1.0e-5;

  for (unsigned int itr = 0; itr < 10000; ++itr) {
    // make random rest shape of a tetrahedron
    const dfm2::CVec3d Pos0[4] = {
        delfem2::RandomVec3(dist_01, rndeng),
        delfem2::RandomVec3(dist_01, rndeng),
        delfem2::RandomVec3(dist_01, rndeng),
        delfem2::RandomVec3(dist_01, rndeng)};
    double Vol = dfm2::Volume_Tet(Pos0[0], Pos0[1], Pos0[2], Pos0[3]);
    if (Vol < 0.01) { continue; }

    // make random deformed shape of a tetrahedron
    const dfm2::CVec3d pos0[4] = {
        Pos0[0] + dfm2::CVec3d(dfm2::RandomVec3(dist_01, rndeng)) * 0.2,
        Pos0[1] + dfm2::CVec3d(dfm2::RandomVec3(dist_01, rndeng)) * 0.2,
        Pos0[2] + dfm2::CVec3d(dfm2::RandomVec3(dist_01, rndeng)) * 0.2,
        Pos0[3] + dfm2::CVec3d(dfm2::RandomVec3(dist_01, rndeng)) * 0.2};

    const dfm2::CMat3d Basis0 = dfm2::Mat3_From3Bases(
        Pos0[1] - Pos0[0],
        Pos0[2] - Pos0[0],
        Pos0[3] - Pos0[0]);
    const dfm2::CMat3d basis0 = dfm2::Mat3_From3Bases(
        pos0[1] - pos0[0],
        pos0[2] - pos0[0],
        pos0[3] - pos0[0]);
    const dfm2::CMat3d F0 = basis0 * Basis0.Inverse(); // deformation gradient tensor

    const std::array<dfm2::CVec3d, 4> dFdu = DiffDeformationGradient(
        Pos0[0], Pos0[1], Pos0[2], Pos0[3]);

    // check dFdu
    for (unsigned int ino = 0; ino < 4; ++ino) {
      for (unsigned int idim = 0; idim < 3; ++idim) {
        const dfm2::CMat3d A = dfm2::RandomMat3(dist_01, rndeng);
        dfm2::CVec3d pos1[4] = {pos0[0], pos0[1], pos0[2], pos0[3]};
        pos1[ino][idim] += eps;
        const dfm2::CMat3d basis1 = dfm2::Mat3_From3Bases(
            pos1[1] - pos1[0],
            pos1[2] - pos1[0],
            pos1[3] - pos1[0]);
        const dfm2::CMat3d F1 = basis1 * Basis0.Inverse();
        const double v0 = (A.transpose() * (F1 - F0)).trace() / eps;
        const double v1 = (A * dFdu[ino])[idim];
        EXPECT_NEAR(v0, v1, 1.0e-5);
      }
    }

    const auto[U0, S0, V0] = dfm2::Svd3<dfm2::CMat3d>(F0, 30);
    EXPECT_LT((U0 * (S0 * V0.transpose()) - F0).squaredNorm(), 1.0e-20);
    double diff[3][3][9];
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
          EXPECT_NEAR(dU(1, 2), diff[i][j][0], 2.0e-3 * (1.0 + fabs(diff[i][j][0])));
          EXPECT_NEAR(dU(2, 0), diff[i][j][1], 2.0e-3 * (1.0 + fabs(diff[i][j][1])));
          EXPECT_NEAR(dU(0, 1), diff[i][j][2], 2.0e-3 * (1.0 + fabs(diff[i][j][2])));
        }
        {  // check diff. of singular value
          double v0 = (S1(0, 0) - S0(0, 0)) / eps;
          double v1 = (S1(1, 1) - S0(1, 1)) / eps;
          double v2 = (S1(2, 2) - S0(2, 2)) / eps;
          EXPECT_NEAR(v0, diff[i][j][3], 5.0e-4 * (1.0 + fabs(diff[i][j][3])));
          EXPECT_NEAR(v1, diff[i][j][4], 5.0e-4 * (1.0 + fabs(diff[i][j][4])));
          EXPECT_NEAR(v2, diff[i][j][5], 5.0e-4 * (1.0 + fabs(diff[i][j][5])));
        }
        {   // check diff. of V. Asym(diff) = V^T V'
          dfm2::CMat3d dV = V0.transpose() * (V1 - V0) / eps;
          EXPECT_NEAR(dV(1, 2), diff[i][j][6], 2.0e-3 * (1.0 + fabs(diff[i][j][6])));
          EXPECT_NEAR(dV(2, 0), diff[i][j][7], 2.0e-3 * (1.0 + fabs(diff[i][j][7])));
          EXPECT_NEAR(dV(0, 1), diff[i][j][8], 2.0e-3 * (1.0 + fabs(diff[i][j][8])));
        }
      }
    }
  }
}