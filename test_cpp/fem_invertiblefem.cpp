//
// Created by Nobuyuki Umetani on 2022/02/01.
//

#include "gtest/gtest.h"

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/geosolidelm_v3.h"
#include "delfem2/sampling.h"
#include "delfem2/svd3.h"


TEST(invertible_fem, dF) {
  namespace dfm2 = delfem2;

  std::mt19937 rndeng(0);
  std::uniform_real_distribution<double> dist_01(0, 1);
  double eps = 1.0e-5;

  for (unsigned int itr = 0; itr < 1; ++itr) {
    // make random rest shape of a tetrahedron
    const dfm2::CVec3d Pos0[4] = {
        delfem2::RandomVec<3>(dist_01, rndeng),
        delfem2::RandomVec<3>(dist_01, rndeng),
        delfem2::RandomVec<3>(dist_01, rndeng),
        delfem2::RandomVec<3>(dist_01, rndeng)};
    double Vol = dfm2::Volume_Tet(Pos0[0], Pos0[1], Pos0[2], Pos0[3]);
    if (Vol < 0.01) { continue; }

    // make random deformed shape of a tetrahedron
    const dfm2::CVec3d pos0[4] = {
        Pos0[0] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2,
        Pos0[1] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2,
        Pos0[2] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2,
        Pos0[3] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2};

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
  }
}

TEST(invertible_fem, test0) {
  namespace dfm2 = delfem2;

  std::mt19937 rndeng(0);
  std::uniform_real_distribution<double> dist_01(0, 1);
  double eps = 1.0e-5;

  for (unsigned int itr = 0; itr < 1; ++itr) {
    // make random rest shape of a tetrahedron
    const dfm2::CVec3d Pos0[4] = {
        delfem2::RandomVec<3>(dist_01, rndeng),
        delfem2::RandomVec<3>(dist_01, rndeng),
        delfem2::RandomVec<3>(dist_01, rndeng),
        delfem2::RandomVec<3>(dist_01, rndeng)};
    double Vol = dfm2::Volume_Tet(Pos0[0], Pos0[1], Pos0[2], Pos0[3]);
    if (Vol < 0.01) { continue; }

    // make random deformed shape of a tetrahedron
    const dfm2::CVec3d pos0[4] = {
        Pos0[0] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2,
        Pos0[1] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2,
        Pos0[2] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2,
        Pos0[3] + dfm2::CVec3d(dfm2::RandomVec<3>(dist_01, rndeng)) * 0.2};

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

    const auto[U0, S0, V0] = dfm2::Svd3<dfm2::CMat3d>(F0, 30);
    EXPECT_LT((U0 * (S0 * V0.transpose()) - F0).squaredNorm(), 1.0e-20);
    double diff[9][3][3];
    Svd3Differential(
        diff,
        U0, S0, V0);

    auto neohook = [](
        double dW[3], double ddW[3][3],
        double l0, double l1, double l2) -> double {
      dW[0] = (l0 - 1);
      dW[1] = (l1 - 1);
      dW[2] = (l2 - 1);
      ddW[0][0] = 1;
      ddW[0][1] = 0;
      ddW[0][2] = 0;
      ddW[1][0] = 0;
      ddW[1][1] = 1;
      ddW[1][2] = 0;
      ddW[2][0] = 0;
      ddW[2][1] = 0;
      ddW[2][2] = 1;
      return (l0 - 1) * (l0 - 1) * 0.5
          + (l1 - 1) * (l1 - 1) * 0.5
          + (l2 - 1) * (l2 - 1) * 0.5;
    };

    // energy
    double dWdl0[3], ddWddl0[3][3];
    const double W0 = neohook(dWdl0, ddWddl0, S0(0, 0), S0(1, 1), S0(2, 2));
    const dfm2::CMat3d T0(dWdl0[0], dWdl0[1], dWdl0[2]);
    const dfm2::CMat3d P0 = U0 * (T0 * V0.transpose());
    double dPdF[3][3][3][3];  // dPdF[k][l][i][j] = dP_kl / dF_ij
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        const dfm2::CMat3d dU = -U0 * dfm2::CMat3d::Spin(
            std::array<double, 3>{diff[0][i][j], diff[1][i][j], diff[2][i][j]}.data());
        const dfm2::CMat3d dV = -V0 * dfm2::CMat3d::Spin(
            std::array<double, 3>{diff[6][i][j], diff[7][i][j], diff[8][i][j]}.data());
        const dfm2::CMat3d dS(diff[3][i][j], diff[4][i][j], diff[5][i][j]);
        double t0 = ddWddl0[0][0] * dS(0,0) + ddWddl0[0][1] * dS(1,1) + ddWddl0[0][2] * dS(2,2);
        double t1 = ddWddl0[1][0] * dS(0,0) + ddWddl0[1][1] * dS(1,1) + ddWddl0[1][2] * dS(2,2);
        double t2 = ddWddl0[2][0] * dS(0,0) + ddWddl0[2][1] * dS(1,1) + ddWddl0[2][2] * dS(2,2);
        const dfm2::CMat3d dT(t0,t1,t2);
        const dfm2::CMat3d dP = dU * T0 * V0.transpose()
            + U0 * dT * V0.transpose()
            + U0 * (T0 * dV.transpose());
        for (unsigned k = 0; k < 3; ++k) {
          for (unsigned l = 0; l < 3; ++l) {
            dPdF[k][l][i][j] = dP(k, l);
          }
        }
      }
    }

    // check SVD3 differential
    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        dfm2::CMat3d F1 = F0;
        F1(i, j) += eps;
        const auto[U1, S1, V1] = dfm2::Svd3<dfm2::CMat3d>(F1, 30);
        double dWdl1[3], ddWddl1[3][3];
        const double W1 = neohook(dWdl1, ddWddl1, S1(0, 0), S1(1, 1), S1(2, 2));
        EXPECT_NEAR((W1 - W0) / eps, P0(i, j), 1.0e-3);
        {
          const double dWdF = dWdl0[0] * diff[3][i][j]
              + dWdl0[1] * diff[4][i][j]
              + dWdl0[2] * diff[5][i][j];
          EXPECT_NEAR(dWdF, P0(i, j), 1.0e-8);
        }
        dfm2::CMat3d T1(dWdl1[0], dWdl1[1], dWdl1[2]);
        dfm2::CMat3d P1 = U1 * (T1 * V1.transpose());
        for (unsigned int k = 0; k < 3; ++k) {
          for (unsigned int l = 0; l < 3; ++l) {
            double v0 = (P1 - P0)(k, l) / eps;
            EXPECT_NEAR(v0, dPdF[k][l][i][j], 1.0e-5);
            EXPECT_NEAR(dPdF[k][l][i][j], dPdF[i][j][k][l], 1.0e-5);  // dPdF is symmetric
          }
        }
      }
    }
  }
}

