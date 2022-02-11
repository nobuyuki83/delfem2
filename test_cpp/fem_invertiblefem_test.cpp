//
// Created by Nobuyuki Umetani on 2022/02/01.
//

#include "gtest/gtest.h"

#include "delfem2/fem_invertiblefem.h"
#include "delfem2/sampling.h"
#include "delfem2/mat3_funcs.h"

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

    const dfm2::CMat3d F0 = DeformationGradientOfTet(
        Pos0[0], Pos0[1], Pos0[2], Pos0[3],
        pos0[0], pos0[1], pos0[2], pos0[3]);

    double dFdu[4][3];
    DiffDeformationGradientOfTet(
        dFdu,
        Pos0[0], Pos0[1], Pos0[2], Pos0[3]);

    const dfm2::CMat3d Basis0 = dfm2::Mat3_From3Bases(
        Pos0[1] - Pos0[0],
        Pos0[2] - Pos0[0],
        Pos0[3] - Pos0[0]);

    // check dFdu
    for (unsigned int ino = 0; ino < 4; ++ino) {
      for (unsigned int idim = 0; idim < 3; ++idim) {
        const dfm2::CMat3d A = dfm2::RandomVec<9>(dist_01, rndeng);
        dfm2::CVec3d pos1[4] = {pos0[0], pos0[1], pos0[2], pos0[3]};
        pos1[ino][idim] += eps;
        const dfm2::CMat3d basis1 = dfm2::Mat3_From3Bases(
            pos1[1] - pos1[0],
            pos1[2] - pos1[0],
            pos1[3] - pos1[0]);
        const dfm2::CMat3d F1 = basis1 * Basis0.Inverse();
        const double v0 = (A.transpose() * (F1 - F0)).trace() / eps;
        const double v1 = A(idim, 0) * dFdu[ino][0] + A(idim, 1) * dFdu[ino][1] + A(idim, 2) * dFdu[ino][2];
        EXPECT_NEAR(v0, v1, 1.0e-5);
      }
    }
  }
}

TEST(invertible_fem, dPdF) {
  namespace dfm2 = delfem2;

  std::mt19937 rndeng(0);
  std::uniform_real_distribution<double> dist_01(0, 1);
  double eps = 1.0e-5;

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

  for (unsigned int itr = 0; itr < 1000; ++itr) {
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

    const dfm2::CMat3d F0 = DeformationGradientOfTet(
        Pos0[0], Pos0[1], Pos0[2], Pos0[3],
        pos0[0], pos0[1], pos0[2], pos0[3]);

    double dFdu[4][3];
    DiffDeformationGradientOfTet(
        dFdu,
        Pos0[0], Pos0[1], Pos0[2], Pos0[3]);

    const auto[U0, S0, V0] = dfm2::Svd3<dfm2::CMat3d>(F0, 30);
    EXPECT_LT((U0 * (S0 * V0.transpose()) - F0).squaredNorm(), 1.0e-20);
    double diff[9][3][3];
    Svd3Differential(
        diff,
        U0, S0, V0);

    dfm2::CMat3d P0;
    double dPdF[3][3][3][3];
    double W0 = DiffPiolaKirchhoff1st(
        P0, dPdF,
        U0, S0, V0, diff, neohook);

    double dWdl0[3];
    {
      double ddWddl[3][3];
      const double Wt = neohook(dWdl0, ddWddl, S0(0, 0), S0(1, 1), S0(2, 2));
    }

    // check dPdF
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
            EXPECT_NEAR(v0, dPdF[k][l][i][j], 2.0e-5);
            EXPECT_NEAR(dPdF[k][l][i][j], dPdF[i][j][k][l], 1.0e-5);  // dPdF is symmetric
          }
        }
      }
    }
  }
}

TEST(invertible_fem, energy) {
  namespace dfm2 = delfem2;

  std::mt19937 rndeng(0);
  std::uniform_real_distribution<double> dist_01(0, 1);
  double eps = 1.0e-5;

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

  for (unsigned int itr = 0; itr < 10000; ++itr) {
    // make random rest shape of a tetrahedron
    double Pos0[4][3];
    delfem2::Fill2dArrayWithRandomValue<4,3>(Pos0, dist_01, rndeng);
    double Vol = dfm2::Volume_Tet<double [3], double>(Pos0[0], Pos0[1], Pos0[2], Pos0[3]);
    if (Vol < 0.01) { continue; }

    // make random deformed shape of a tetrahedron
    double pos0[4][3];
    delfem2::Fill2dArrayWithRandomValue<4,3>(pos0, dist_01, rndeng);

    double dW0[4][3], ddW0[4][4][3][3];
    const double W0 = delfem2::WdWddW_InvertibleFEM(
        dW0, ddW0,
        Pos0, pos0, neohook);

    for (unsigned int ino = 0; ino < 4; ++ino) {
      for (unsigned int jno = 0; jno < 4; ++jno) {
        for (unsigned int idim = 0; idim < 3; ++idim) {
          for (unsigned int jdim = 0; jdim < 3; ++jdim) {
            const double v0 = ddW0[ino][jno][idim][jdim];
            const double v1 = ddW0[jno][ino][jdim][idim];
            EXPECT_NEAR(v0, v1, 1.0e-5);
          }
        }
      }
    }

    for (unsigned int ino = 0; ino < 4; ++ino) {
      for (unsigned int idim = 0; idim < 3; ++idim) {
        double pos1[4][3] = {
            {pos0[0][0], pos0[0][1], pos0[0][2]},
            {pos0[1][0], pos0[1][1], pos0[1][2]},
            {pos0[2][0], pos0[2][1], pos0[2][2]},
            {pos0[3][0], pos0[3][1], pos0[3][2]} };
        pos1[ino][idim] += eps;
        //
        double dW1[4][3], ddW1[4][4][3][3];
        const double W1 = delfem2::WdWddW_InvertibleFEM(
            dW1, ddW1,
            Pos0, pos1, neohook);
        {
          double v0 = (W1 - W0) / eps;
          double v1 = dW0[ino][idim];
          EXPECT_NEAR(v0, v1, 2.0e-3 * (1.0 + fabs(v1)));
        }
        for (unsigned int jno = 0; jno < 4; ++jno) {
          for (unsigned int jdim = 0; jdim < 3; ++jdim) {
            double v0 = (dW1[jno][jdim] - dW0[jno][jdim]) / eps;
            double v1 = ddW0[jno][ino][jdim][idim];
            EXPECT_NEAR(v0, v1, 5.0e-2 * (1.0 + fabs(v1)));
          }
        }
      }
    }

  }
}