#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/fem_rod3_darboux.h"
#include "delfem2/geo3_v23m34q.h"

namespace dfm2 = delfem2;

bool GenRandomConfigRod(
    dfm2::CVec3d P0[3],
    dfm2::CVec3d S0[2],
    std::mt19937 &randomEng,
    std::uniform_real_distribution<double> &dist_m1p1) {
  P0[0] = dfm2::CVec3d::Random(dist_m1p1, randomEng);
  P0[1] = dfm2::CVec3d::Random(dist_m1p1, randomEng);
  P0[2] = dfm2::CVec3d::Random(dist_m1p1, randomEng);
  if ((P0[1] - P0[0]).norm() < 0.1) { return false; }
  if ((P0[2] - P0[1]).norm() < 0.1) { return false; }
  if ((P0[1] - P0[0]).normalized().dot((P0[2] - P0[1]).normalized()) < -0.1) { return false; }
  {
    S0[0].SetRandom(dist_m1p1, randomEng);
    const dfm2::CVec3d U0 = (P0[1] - P0[0]).normalized();
    S0[0] -= (S0[0].dot(U0)) * U0;
    S0[0].normalize();
  }
  {
    S0[1].SetRandom(dist_m1p1, randomEng);
    const dfm2::CVec3d U1 = (P0[2] - P0[1]).normalized();
    S0[1] -= (S0[1].dot(U1)) * U1;
    S0[1].normalize();
  }
  if (S0[0].dot(S0[1]) < -0.1) { return false; }
  return true;
}

TEST(fem_rod, dWddW_RodFrameTrans) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0.0, 1.0);
  std::uniform_real_distribution<double> dist_m1p1(-1.0, +1.0);
  const double eps = 1.0e-6;
//
  for (int itr = 0; itr < 10000; ++itr) {
    const dfm2::CVec3d V01 = dfm2::CVec3d::Random(dist_01, randomEng);
    if (V01.norm() < 0.1) { continue; }
    dfm2::CVec3d Frm[3];
    {
      Frm[2] = V01;
      Frm[2].normalize();
      Frm[0].SetRandom(dist_01, randomEng);
      Frm[0] -= (Frm[0].dot(Frm[2])) * Frm[2];
      Frm[0].normalize();
      Frm[1] = (Frm[2] ^ Frm[0]);
    }
    dfm2::CVec3d Q;
    Q.SetRandom(dist_01, randomEng);
//  Q = Frm[2];
// --------------------------------
    double W[3] = {Q.dot(Frm[0]), Q.dot(Frm[1]), Q.dot(Frm[2])};
    dfm2::CVec3d DW_Dv[3];
    double DW_Dt[3];
    {
      dfm2::CMat3d dF_dv[3];
      dfm2::CVec3d dF_dt[3];
      dfm2::DiffFrameRod(dF_dv, dF_dt,
                         V01.norm(), Frm);
      DW_Dv[0] = Q * dF_dv[0];
      DW_Dv[1] = Q * dF_dv[1];
      DW_Dv[2] = Q * dF_dv[2];
      DW_Dt[0] = Q.dot(dF_dt[0]);
      DW_Dt[1] = Q.dot(dF_dt[1]);
      DW_Dt[2] = Q.dot(dF_dt[2]);
    }
// ---------
    const dfm2::CVec3d du = dfm2::CVec3d::Random(dist_01, randomEng) * eps;
    const double dtheta = dist_m1p1(randomEng) * eps;
// ------
    dfm2::CVec3d frm[3];
    RodFrameTrans(frm,
                  Frm[0], V01,
                  du, dtheta);
    const double w[3] = {Q.dot(frm[0]), Q.dot(frm[1]), Q.dot(frm[2])};
    dfm2::CVec3d dw_dv[3];
    double dw_dt[3];
    {
      dfm2::CMat3d df_dv[3];
      dfm2::CVec3d df_dt[3];
      DiffFrameRod(df_dv, df_dt,
                   (V01 + du).norm(), frm);
      dw_dv[0] = Q * df_dv[0];
      dw_dv[1] = Q * df_dv[1];
      dw_dv[2] = Q * df_dv[2];
      dw_dt[0] = Q.dot(df_dt[0]);
      dw_dt[1] = Q.dot(df_dt[1]);
      dw_dt[2] = Q.dot(df_dt[2]);
    }
    for (int i = 0; i < 3; ++i) {
      double val0 = (w[i] - W[i]) / eps;
      double val1 = (DW_Dt[i] * dtheta + DW_Dv[i].dot(du)) / eps;
      EXPECT_NEAR(val0, val1, 1.0e-3);
    }
//
    for (int i = 0; i < 3; ++i) {
      dfm2::CMat3d DDW_DDv;
      dfm2::CVec3d DDW_DvDt;
      double DDW_DDt;
      DifDifFrameRod(DDW_DDv, DDW_DvDt, DDW_DDt,
                     i, V01.norm(), Q, Frm);
      double val0 = (dw_dt[i] - DW_Dt[i]) / eps;
      double val1 = (DDW_DDt * dtheta + DDW_DvDt.dot(du)) / eps;
      EXPECT_NEAR(val0, val1, 1.0e-3);
    }
    for (int i = 0; i < 3; ++i) {
      dfm2::CMat3d DDW_DDv;
      dfm2::CVec3d DDW_DvDt;
      double DDW_DDt;
      DifDifFrameRod(DDW_DDv, DDW_DvDt, DDW_DDt,
                     i, V01.norm(), Q, Frm);
      const dfm2::CVec3d vec0 = (dw_dv[i] - DW_Dv[i]) / eps;
      const dfm2::CVec3d vec1 = (DDW_DvDt * dtheta + DDW_DDv * du) / eps;
      EXPECT_LT((vec0 - vec1).norm(), 3.5e-3);
    }
  }
}

TEST(fem_rod, WdWddW_DotFrame) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0, 1);
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  const double eps = 1.0e-5;
//
  for (int itr = 0; itr < 10000; ++itr) {
    dfm2::CVec3d P[3] = {
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng)};
    if ((P[1] - P[0]).norm() < 0.01) { continue; }
    if ((P[2] - P[1]).norm() < 0.01) { continue; }
    if ((P[1] - P[0]).normalized().dot((P[2] - P[1]).normalized()) < -0.1) { continue; }
    dfm2::CVec3d S[2];
    {
      S[0].SetRandom(dist_01, randomEng);
      const dfm2::CVec3d U0 = (P[1] - P[0]).normalized();
      S[0] -= (S[0].dot(U0)) * U0;
      S[0].normalize();
    }
    {
      S[1].SetRandom(dist_01, randomEng);
      const dfm2::CVec3d U1 = (P[2] - P[1]).normalized();
      S[1] -= (S[1].dot(U1)) * U1;
      S[1].normalize();
    }
    const double off[3] = {
        dist_m1p1(randomEng),
        dist_m1p1(randomEng),
        dist_m1p1(randomEng)};
// ------------------------
    dfm2::CVec3d dW_dP[3];
    double dW_dt[2];
    dfm2::CMat3d ddW_ddP[3][3];
    dfm2::CVec3d ddW_dtdP[2][3];
    double ddW_ddt[2][2];
    double W = WdWddW_DotFrame(dW_dP, dW_dt,
                               ddW_ddP, ddW_dtdP, ddW_ddt,
                               P, S, off);
// -----------------------
    const dfm2::CVec3d dP[3] = {
        dfm2::CVec3d::Random(dist_01, randomEng).normalized() * eps,
        dfm2::CVec3d::Random(dist_01, randomEng).normalized() * eps,
        dfm2::CVec3d::Random(dist_01, randomEng).normalized() * eps};
    const double dT[2] = {
        dist_m1p1(randomEng) * eps,
        dist_m1p1(randomEng) * eps};
    dfm2::CVec3d frm0[3], frm1[3];
    RodFrameTrans(frm0,
                  S[0], P[1] - P[0], dP[1] - dP[0], dT[0]);
    RodFrameTrans(frm1,
                  S[1], P[2] - P[1], dP[2] - dP[1], dT[1]);
    const dfm2::CVec3d p[3] = {P[0] + dP[0], P[1] + dP[1], P[2] + dP[2]};
    const dfm2::CVec3d s[2] = {frm0[0], frm1[0]};
    dfm2::CVec3d dw_dP[3];
    double dw_dt[2];
    double w = 0;
    {
      dfm2::CMat3d ddw_ddP[3][3];
      dfm2::CVec3d ddw_dtdP[2][3];
      double ddw_ddt[2][2];
      w = WdWddW_DotFrame(dw_dP, dw_dt,
                          ddw_ddP, ddw_dtdP, ddw_ddt,
                          p, s, off);
    }
    {
      const double val0 = (w - W) / eps;
      const double val1 = (+dW_dt[0] * dT[0]
          + dW_dt[1] * dT[1]
          + dW_dP[0].dot(dP[0])
          + dW_dP[1].dot(dP[1])
          + dW_dP[2].dot(dP[2])) / eps;
      EXPECT_NEAR(val0, val1, 3.5e-2 * (1.0 + fabs(val1)));
    }
    {
      const double val0 = (dw_dt[0] - dW_dt[0]) / eps;
      const double val1 = (+ddW_ddt[0][0] * dT[0]
          + ddW_ddt[0][1] * dT[1]
          + ddW_dtdP[0][0].dot(dP[0])
          + ddW_dtdP[0][1].dot(dP[1])
          + ddW_dtdP[0][2].dot(dP[2])) / eps;
      EXPECT_NEAR(val0, val1, 3.0e-2 * (1.0 + fabs(val1)));
    }
    {
      const double val0 = (dw_dt[1] - dW_dt[1]) / eps;
      const double val1 = (+ddW_ddt[1][0] * dT[0]
          + ddW_ddt[1][1] * dT[1]
          + ddW_dtdP[1][0].dot(dP[0])
          + ddW_dtdP[1][1].dot(dP[1])
          + ddW_dtdP[1][2].dot(dP[2])) / eps;
      EXPECT_NEAR(val0, val1, 3.5e-2 * (1.0 + fabs(val1)));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[0] - dW_dP[0]) / eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][0] * dT[0]
          + ddW_dtdP[1][0] * dT[1]
          + ddW_ddP[0][0] * dP[0]
          + ddW_ddP[0][1] * dP[1]
          + ddW_ddP[0][2] * dP[2]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-2 * (1 + val1.norm()));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[1] - dW_dP[1]) / eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][1] * dT[0]
          + ddW_dtdP[1][1] * dT[1]
          + ddW_ddP[1][0] * dP[0]
          + ddW_ddP[1][1] * dP[1]
          + ddW_ddP[1][2] * dP[2]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-2 * (1 + val1.norm()));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[2] - dW_dP[2]) / eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][2] * dT[0]
          + ddW_dtdP[1][2] * dT[1]
          + ddW_ddP[2][0] * dP[0]
          + ddW_ddP[2][1] * dP[1]
          + ddW_ddP[2][2] * dP[2]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-2 * (1 + val1.norm()));
    }
  }
}

TEST(fem_rod, WdWddW_Rod) {
  std::random_device rd;
  std::mt19937 randomEng(rd());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  std::uniform_real_distribution<double> dist_12(+1, +2);
  std::uniform_real_distribution<double> dist_01(+0, +1);
  double eps = 1.0e-7;
//
  for (int itr = 0; itr < 1000; ++itr) {
    dfm2::CVec3d P[3], S[2];
    bool res = GenRandomConfigRod(P, S, randomEng, dist_m1p1);
    if (!res) { continue; }
    const double stiff_bendtwist[3] = {
        dist_12(randomEng),
        dist_12(randomEng),
        dist_12(randomEng)};
    const double off[3] = {
        dist_m1p1(randomEng),
        dist_m1p1(randomEng),
        dist_m1p1(randomEng)};
    // ------------------------
    dfm2::CVec3d dW_dP[3];
    double dW_dt[2];
    dfm2::CMat3d ddW_ddP[3][3];
    dfm2::CVec3d ddW_dtdP[2][3];
    double ddW_ddt[2][2];
    double W = WdWddW_Rod3Exact(dW_dP, dW_dt,
                          ddW_ddP, ddW_dtdP, ddW_ddt,
                          stiff_bendtwist, P, S,
                          dfm2::CVec3d(off));
// -----------------------
    const dfm2::CVec3d dP[3] = {
        dfm2::CVec3d::Random(dist_01, randomEng) * eps,
        dfm2::CVec3d::Random(dist_01, randomEng) * eps,
        dfm2::CVec3d::Random(dist_01, randomEng) * eps};
    const double dT[2] = {
        dist_m1p1(randomEng) * eps,
        dist_m1p1(randomEng) * eps};
    dfm2::CVec3d frm0[3], frm1[3];
    RodFrameTrans(frm0,
                  S[0], P[1] - P[0], dP[1] - dP[0], dT[0]);
    RodFrameTrans(frm1,
                  S[1], P[2] - P[1], dP[2] - dP[1], dT[1]);
    const dfm2::CVec3d p[3] = {P[0] + dP[0], P[1] + dP[1], P[2] + dP[2]};
    const dfm2::CVec3d s[2] = {frm0[0], frm1[0]};
    dfm2::CVec3d dw_dP[3];
    double dw_dt[2];
    double w = 0;
    {
      dfm2::CMat3d ddw_ddP[3][3];
      dfm2::CVec3d ddw_dtdP[2][3];
      double ddw_ddt[2][2];
      w = WdWddW_Rod3Exact(dw_dP, dw_dt,
                     ddw_ddP, ddw_dtdP, ddw_ddt,
                     stiff_bendtwist, p, s,
                     dfm2::CVec3d(off));
    }
    {
      const double val0 = (w - W) / eps;
      const double val1 = (+dW_dt[0] * dT[0]
          + dW_dt[1] * dT[1]
          + dW_dP[0].dot(dP[0])
          + dW_dP[1].dot(dP[1])
          + dW_dP[2].dot(dP[2])) / eps;
      EXPECT_NEAR(val0, val1, 5.0e-3 * (1 + fabs(val1)));
    }
    {
      const double val0 = (dw_dt[0] - dW_dt[0]) / eps;
      const double val1 = (+ddW_ddt[0][0] * dT[0]
          + ddW_ddt[0][1] * dT[1]
          + ddW_dtdP[0][0].dot(dP[0])
          + ddW_dtdP[0][1].dot(dP[1])
          + ddW_dtdP[0][2].dot(dP[2])) / eps;
      EXPECT_NEAR(val0, val1, 2.5e-3 * (1 + fabs(val1)));
    }
    {
      const double val0 = (dw_dt[1] - dW_dt[1]) / eps;
      const double val1 = (+ddW_ddt[1][0] * dT[0]
          + ddW_ddt[1][1] * dT[1]
          + ddW_dtdP[1][0].dot(dP[0])
          + ddW_dtdP[1][1].dot(dP[1])
          + ddW_dtdP[1][2].dot(dP[2])) / eps;
      EXPECT_NEAR(val0, val1, 1.0e-3 * (1 + fabs(val1)));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[0] - dW_dP[0]) / eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][0] * dT[0]
          + ddW_dtdP[1][0] * dT[1]
          + ddW_ddP[0][0] * dP[0]
          + ddW_ddP[0][1] * dP[1]
          + ddW_ddP[0][2] * dP[2]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-3 * (1 + val1.norm()));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[1] - dW_dP[1]) / eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][1] * dT[0]
          + ddW_dtdP[1][1] * dT[1]
          + ddW_ddP[1][0] * dP[0]
          + ddW_ddP[1][1] * dP[1]
          + ddW_ddP[1][2] * dP[2]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-3 * (1 + val1.norm()));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[2] - dW_dP[2]) / eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][2] * dT[0]
          + ddW_dtdP[1][2] * dT[1]
          + ddW_ddP[2][0] * dP[0]
          + ddW_ddP[2][1] * dP[1]
          + ddW_ddP[2][2] * dP[2]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-3 * (1 + val1.norm()));
    }
  }
}

TEST(fem_rod, CdC_Rod) {
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  const double eps = 1.0e-5;
  for (int itr = 0; itr < 10000; ++itr) {
    dfm2::CVec3d P0[3], S0[2];
    bool res = GenRandomConfigRod(P0, S0, randomEng, dist_m1p1);
    if (!res) { continue; }
    //
    double C0[3];
    dfm2::CVec3d dC0_dP[3][3];
    double dC0_dt[3][2];
    dfm2::CdC_Rod(
        C0, dC0_dP, dC0_dt,
        P0, S0);
    //
    const dfm2::CVec3d dP[3] = {
        dfm2::CVec3d::Random(dist_m1p1, randomEng) * eps,
        dfm2::CVec3d::Random(dist_m1p1, randomEng) * eps,
        dfm2::CVec3d::Random(dist_m1p1, randomEng) * eps};
    const double dT[2] = {
        dist_m1p1(randomEng) * eps,
        dist_m1p1(randomEng) * eps};
    dfm2::CVec3d P1[3], S1[2];
    {
      dfm2::CVec3d frm0[3], frm1[3];
      RodFrameTrans(
          frm0,
          S0[0], P0[1] - P0[0], dP[1] - dP[0], dT[0]);
      RodFrameTrans(
          frm1,
          S0[1], P0[2] - P0[1], dP[2] - dP[1], dT[1]);
      P1[0] = P0[0] + dP[0];
      P1[1] = P0[1] + dP[1];
      P1[2] = P0[2] + dP[2];
      S1[0] = frm0[0];
      S1[1] = frm1[0];
    }
    double C1[3];
    {
      dfm2::CVec3d dC1_dP[3][3];
      double dC1_dt[3][2];
      dfm2::CdC_Rod(
          C1, dC1_dP, dC1_dt,
          P1, S1);
    }
    for (unsigned int i = 0; i < 3; ++i) {
      const double val0 = (C1[i] - C0[i]) / eps;
      const double val1 = (
          +dC0_dt[i][0] * dT[0]
              + dC0_dt[i][1] * dT[1]
              + dC0_dP[i][0].dot(dP[0])
              + dC0_dP[i][1].dot(dP[1])
              + dC0_dP[i][2].dot(dP[2])) / eps;
      EXPECT_NEAR(val0, val1, 1.0e-2 * (1 + fabs(val1)));
    }
  }
}

TEST(fem_rod, CW_Rod) {
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  std::uniform_real_distribution<double> dist_12(+1, +2);
  const double eps = 1.0e-5;
  for (int itr = 0; itr < 100; ++itr) {
    dfm2::CVec3d P0[3], S0[2];
    bool res = GenRandomConfigRod(P0, S0, randomEng, dist_m1p1);
    if (!res) { continue; }
    double C[3];
    dfm2::CVec3d dC_dP[3][3];
    double dC_dt[3][2];
    dfm2::CdC_Rod(
        C, dC_dP, dC_dt,
        P0, S0);
    //
    const double stiff_bendtwist[3] = {
        dist_12(randomEng),
        dist_12(randomEng),
        dist_12(randomEng)};
    const double off[3] = {
        dist_m1p1(randomEng),
        dist_m1p1(randomEng),
        dist_m1p1(randomEng)};
    double W;
    dfm2::CVec3d dW_dP[3];
    double dW_dt[2];
    {
      dfm2::CMat3d ddW_ddP[3][3];
      dfm2::CVec3d ddW_dtdP[2][3];
      double ddW_ddt[2][2];
      W = WdWddW_Rod3Exact(dW_dP, dW_dt,
                     ddW_ddP, ddW_dtdP, ddW_ddt,
                     stiff_bendtwist, P0, S0,
                     dfm2::CVec3d(off));
    }
    EXPECT_NEAR(W,
                0.5 * (
                      stiff_bendtwist[0] * (C[0]-off[0]) * (C[0]-off[0])
                    + stiff_bendtwist[1] * (C[1]-off[1]) * (C[1]-off[1])
                    + stiff_bendtwist[2] * (C[2]-off[2]) * (C[2]-off[2])),
                1.0e-2 * (1 + W));
    for (int ino = 0; ino < 3; ++ino) {
      dfm2::CVec3d v0 = dW_dP[ino];
      dfm2::CVec3d v1 =
            stiff_bendtwist[0] * (C[0]-off[0]) * dC_dP[0][ino]
          + stiff_bendtwist[1] * (C[1]-off[1]) * dC_dP[1][ino]
          + stiff_bendtwist[2] * (C[2]-off[2]) * dC_dP[2][ino];
      EXPECT_LT((v0 - v1).norm(), v0.norm() * 1.0e-5);
    }
    for (int ino = 0; ino < 2; ++ino) {
      double v0 = dW_dt[ino];
      double v1 =
            stiff_bendtwist[0] * (C[0]-off[0]) * dC_dt[0][ino]
          + stiff_bendtwist[1] * (C[1]-off[1]) * dC_dt[1][ino]
          + stiff_bendtwist[2] * (C[2]-off[2]) * dC_dt[2][ino];
      EXPECT_NEAR(v0, v1, (fabs(v0) + fabs(v1)) * 1.0e-5);
    }
  }
}