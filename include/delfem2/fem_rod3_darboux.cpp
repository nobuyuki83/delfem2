/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/fem_rod3_darboux.h"

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/mat3.h"
#include "delfem2/mshuni.h"

namespace delfem2::femrod {

/**
 * @brief add derivatie of dot( Frm0[i], Frm1[j] ) with respect to the 3 points and 2 rotations of the rod element
 */
DFM2_INLINE void AddDiff_DotFrameRod(
    CVec3d dV_dP[3],
    double dV_dt[2],
    //
    double c,
    unsigned int i0,
    const CVec3d Frm0[3],
    const CMat3d dF0_dv[3],
    const CVec3d dF0_dt[3],
    unsigned int i1,
    const CVec3d Frm1[3],
    const CMat3d dF1_dv[3],
    const CVec3d dF1_dt[3]) {
  dV_dt[0] += c * Frm1[i1].dot(dF0_dt[i0]);
  dV_dt[1] += c * Frm0[i0].dot(dF1_dt[i1]);
  dV_dP[0] -= c * Frm1[i1] * dF0_dv[i0];
  dV_dP[1] += c * Frm1[i1] * dF0_dv[i0];
  dV_dP[1] -= c * Frm0[i0] * dF1_dv[i1];
  dV_dP[2] += c * Frm0[i0] * dF1_dv[i1];
}

DFM2_INLINE void AddDiffDiff_DotFrameRod(
    CMat3d ddV_ddP[3][3],
    CVec3d ddV_dtdP[2][3],
    double ddV_ddt[2][2],
    //
    double c,
    const CVec3d P[3],
    unsigned int i0,
    const CVec3d F0[3],
    const CMat3d dF0_dv[3],
    const CVec3d dF0_dt[3],
    unsigned int i1,
    const CVec3d F1[3],
    const CMat3d dF1_dv[3],
    const CVec3d dF1_dt[3]) {
  {
    CMat3d ddW_ddv;
    CVec3d ddW_dvdt;
    double ddW_ddt;
    DifDifFrameRod(ddW_ddv, ddW_dvdt, ddW_ddt, i0, (P[1] - P[0]).norm(), F1[i1], F0);
    ddV_dtdP[0][0] += c * (-ddW_dvdt);
    ddV_dtdP[0][1] += c * (+ddW_dvdt - dF0_dt[i0] * dF1_dv[i1]);
    ddV_dtdP[0][2] += c * (+dF0_dt[i0] * dF1_dv[i1]);
    ddV_ddt[0][0] += c * ddW_ddt;
    ddV_ddt[0][1] += c * dF0_dt[i0].dot(dF1_dt[i1]);
    const CMat3d T = dF0_dv[i0].transpose() * dF1_dv[i1];
    ddV_ddP[0][0] += c * ddW_ddv;
    ddV_ddP[0][1] += c * (-ddW_ddv + T);
    ddV_ddP[0][2] += c * (-T);
    ddV_ddP[1][0] += c * (-ddW_ddv);
    ddV_ddP[1][1] += c * (+ddW_ddv - T);
    ddV_ddP[1][2] += c * (+T);
  }
  {
    CMat3d ddW_ddv;
    CVec3d ddW_dvdt;
    double ddW_ddt;
    DifDifFrameRod(ddW_ddv, ddW_dvdt, ddW_ddt, i1, (P[2] - P[1]).norm(), F0[i0], F1);
    ddV_dtdP[1][0] += c * -dF1_dt[i1] * dF0_dv[i0];
    ddV_dtdP[1][1] += c * (-ddW_dvdt + dF1_dt[i1] * dF0_dv[i0]);
    ddV_dtdP[1][2] += c * +ddW_dvdt;
    ddV_ddt[1][0] += c * dF0_dt[i0].dot(dF1_dt[i1]);
    ddV_ddt[1][1] += c * ddW_ddt;
    const CMat3d T = dF1_dv[i1].transpose() * dF0_dv[i0];
    ddV_ddP[1][0] += c * +T;
    ddV_ddP[1][1] += c * (+ddW_ddv - T);
    ddV_ddP[1][2] += c * (-ddW_ddv);
    ddV_ddP[2][0] += c * (-T);
    ddV_ddP[2][1] += c * (-ddW_ddv + T);
    ddV_ddP[2][2] += c * (+ddW_ddv);
  }
}

/*
 DFM2_INLINE void AddDiffDiff_DotFrameRodSym
 (CMat3d ddV_ddP[3][3],
 CVec3d ddV_dtdP[2][3],
 double ddV_ddt[2][2],
 //
 double c,
 const CVec3d P[3],
 unsigned int i0,
 const CVec3d F0[3],
 const CMat3d dF0_dv[3],
 const CVec3d dF0_dt[3],
 unsigned int i1,
 const CVec3d F1[3],
 const CMat3d dF1_dv[3],
 const CVec3d dF1_dt[3])
 {
 CMat3d ddV0_ddP[3][3];
 double ddV0_ddt[2][2] = {{0,0},{0,0}};
 AddDiffDiff_DotFrameRod(ddV0_ddP, ddV_dtdP, ddV0_ddt,
 c, P,
 i0, F0, dF0_dv, dF0_dt,
 i1, F1, dF1_dv, dF1_dt);
 for(int i=0;i<2;++i){
 for(int j=0;j<2;++j){
 ddV_ddt[i][j] += 0.5*(ddV0_ddt[i][j] + ddV0_ddt[j][i]);
 //      ddV_ddt[i][j] += ddV0_ddt[i][j];
 }
 }
 for(int i=0;i<3;++i){
 for(int j=0;j<3;++j){
 //      ddV_ddP[i][j] += ddV0_ddP[i][j];
 ddV_ddP[i][j] += 0.5*(ddV0_ddP[i][j] + ddV0_ddP[j][i].Trans());
 }
 }
 }
 */

DFM2_INLINE void AddOuterProduct_FrameRod(
    CMat3d ddV_ddP[3][3],
    CVec3d ddV_dtdP[2][3],
    double ddV_ddt[2][2],
    //
    double c,
    const CVec3d dA_dP[3],
    const double dA_dt[2],
    const CVec3d dB_dP[3],
    const double dB_dt[2]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ddV_ddP[i][j] += c * Mat3_OuterProduct(dA_dP[i], dB_dP[j]);
    }
  }
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      ddV_dtdP[i][j] += c * dA_dt[i] * dB_dP[j];
    }
  }
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      ddV_ddt[i][j] += c * dA_dt[i] * dB_dt[j];
    }
  }
}

}  // namespace delfem2


// -------------------------------------------------

/**
 * @brief energy W and its derivative dW and second derivative ddW
 * where W = a^T R(dn) b(theta)
 */
DFM2_INLINE void delfem2::RodFrameTrans(
    CVec3d frm[3],
    const CVec3d &S0,
    const CVec3d &V01,
    const CVec3d &du,
    double dtheta) {
  //  std::cout << "      "  << S0.Length() << std::endl;
  assert(fabs(S0.norm() - 1.0) < 1.0e-3);
  assert(fabs(S0.dot(V01)) < 1.0e-3);
  const CVec3d &U0 = V01.normalized();
  const CVec3d &T0 = U0 ^ S0;
  frm[2] = (V01 + du).normalized();
  const CMat3d &R = Mat3_MinimumRotation(U0, frm[2]);
  frm[0] = R * (cos(dtheta) * S0 + sin(dtheta) * T0);
  frm[1] = R * (cos(dtheta) * T0 - sin(dtheta) * S0);
}

DFM2_INLINE void delfem2::DiffFrameRod(
    CMat3d dF_dv[3],
    CVec3d dF_dt[3],
    //
    double l01,
    const CVec3d Frm[3]) {
  dF_dt[0] = +Frm[1];
  dF_dt[1] = -Frm[0];
  dF_dt[2].setZero();
  dF_dv[0] = (-1.0 / l01) * Mat3_OuterProduct(Frm[2], Frm[0]);
  dF_dv[1] = (-1.0 / l01) * Mat3_OuterProduct(Frm[2], Frm[1]);
  dF_dv[2] = (+1.0 / l01) * (Mat3_Identity(1.0) - Mat3_OuterProduct(Frm[2], Frm[2]));
}

DFM2_INLINE void delfem2::DifDifFrameRod(
    CMat3d &ddW_ddv,
    CVec3d &ddW_dvdt,
    double &ddW_dtt,
    //
    unsigned int iaxis,
    double l01,
    const CVec3d &Q,
    const CVec3d Frm[3]) {
  if (iaxis == 0) {
    ddW_dtt = -Frm[0].dot(Q);
    ddW_dvdt = -(Q.dot(Frm[2])) * Frm[1] / l01;
  } else if (iaxis == 1) {
    ddW_dtt = -Frm[1].dot(Q);
    ddW_dvdt = +(Q.dot(Frm[2])) * Frm[0] / l01;
  } else if (iaxis == 2) {
    ddW_dtt = 0.0;
    ddW_dvdt = CVec3d(0, 0, 0);
  }
  {
    CMat3d S = Mat3_Spin(Frm[2]);
    CMat3d A = Mat3_Spin(Frm[iaxis]) * Mat3_Spin(Q);
    CMat3d M0a = -S * (A * S);
    CVec3d b0 = (-A + A.transpose()) * Frm[2];
    CMat3d M1 = Mat3_OuterProduct(Frm[2], b0);
    CMat3d M3 = (b0.dot(Frm[2])) * (3 * Mat3_OuterProduct(Frm[2], Frm[2]) - Mat3_Identity(1.0));
    ddW_ddv = (1.0 / (l01 * l01)) * (M0a + M1 + M1.transpose() + M3);
  }
}

DFM2_INLINE double delfem2::WdWddW_DotFrame(
    CVec3d dV_dP[3],
    double dV_dt[2],
    CMat3d ddV_ddP[3][3],
    CVec3d ddV_dtdP[2][3],
    double ddV_ddt[2][2],
    //
    const CVec3d P[3],
    const CVec3d S[2],
    [[maybe_unused]] const double off[3]) {
  assert(fabs(S[0].norm() - 1.0) < 1.0e-10);
  assert(fabs(S[0].dot((P[1] - P[0]).normalized())) < 1.0e-10);
  assert(fabs(S[1].norm() - 1.0) < 1.0e-10);
  assert(fabs(S[1].dot((P[2] - P[1]).normalized())) < 1.0e-10);
  CVec3d Frm0[3];
  {
    Frm0[2] = (P[1] - P[0]).normalized();
    Frm0[0] = S[0];
    Frm0[1] = Cross(Frm0[2], Frm0[0]);
  }
  CVec3d Frm1[3];
  {
    Frm1[2] = (P[2] - P[1]).normalized();
    Frm1[0] = S[1];
    Frm1[1] = Cross(Frm1[2], Frm1[0]);
  }
  // ----------
  CMat3d dF0_dv[3];
  CVec3d dF0_dt[3];
  DiffFrameRod(dF0_dv, dF0_dt,
               (P[1] - P[0]).norm(), Frm0);
  CMat3d dF1_dv[3];
  CVec3d dF1_dt[3];
  DiffFrameRod(dF1_dv, dF1_dt,
               (P[2] - P[1]).norm(), Frm1);
  double V = 0;
  for (int i = 0; i < 3; ++i) { dV_dP[i].setZero(); }
  for (int i = 0; i < 2; ++i) { dV_dt[i] = 0.0; }
  for (int i = 0; i < 4; ++i) { (&ddV_ddt[0][0])[i] = 0.0; }
  for (int i = 0; i < 6; ++i) { (&ddV_dtdP[0][0])[i].setZero(); }
  for (int i = 0; i < 9; ++i) { (&ddV_ddP[0][0])[i].setZero(); }
  // ---------------------
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double c = i * 3 + j * 5 + 7;
      V += c * Frm0[i].dot(Frm1[j]);
      femrod::AddDiff_DotFrameRod(dV_dP, dV_dt,
                                  c,
                                  i, Frm0, dF0_dv, dF0_dt,
                                  j, Frm1, dF1_dv, dF1_dt);
      femrod::AddDiffDiff_DotFrameRod(ddV_ddP, ddV_dtdP, ddV_ddt,
                                      c, P,
                                      i, Frm0, dF0_dv, dF0_dt,
                                      j, Frm1, dF1_dv, dF1_dt);
    }
  }
  return V;
}

DFM2_INLINE delfem2::CVec3d delfem2::Darboux_Rod
    (const CVec3d P[3],
     const CVec3d S[2]) {
  assert(fabs(S[0].norm() - 1.0) < 1.0e-5);
  assert(fabs(S[0].dot((P[1] - P[0]).normalized())) < 1.0e-5);
  assert(fabs(S[1].norm() - 1.0) < 1.0e-5);
  assert(fabs(S[1].dot((P[2] - P[1]).normalized())) < 1.0e-5);
  CVec3d F0[3];
  { // frame on line segment 01
    F0[2] = (P[1] - P[0]).normalized();
    F0[0] = S[0];
    F0[1] = Cross(F0[2], F0[0]);
  }
  CVec3d F1[3];
  { // frame on line segment 12
    F1[2] = (P[2] - P[1]).normalized();
    F1[0] = S[1];
    F1[1] = Cross(F1[2], F1[0]);
  }
  const double Y = 1 + F0[0].dot(F1[0]) + F0[1].dot(F1[1]) + F0[2].dot(F1[2]);
  const double X[3] = {
      F0[1].dot(F1[2]) - F0[2].dot(F1[1]),
      F0[2].dot(F1[0]) - F0[0].dot(F1[2]),
      F0[0].dot(F1[1]) - F0[1].dot(F1[0])};
  return CVec3d(X[0] / Y, X[1] / Y, X[2] / Y);
}

// TODO: stiffness
DFM2_INLINE double delfem2::WdWddW_Rod(
    CVec3d dW_dP[3],
    double dW_dt[2],
    CMat3d ddW_ddP[3][3],
    CVec3d ddW_dtdP[2][3],
    double ddW_ddt[2][2],
    const double stiff_bendtwist[3],
    const CVec3d P[3],
    const CVec3d S[2],
    const CVec3d &darboux0,
    bool is_exact) {
  assert(fabs(S[0].norm() - 1.0) < 1.0e-5);
  assert(fabs(S[0].dot((P[1] - P[0]).normalized())) < 1.0e-5);
  assert(fabs(S[1].norm() - 1.0) < 1.0e-5);
  assert(fabs(S[1].dot((P[2] - P[1]).normalized())) < 1.0e-5);
  CVec3d F[3];
  { // compute frame on 01 segment
    F[2] = (P[1] - P[0]).normalized();
    F[0] = S[0];
    F[1] = Cross(F[2], F[0]);
  }
  CVec3d G[3];
  { // compute frame on 12 segment
    G[2] = (P[2] - P[1]).normalized();
    G[0] = S[1];
    G[1] = Cross(G[2], G[0]);
  }
  // ----------
  CMat3d dF_dv[3];
  CVec3d dF_dt[3];
  DiffFrameRod(
      dF_dv, dF_dt,
      (P[1] - P[0]).norm(), F);
  CMat3d dG_dv[3];
  CVec3d dG_dt[3];
  DiffFrameRod(
      dG_dv, dG_dt,
      (P[2] - P[1]).norm(), G);
  // ------------
  for (int i = 0; i < 3; ++i) { dW_dP[i].setZero(); }
  for (int i = 0; i < 2; ++i) { dW_dt[i] = 0.0; }
  for (int i = 0; i < 4; ++i) { (&ddW_ddt[0][0])[i] = 0.0; }
  for (int i = 0; i < 6; ++i) { (&ddW_dtdP[0][0])[i].setZero(); }
  for (int i = 0; i < 9; ++i) { (&ddW_ddP[0][0])[i].setZero(); }
  // ------------
  const double Y = 1 + F[0].dot(G[0]) + F[1].dot(G[1]) + F[2].dot(G[2]);
  CVec3d dY_dp[3]; // how Y changes w.r.t. the position
  double dY_dt[2]; // how Y changes w.r.t. the twist
  { // making derivative of Y
    dY_dp[0].setZero();
    dY_dp[1].setZero();
    dY_dp[2].setZero();
    dY_dt[0] = 0.0;
    dY_dt[1] = 0.0;
    femrod::AddDiff_DotFrameRod(
        dY_dp, dY_dt,
        +1,
        0, F, dF_dv, dF_dt,
        0, G, dG_dv, dG_dt);
    femrod::AddDiff_DotFrameRod(
        dY_dp, dY_dt,
        +1,
        1, F, dF_dv, dF_dt,
        1, G, dG_dv, dG_dt);
    femrod::AddDiff_DotFrameRod(
        dY_dp, dY_dt,
        +1,
        2, F, dF_dv, dF_dt,
        2, G, dG_dv, dG_dt);
  }
  // ---------------------
  const double X[3] = {
      F[1].dot(G[2]) - F[2].dot(G[1]),
      F[2].dot(G[0]) - F[0].dot(G[2]),
      F[0].dot(G[1]) - F[1].dot(G[0])};
  const double R[3] = {
      X[0] / Y - darboux0.x,
      X[1] / Y - darboux0.y,
      X[2] / Y - darboux0.z};
  for (unsigned int iaxis = 0; iaxis < 3; ++iaxis) {
    const double stfa = stiff_bendtwist[iaxis];
    const unsigned int jaxis = (iaxis + 1) % 3;
    const unsigned int kaxis = (iaxis + 2) % 3;
    CVec3d dX_dp[3];
    double dX_dt[2];
    {
      dX_dp[0].setZero();
      dX_dp[1].setZero();
      dX_dp[2].setZero();
      dX_dt[0] = 0.0;
      dX_dt[1] = 0.0;
      femrod::AddDiff_DotFrameRod(
          dX_dp, dX_dt,
          +1,
          jaxis, F, dF_dv, dF_dt,
          kaxis, G, dG_dv, dG_dt);
      femrod::AddDiff_DotFrameRod(
          dX_dp, dX_dt,
          -1,
          kaxis, F, dF_dv, dF_dt,
          jaxis, G, dG_dv, dG_dt);
    }
    {
      CVec3d dR_dp[3];
      double dR_dt[2];
      {
        const double t0 = 1.0 / Y;
        const double t1 = -X[iaxis] / (Y * Y);
        dR_dp[0] = t0 * dX_dp[0] + t1 * dY_dp[0];
        dR_dp[1] = t0 * dX_dp[1] + t1 * dY_dp[1];
        dR_dp[2] = t0 * dX_dp[2] + t1 * dY_dp[2];
        dR_dt[0] = t0 * dX_dt[0] + t1 * dY_dt[0];
        dR_dt[1] = t0 * dX_dt[1] + t1 * dY_dt[1];
      }
      dW_dP[0] += stfa * R[iaxis] * dR_dp[0];
      dW_dP[1] += stfa * R[iaxis] * dR_dp[1];
      dW_dP[2] += stfa * R[iaxis] * dR_dp[2];
      dW_dt[0] += stfa * R[iaxis] * dR_dt[0];
      dW_dt[1] += stfa * R[iaxis] * dR_dt[1];
      // [dR/dp][dR/dq]
      femrod::AddOuterProduct_FrameRod(
          ddW_ddP, ddW_dtdP, ddW_ddt,
          stfa, dR_dp, dR_dt, dR_dp, dR_dt);
    }
    { // -Ri/(Y*Y) { [dY/dq][dXi/dp] + [dY/dp][dXi/dq] }
      const double t0 = -stfa * R[iaxis] / (Y * Y);
      femrod::AddOuterProduct_FrameRod(
          ddW_ddP, ddW_dtdP, ddW_ddt,
          t0, dX_dp, dX_dt, dY_dp, dY_dt);
      femrod::AddOuterProduct_FrameRod(
          ddW_ddP, ddW_dtdP, ddW_ddt,
          t0, dY_dp, dY_dt, dX_dp, dX_dt);
    }
    // -------------
    if (is_exact) { // (Ri/Y) [[dXi/dpdq]]
      const double t0 = stfa * R[iaxis] / Y;
      femrod::AddDiffDiff_DotFrameRod(
          ddW_ddP, ddW_dtdP, ddW_ddt,
          +t0, P,
          jaxis, F, dF_dv, dF_dt,
          kaxis, G, dG_dv, dG_dt);
      femrod::AddDiffDiff_DotFrameRod(
          ddW_ddP, ddW_dtdP, ddW_ddt,
          -t0, P,
          kaxis, F, dF_dv, dF_dt,
          jaxis, G, dG_dv, dG_dt);
    }
  }
  // ---------------
  if (is_exact) { // -(R0*X0+R1*X1+R2*X2)/(Y*Y) [[ddY/dpdq]]
    const double stf0 = stiff_bendtwist[0];
    const double stf1 = stiff_bendtwist[1];
    const double stf2 = stiff_bendtwist[2];
    const double t0 = -(stf0 * R[0] * X[0] + stf1 * R[1] * X[1] + stf2 * R[2] * X[2]) / (Y * Y);
    femrod::AddDiffDiff_DotFrameRod(
        ddW_ddP, ddW_dtdP, ddW_ddt,
        t0, P,
        0, F, dF_dv, dF_dt,
        0, G, dG_dv, dG_dt);
    femrod::AddDiffDiff_DotFrameRod(
        ddW_ddP, ddW_dtdP, ddW_ddt,
        t0, P,
        1, F, dF_dv, dF_dt,
        1, G, dG_dv, dG_dt);
    femrod::AddDiffDiff_DotFrameRod(
        ddW_ddP, ddW_dtdP, ddW_ddt,
        t0, P,
        2, F, dF_dv, dF_dt,
        2, G, dG_dv, dG_dt);
  }
  { // 2*(R0*X0+R1*X1+R2*X2)/(Y*Y*Y) [dY/dp] * [dY/dq]
    const double stf0 = stiff_bendtwist[0];
    const double stf1 = stiff_bendtwist[1];
    const double stf2 = stiff_bendtwist[2];
    const double t0 = +(stf0 * R[0] * X[0] + stf1 * R[1] * X[1] + stf2 * R[2] * X[2]) * 2.0 / (Y * Y * Y);
    femrod::AddOuterProduct_FrameRod(
        ddW_ddP, ddW_dtdP, ddW_ddt,
        t0, dY_dp, dY_dt, dY_dp, dY_dt);
  }
  return
      0.5 * stiff_bendtwist[0] * R[0] * R[0] +
          0.5 * stiff_bendtwist[1] * R[1] * R[1] +
          0.5 * stiff_bendtwist[2] * R[2] * R[2];
}
