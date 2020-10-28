/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femrod.h"
#include "delfem2/mat3.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/mshtopo.h"

namespace delfem2 {
namespace femrod{

/**
 * @brief add derivatie of dot( Frm0[i], Frm1[j] ) with respect to the 3 points and 2 rotations of the rod element
 */
DFM2_INLINE void AddDiff_DotFrameRod
(CVec3d dV_dP[3],
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
 const CVec3d dF1_dt[3])
{
  dV_dt[0] += c*Frm1[i1]*dF0_dt[i0];
  dV_dt[1] += c*Frm0[i0]*dF1_dt[i1];
  dV_dP[0] -= c*Frm1[i1]*dF0_dv[i0];
  dV_dP[1] += c*Frm1[i1]*dF0_dv[i0];
  dV_dP[1] -= c*Frm0[i0]*dF1_dv[i1];
  dV_dP[2] += c*Frm0[i0]*dF1_dv[i1];
}

DFM2_INLINE void AddDiffDiff_DotFrameRod
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
  {
    CMat3d ddW_ddv;
    CVec3d ddW_dvdt;
    double ddW_ddt;
    DifDifFrameRod(ddW_ddv, ddW_dvdt, ddW_ddt, i0, (P[1]-P[0]).Length(), F1[i1], F0);
    ddV_dtdP[0][0] += c*(-ddW_dvdt);
    ddV_dtdP[0][1] += c*(+ddW_dvdt - dF0_dt[i0]*dF1_dv[i1]);
    ddV_dtdP[0][2] += c*(+dF0_dt[i0]*dF1_dv[i1]);
    ddV_ddt[0][0] += c*ddW_ddt;
    ddV_ddt[0][1] += c*dF0_dt[i0]*dF1_dt[i1];
    const CMat3d T = dF0_dv[i0].Trans()*dF1_dv[i1];
    ddV_ddP[0][0] += c*ddW_ddv;
    ddV_ddP[0][1] += c*(-ddW_ddv + T);
    ddV_ddP[0][2] += c*(-T);
    ddV_ddP[1][0] += c*(-ddW_ddv);
    ddV_ddP[1][1] += c*(+ddW_ddv - T);
    ddV_ddP[1][2] += c*(+T);
  }
  {
    CMat3d ddW_ddv;
    CVec3d ddW_dvdt;
    double ddW_ddt;
    DifDifFrameRod(ddW_ddv, ddW_dvdt, ddW_ddt, i1, (P[2]-P[1]).Length(), F0[i0], F1);
    ddV_dtdP[1][0] += c*-dF1_dt[i1]*dF0_dv[i0];
    ddV_dtdP[1][1] += c*(-ddW_dvdt + dF1_dt[i1]*dF0_dv[i0]);
    ddV_dtdP[1][2] += c*+ddW_dvdt;
    ddV_ddt[1][0] += c*dF0_dt[i0]*dF1_dt[i1];
    ddV_ddt[1][1] += c*ddW_ddt;
    const CMat3d T = dF1_dv[i1].Trans()*dF0_dv[i0];
    ddV_ddP[1][0] += c*+T;
    ddV_ddP[1][1] += c*(+ddW_ddv - T);
    ddV_ddP[1][2] += c*(-ddW_ddv);
    ddV_ddP[2][0] += c*(-T);
    ddV_ddP[2][1] += c*(-ddW_ddv + T);
    ddV_ddP[2][2] += c*(+ddW_ddv);
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

DFM2_INLINE void AddOuterProduct_FrameRod
(CMat3d ddV_ddP[3][3],
 CVec3d ddV_dtdP[2][3],
 double ddV_ddt[2][2],
 //
 double c,
 const CVec3d dA_dP[3],
 const double dA_dt[2],
 const CVec3d dB_dP[3],
 const double dB_dt[2])
{
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      ddV_ddP[i][j] += c*Mat3_OuterProduct(dA_dP[i], dB_dP[j]);
    }
  }
  for(int i=0;i<2;++i){
    for(int j=0;j<3;++j){
      ddV_dtdP[i][j] += c*dA_dt[i]*dB_dP[j];
    }
  }
  for(int i=0;i<2;++i){
    for(int j=0;j<2;++j){
      ddV_ddt[i][j] += c*dA_dt[i]*dB_dt[j];
    }
  }
}

}
}


/**
 * @brief energy W and its derivative dW and second derivative ddW
 * where W = a^T R(dn) b(theta)
 */
DFM2_INLINE void delfem2::RodFrameTrans
(CVec3d frm[3],
 const CVec3d& S0,
 const CVec3d& V01,
 const CVec3d& du,
 double dtheta)
{
  //  std::cout << "      "  << S0.Length() << std::endl;
  assert( fabs(S0.Length() - 1.0) < 1.0e-3 );
  assert( fabs(S0*V01) < 1.0e-3 );
  const CVec3d U0 = V01.Normalize();
  const CVec3d T0 = U0^S0;
  frm[2] = (V01+du).Normalize();
  CMat3d R = Mat3_MinimumRotation(U0, frm[2]);
  frm[0] = R*(cos(dtheta)*S0 + sin(dtheta)*T0);
  frm[1] = R*(cos(dtheta)*T0 - sin(dtheta)*S0);
}


DFM2_INLINE void delfem2::DiffFrameRod
(CMat3d dF_dv[3],
 CVec3d dF_dt[3],
 //
 double l01,
 const CVec3d Frm[3])
{
  dF_dt[0] = +Frm[1];
  dF_dt[1] = -Frm[0];
  dF_dt[2].SetZero();
  dF_dv[0] = (-1.0/l01)*Mat3_OuterProduct(Frm[2],Frm[0]);
  dF_dv[1] = (-1.0/l01)*Mat3_OuterProduct(Frm[2],Frm[1]);
  dF_dv[2] = (+1.0/l01)*( Mat3_Identity(1.0) - Mat3_OuterProduct(Frm[2],Frm[2]) );
}

DFM2_INLINE void delfem2::DifDifFrameRod
(CMat3d& ddW_ddv,
 CVec3d& ddW_dvdt,
 double& ddW_dtt,
 //
 unsigned int iaxis,
 double l01,
 const CVec3d& Q,
 const CVec3d Frm[3])
{
  if( iaxis == 0 ){
    ddW_dtt = -Frm[0]*Q;
    ddW_dvdt = -(Q*Frm[2])*Frm[1]/l01;
  }
  else if( iaxis == 1 ){
    ddW_dtt = -Frm[1]*Q;
    ddW_dvdt = +(Q*Frm[2])*Frm[0]/l01;
  }
  else if( iaxis == 2 ){
    ddW_dtt = 0.0;
    ddW_dvdt = CVec3d(0,0,0);
  }
  {
    CMat3d S = Mat3_Spin(Frm[2]);
    CMat3d A = Mat3_Spin(Frm[iaxis])*Mat3_Spin(Q);
    CMat3d M0a = -S*(A*S);
    CVec3d b0 = (-A+A.Trans())*Frm[2];
    CMat3d M1 = Mat3_OuterProduct(Frm[2], b0);
    CMat3d M3 = (b0*Frm[2])*(3*Mat3_OuterProduct(Frm[2], Frm[2])-Mat3_Identity(1.0));
    ddW_ddv = (1.0/(l01*l01))*( M0a + M1 + M1.Trans() + M3 );
  }
}

DFM2_INLINE double delfem2::WdWddW_DotFrame
(CVec3d dV_dP[3],
 double dV_dt[2],
 CMat3d ddV_ddP[3][3],
 CVec3d ddV_dtdP[2][3],
 double ddV_ddt[2][2],
 //
 const CVec3d P[3],
 const CVec3d S[2],
 const double off[3])
{
  assert( fabs(S[0].Length() - 1.0) < 1.0e-10 );
  assert( fabs(S[0]*(P[1]-P[0]).Normalize()) < 1.0e-10 );
  assert( fabs(S[1].Length() - 1.0) < 1.0e-10 );
  assert( fabs(S[1]*(P[2]-P[1]).Normalize()) < 1.0e-10 );
  CVec3d Frm0[3];
  {
    Frm0[2] = (P[1]-P[0]).Normalize();
    Frm0[0] = S[0];
    Frm0[1] = Cross(Frm0[2],Frm0[0]);
  }
  CVec3d Frm1[3];
  {
    Frm1[2] = (P[2]-P[1]).Normalize();
    Frm1[0] = S[1];
    Frm1[1] = Cross(Frm1[2],Frm1[0]);
  }
  // ----------
  CMat3d dF0_dv[3];
  CVec3d dF0_dt[3];
  DiffFrameRod(dF0_dv, dF0_dt,
               (P[1]-P[0]).Length(), Frm0);
  CMat3d dF1_dv[3];
  CVec3d dF1_dt[3];
  DiffFrameRod(dF1_dv, dF1_dt,
               (P[2]-P[1]).Length(), Frm1);
  double V = 0;
  for(int i=0;i<3;++i){ dV_dP[i].SetZero(); }
  for(int i=0;i<2;++i){ dV_dt[i] = 0.0; }
  for(int i=0;i<4;++i){ (&ddV_ddt[0][0])[i] = 0.0; }
  for(int i=0;i<6;++i){ (&ddV_dtdP[0][0])[i].SetZero(); }
  for(int i=0;i<9;++i){ (&ddV_ddP[0][0])[i].SetZero(); }
  // ---------------------
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      const double c = i*3+j*5+7;
      V += c*Frm0[i]*Frm1[j];
      femrod::AddDiff_DotFrameRod(dV_dP, dV_dt,
                                c,
                                i, Frm0, dF0_dv, dF0_dt,
                                j, Frm1, dF1_dv, dF1_dt);
      femrod::AddDiffDiff_DotFrameRod(ddV_ddP, ddV_dtdP,ddV_ddt,
                                    c, P,
                                    i, Frm0, dF0_dv, dF0_dt,
                                    j, Frm1, dF1_dv, dF1_dt);
    }
  }
  return V;
}


DFM2_INLINE delfem2::CVec3d delfem2::Darboux_Rod
 (const CVec3d P[3],
  const CVec3d S[2])
{
  assert( fabs(S[0].Length() - 1.0) < 1.0e-5 );
  assert( fabs(S[0]*(P[1]-P[0]).Normalize()) < 1.0e-5 );
  assert( fabs(S[1].Length() - 1.0) < 1.0e-5 );
  assert( fabs(S[1]*(P[2]-P[1]).Normalize()) < 1.0e-5 );
  CVec3d F0[3];
  {
    F0[2] = (P[1]-P[0]).Normalize();
    F0[0] = S[0];
    F0[1] = Cross(F0[2],F0[0]);
  }
  CVec3d F1[3];
  {
    F1[2] = (P[2]-P[1]).Normalize();
    F1[0] = S[1];
    F1[1] = Cross(F1[2],F1[0]);
  }
  const double Y = 1 + F0[0]*F1[0] + F0[1]*F1[1] + F0[2]*F1[2];
  const double X[3] = {
    F0[1]*F1[2] - F0[2]*F1[1],
    F0[2]*F1[0] - F0[0]*F1[2],
    F0[0]*F1[1] - F0[1]*F1[0] };
  return CVec3d(X[0]/Y,X[1]/Y,X[2]/Y);
}

// TODO: stiffness
DFM2_INLINE double delfem2::WdWddW_Rod
(CVec3d dW_dP[3],
 double dW_dt[2],
 CMat3d ddW_ddP[3][3],
 CVec3d ddW_dtdP[2][3],
 double ddW_ddt[2][2],
 const double stiff_bendtwist[3],
 const CVec3d P[3],
 const CVec3d S[2],
 const CVec3d& darboux0,
 bool is_exact )
{
  assert( fabs(S[0].Length() - 1.0) < 1.0e-5 );
  assert( fabs(S[0]*(P[1]-P[0]).Normalize()) < 1.0e-5 );
  assert( fabs(S[1].Length() - 1.0) < 1.0e-5 );
  assert( fabs(S[1]*(P[2]-P[1]).Normalize()) < 1.0e-5 );
  CVec3d F0[3];
  {
    F0[2] = (P[1]-P[0]).Normalize();
    F0[0] = S[0];
    F0[1] = Cross(F0[2],F0[0]);
  }
  CVec3d F1[3];
  {
    F1[2] = (P[2]-P[1]).Normalize();
    F1[0] = S[1];
    F1[1] = Cross(F1[2],F1[0]);
  }
  // ----------
  CMat3d dF0_dv[3];
  CVec3d dF0_dt[3];
  DiffFrameRod(dF0_dv, dF0_dt,
               (P[1]-P[0]).Length(), F0);
  CMat3d dF1_dv[3];
  CVec3d dF1_dt[3];
  DiffFrameRod(dF1_dv, dF1_dt,
               (P[2]-P[1]).Length(), F1);
  // ------------
  for(int i=0;i<3;++i){ dW_dP[i].SetZero(); }
  for(int i=0;i<2;++i){ dW_dt[i] = 0.0; }
  for(int i=0;i<4;++i){ (&ddW_ddt[0][0])[i] = 0.0; }
  for(int i=0;i<6;++i){ (&ddW_dtdP[0][0])[i].SetZero(); }
  for(int i=0;i<9;++i){ (&ddW_ddP[0][0])[i].SetZero(); }
  // ------------
  const double Y = 1 + F0[0]*F1[0] + F0[1]*F1[1] + F0[2]*F1[2];
  CVec3d dY_dp[3];
  double dY_dt[2];
  { // making derivative of Y
    dY_dp[0].SetZero();
    dY_dp[1].SetZero();
    dY_dp[2].SetZero();
    dY_dt[0] = 0.0;
    dY_dt[1] = 0.0;
    femrod::AddDiff_DotFrameRod(dY_dp, dY_dt,
                              +1,
                              0, F0, dF0_dv, dF0_dt,
                              0, F1, dF1_dv, dF1_dt);
    femrod::AddDiff_DotFrameRod(dY_dp, dY_dt,
                              +1,
                              1, F0, dF0_dv, dF0_dt,
                              1, F1, dF1_dv, dF1_dt);
    femrod::AddDiff_DotFrameRod(dY_dp, dY_dt,
                              +1,
                              2, F0, dF0_dv, dF0_dt,
                              2, F1, dF1_dv, dF1_dt);
  }
  // ---------------------
  const double X[3] = {
    F0[1]*F1[2] - F0[2]*F1[1],
    F0[2]*F1[0] - F0[0]*F1[2],
    F0[0]*F1[1] - F0[1]*F1[0] };
  const double R[3] = {
    X[0]/Y-darboux0.x(),
    X[1]/Y-darboux0.y(),
    X[2]/Y-darboux0.z() };
  for(unsigned int iaxis=0;iaxis<3;++iaxis){
    const double stfa = stiff_bendtwist[iaxis];
    const unsigned int jaxis = (iaxis+1)%3;
    const unsigned int kaxis = (iaxis+2)%3;
    CVec3d dX_dp[3];
    double dX_dt[2];
    {
      dX_dp[0].SetZero();  dX_dp[1].SetZero();  dX_dp[2].SetZero();
      dX_dt[0] = 0.0;  dX_dt[1] = 0.0;
      femrod::AddDiff_DotFrameRod(dX_dp, dX_dt,
                                +1,
                                jaxis, F0, dF0_dv, dF0_dt,
                                kaxis, F1, dF1_dv, dF1_dt);
      femrod::AddDiff_DotFrameRod(dX_dp, dX_dt,
                                -1,
                                kaxis, F0, dF0_dv, dF0_dt,
                                jaxis, F1, dF1_dv, dF1_dt);
    }
    {
      CVec3d dR_dp[3];
      double dR_dt[2];
      {
        const double t0 = 1.0/Y;
        const double t1 = - X[iaxis]/(Y*Y);
        dR_dp[0] = t0*dX_dp[0] + t1*dY_dp[0];
        dR_dp[1] = t0*dX_dp[1] + t1*dY_dp[1];
        dR_dp[2] = t0*dX_dp[2] + t1*dY_dp[2];
        dR_dt[0] = t0*dX_dt[0] + t1*dY_dt[0];
        dR_dt[1] = t0*dX_dt[1] + t1*dY_dt[1];
      }
      dW_dP[0] += stfa*R[iaxis]*dR_dp[0];
      dW_dP[1] += stfa*R[iaxis]*dR_dp[1];
      dW_dP[2] += stfa*R[iaxis]*dR_dp[2];
      dW_dt[0] += stfa*R[iaxis]*dR_dt[0];
      dW_dt[1] += stfa*R[iaxis]*dR_dt[1];
      // [dR/dp][dR/dq]
      femrod::AddOuterProduct_FrameRod(ddW_ddP, ddW_dtdP, ddW_ddt,
                                       stfa, dR_dp, dR_dt, dR_dp, dR_dt);
    }
    { // -Ri/(Y*Y) { [dY/dq][dXi/dp] + [dY/dp][dXi/dq] }
      const double t0 = -stfa*R[iaxis]/(Y*Y);
      femrod::AddOuterProduct_FrameRod(ddW_ddP, ddW_dtdP, ddW_ddt,
                                       t0, dX_dp, dX_dt, dY_dp, dY_dt);
      femrod::AddOuterProduct_FrameRod(ddW_ddP, ddW_dtdP, ddW_ddt,
                                       t0, dY_dp, dY_dt, dX_dp, dX_dt);
    }
    // -------------
    if( is_exact ){ // (Ri/Y) [[dXi/dpdq]]
      const double t0 = stfa*R[iaxis]/Y;
      femrod::AddDiffDiff_DotFrameRod(ddW_ddP, ddW_dtdP,ddW_ddt,
                                    +t0,P,
                                    jaxis,F0, dF0_dv, dF0_dt,
                                    kaxis,F1, dF1_dv, dF1_dt);
      femrod::AddDiffDiff_DotFrameRod(ddW_ddP, ddW_dtdP,ddW_ddt,
                                    -t0,P,
                                    kaxis,F0, dF0_dv, dF0_dt,
                                    jaxis,F1, dF1_dv, dF1_dt);
    }
  }
  // ---------------
  if( is_exact ){ // -(R0*X0+R1*X1+R2*X2)/(Y*Y) [[ddY/dpdq]]
    const double stf0 = stiff_bendtwist[0];
    const double stf1 = stiff_bendtwist[1];
    const double stf2 = stiff_bendtwist[2];
    const double t0 = -(stf0*R[0]*X[0] + stf1*R[1]*X[1] + stf2*R[2]*X[2])/(Y*Y);
    femrod::AddDiffDiff_DotFrameRod(ddW_ddP, ddW_dtdP,ddW_ddt,
                                  t0,P,
                                  0,F0, dF0_dv, dF0_dt,
                                  0,F1, dF1_dv, dF1_dt);
    femrod::AddDiffDiff_DotFrameRod(ddW_ddP, ddW_dtdP,ddW_ddt,
                                  t0,P,
                                  1,F0,dF0_dv, dF0_dt,
                                  1,F1,dF1_dv, dF1_dt);
    femrod::AddDiffDiff_DotFrameRod(ddW_ddP, ddW_dtdP,ddW_ddt,
                                  t0,P,
                                  2,F0,dF0_dv, dF0_dt,
                                  2,F1,dF1_dv, dF1_dt);
  }
  { // 2*(R0*X0+R1*X1+R2*X2)/(Y*Y*Y) [dY/dp] * [dY/dq]
    const double stf0 = stiff_bendtwist[0];
    const double stf1 = stiff_bendtwist[1];
    const double stf2 = stiff_bendtwist[2];
    const double t0 = +(stf0*R[0]*X[0] + stf1*R[1]*X[1] + stf2*R[2]*X[2])*2.0/(Y*Y*Y);
    femrod::AddOuterProduct_FrameRod(ddW_ddP, ddW_dtdP, ddW_ddt,
                                   t0, dY_dp, dY_dt, dY_dp, dY_dt);
  }
  return
  0.5*stiff_bendtwist[0]*R[0]*R[0] +
  0.5*stiff_bendtwist[1]*R[1]*R[1] +
  0.5*stiff_bendtwist[2]*R[2]*R[2];
}



DFM2_INLINE double delfem2::WdWddW_SquareLengthLineseg3D(
    CVec3d dW_dP[2],
    CMat3d ddW_ddP[2][2],
    const double stiff,
    const CVec3d p[2],
    double L0)
{
  const double l  = sqrt(+ (p[0][0]-p[1][0])*(p[0][0]-p[1][0])
                         + (p[0][1]-p[1][1])*(p[0][1]-p[1][1])
                         + (p[0][2]-p[1][2])*(p[0][2]-p[1][2]));
  const CVec3d v = p[0]-p[1];
  const double R = L0-l;
  dW_dP[0] = (-R*stiff/l)*v;
  dW_dP[1] = (+R*stiff/l)*v;
  const CMat3d m = (stiff*L0/(l*l*l))*Mat3_OuterProduct(v, v) + (stiff*(l-L0)/l)*Mat3_Identity(1.0);
  ddW_ddP[0][0] = m;
  ddW_ddP[0][1] = -m;
  ddW_ddP[1][0] = -m;
  ddW_ddP[1][1] = m;
  return 0.5*stiff*R*R;
}


DFM2_INLINE void delfem2::Solve_DispRotSeparate
 (std::vector<CVec3d>& aP,
  std::vector<CVec3d>& aS,
  CMatrixSparse<double>& mats,
  const double stiff_stretch,
  const double stiff_bendtwist[3],
  const std::vector<CVec3d>& aP0,
  const std::vector<CVec3d>& aS0,
  const std::vector<unsigned int>& aElemSeg,
  const std::vector<unsigned int>& aElemRod,
  const std::vector<int>& aBCFlag)
{
  assert( mats.len_col == 3 );
  assert( mats.len_row == 3 );
  const unsigned int nNode = aBCFlag.size()/3;
  assert(aP.size()+aS.size()==nNode);
  mats.SetZero();
  std::vector<double> vec_r;
  vec_r.assign(nNode*3, 0.0);
  std::vector<int> tmp_buffer;
  double W = 0;
  for(unsigned int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    const unsigned int i0 = aElemSeg[iseg*2+0];
    const unsigned int i1 = aElemSeg[iseg*2+1];
    const unsigned int* aINoel = aElemSeg.data()+iseg*2;
    const double L0 = (aP0[i0]-aP0[i1]).Length();
    const CVec3d aPE[2] = { aP[i0], aP[i1] };
    // --------------
    CVec3d dW_dP[2];
    CMat3d ddW_ddP[2][2];
    W += WdWddW_SquareLengthLineseg3D(dW_dP, ddW_ddP,
                                      stiff_stretch, aPE, L0);
    {
      double eM[2*2*3*3];
      for(int in=0;in<2;++in){
        for(int jn=0;jn<2;++jn){
          ddW_ddP[in][jn].CopyTo(eM+(in*2+jn)*9);
        }
      }
      mats.Mearge(2, aINoel, 2, aINoel, 9, eM, tmp_buffer);
    }
    {
      for (int inoel=0; inoel<2; inoel++){
        const unsigned int ip = aINoel[inoel];
        vec_r[ip*3+0] -= dW_dP[inoel].x();
        vec_r[ip*3+1] -= dW_dP[inoel].y();
        vec_r[ip*3+2] -= dW_dP[inoel].z();
      }
    }
  }
  for(unsigned int irod=0;irod<aElemRod.size()/5;++irod){
    const unsigned int* aINoel = aElemRod.data()+irod*5;
    const unsigned int nP = aP.size();
    const CVec3d aPE[3] = { aP[aINoel[0]], aP[aINoel[1]], aP[aINoel[2]] };
    const CVec3d aSE[2] = { aS[aINoel[3]-nP], aS[aINoel[4]-nP] };
    CVec3d Darboux0;
    {
      const CVec3d aPE0[3] = { aP0[aINoel[0]], aP0[aINoel[1]], aP0[aINoel[2]] };
      const CVec3d aSE0[2] = { aS0[aINoel[3]-nP], aS0[aINoel[4]-nP] };
      Darboux0 = Darboux_Rod(aPE0, aSE0);
    }
    // ------
    CVec3d dW_dP[3];
    double dW_dt[2];
    CMat3d ddW_ddP[3][3];
    CVec3d ddW_dtdP[2][3];
    double ddW_ddt[2][2];
    W +=  WdWddW_Rod(dW_dP,dW_dt,ddW_ddP,ddW_dtdP,ddW_ddt,
                     stiff_bendtwist,
                     aPE,aSE, Darboux0, false);
    {
      double eM[5][5][3][3];
      for(int i=0;i<5*5*3*3;++i){ (&eM[0][0][0][0])[i] = 0.0; }
      for(int in=0;in<3;++in){
        for(int jn=0;jn<3;++jn){
          ddW_ddP[in][jn].CopyTo(&eM[in][jn][0][0]);
        }
      }
      for(int in=0;in<3;++in){
        for(int jn=0;jn<2;++jn){
          eM[3+jn][in][0][0] = eM[in][jn+3][0][0] = ddW_dtdP[jn][in].x();
          eM[3+jn][in][0][1] = eM[in][jn+3][1][0] = ddW_dtdP[jn][in].y();
          eM[3+jn][in][0][2] = eM[in][jn+3][2][0] = ddW_dtdP[jn][in].z();
        }
      }
      for(int in=0;in<2;++in){
        for(int jn=0;jn<2;++jn){
          eM[in+3][jn+3][0][0] = ddW_ddt[in][jn];
        }
      }
      mats.Mearge(5, aINoel, 5, aINoel, 9, &eM[0][0][0][0], tmp_buffer);
    }
    {
      for (int inoel=0; inoel<3; inoel++){
        const unsigned int ip = aINoel[inoel];
        vec_r[ip*3+0] -= dW_dP[inoel].x();
        vec_r[ip*3+1] -= dW_dP[inoel].y();
        vec_r[ip*3+2] -= dW_dP[inoel].z();
      }
      for (int inoel=0; inoel<2; inoel++){
        const unsigned int in0 = aINoel[3+inoel];
        vec_r[in0*3+0] -= dW_dt[inoel];
      }
    }
  }
  //  std::cout << CheckSymmetry(mats) << std::endl;
  //  mats.AddDia(0.00001);
  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << CheckSymmetry(mats) << std::endl;
  mats.SetFixedBC(aBCFlag.data());
  setRHS_Zero(vec_r, aBCFlag,0);
  std::vector<double> vec_x;
  vec_x.assign(nNode*3, 0.0);
  {
    auto aConvHist = Solve_CG(vec_r.data(),vec_x.data(),
                              vec_r.size(), 1.0e-4, 300, mats);
    if( aConvHist.size() > 0 ){
      std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
    }
  }
  /*
   {
   auto aConvHist = Solve_BiCGStab(vec_r,vec_x,
   1.0e-4, 300, mats);
   if( aConvHist.size() > 0 ){
   std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
   }
   }
   */
  assert( aS.size() == aElemSeg.size()/2 );
  for(unsigned int is=0;is<aS.size();++is){
    const unsigned int i0 = aElemSeg[is*2+0];
    const unsigned int i1 = aElemSeg[is*2+1];
    CVec3d V01 = aP[i1]-aP[i0];
    CVec3d du(vec_x[i1*3+0]-vec_x[i0*3+0],
                    vec_x[i1*3+1]-vec_x[i0*3+1],
                    vec_x[i1*3+2]-vec_x[i0*3+2]);
    const unsigned int np = aP.size();
    const double dtheta = vec_x[ np*3 + is*3 ];
    CVec3d frm[3];
    RodFrameTrans(frm,
                  aS[is],V01,du,dtheta);
    aS[is] = frm[0];
  }
  for(unsigned int ip=0;ip<aP.size();++ip){
    aP[ip].p[0] += vec_x[ip*3+0];
    aP[ip].p[1] += vec_x[ip*3+1];
    aP[ip].p[2] += vec_x[ip*3+2];
  }
  for(unsigned int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    const unsigned int i0 = aElemSeg[iseg*2+0];
    const unsigned int i1 = aElemSeg[iseg*2+1];
    const CVec3d& p0 = aP[i0];
    const CVec3d& p1 = aP[i1];
    const CVec3d e01 = (p1-p0).Normalize();
    aS[iseg] -= (aS[iseg]*e01)*e01;
    aS[iseg].SetNormalizedVector();
  }
}


// ========================================
// RodHair

DFM2_INLINE void delfem2::ParallelTransport_RodHair(
    std::vector<CVec3d>& aP0,
    std::vector<CVec3d>& aS0,
    const std::vector<unsigned int>& aIP_HairRoot)
{
  assert( aP0.size() == aS0.size() );
  assert( !aIP_HairRoot.empty() && aIP_HairRoot[0] == 0 );
  assert( aP0.size() == aIP_HairRoot[aIP_HairRoot.size()-1]);
  for(unsigned int ih=0;ih<aIP_HairRoot.size()-1;++ih){
    const unsigned int ip_r = aIP_HairRoot[ih];
    const unsigned int np = aIP_HairRoot[ih+1]-ip_r;
    for(unsigned int ir=0;ir<np-2;++ir){
      const unsigned int ip0 = ip_r+ir+0;
      const unsigned int ip1 = ip_r+ir+1;
      const unsigned int ip2 = ip_r+ir+2;
      const unsigned int is0 = ip_r+ir+0;
      const unsigned int is1 = ip_r+ir+1;
      const CMat3d CMat3 = Mat3_MinimumRotation(aP0[ip1]-aP0[ip0], aP0[ip2]-aP0[ip1]);
      CVec3d s1 = CMat3*aS0[is0] + aS0[is1];
      const CVec3d v = (aP0[ip2] - aP0[ip1]).Normalize();
      aS0[is1] = (s1-(s1*v)*v).Normalize();
    }
  }
}

DFM2_INLINE void delfem2::MakeBCFlag_RodHair(
    std::vector<int>& aBCFlag,
    const std::vector<unsigned int>& aIP_HairRoot)
{
  assert( !aIP_HairRoot.empty() && aIP_HairRoot[0] == 0 );
  const unsigned int np = aIP_HairRoot[aIP_HairRoot.size()-1];
  aBCFlag.assign(np*4, 0);
  for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
    assert( aIP_HairRoot[ihair+1] > aIP_HairRoot[ihair] );
    unsigned int ips0 = aIP_HairRoot[ihair]+0;
    unsigned int ips1 = aIP_HairRoot[ihair]+1;
    unsigned int ipe1 = aIP_HairRoot[ihair+1]-1;
    aBCFlag[ips0*4+0] = 1; aBCFlag[ips0*4+1] = 1; aBCFlag[ips0*4+2] = 1;
    aBCFlag[ips1*4+0] = 1; aBCFlag[ips1*4+1] = 1; aBCFlag[ips1*4+2] = 1;
    aBCFlag[ips0*4+3] = 1;
    aBCFlag[ipe1*4+3] = 1;
  }
}

DFM2_INLINE void delfem2::MakeSparseMatrix_RodHair(
    CMatrixSparse<double>& mats,
    const std::vector<unsigned int>& aIP_HairRoot)
{
  assert( !aIP_HairRoot.empty() && aIP_HairRoot[0] == 0 );
  const unsigned int np = aIP_HairRoot[aIP_HairRoot.size()-1];
  std::vector<unsigned int> psup_ind, psup;
  {
    std::vector<unsigned int> aElemRod;
    for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
      unsigned int ip0 = aIP_HairRoot[ihair];
      unsigned int nr = aIP_HairRoot[ihair+1]-ip0-2; // number of rod elements
      for (unsigned int ir = 0; ir < nr; ++ir) {
        aElemRod.push_back(ip0 + ir + 0);
        aElemRod.push_back(ip0 + ir + 1);
        aElemRod.push_back(ip0 + ir + 2);
      }
    }
    JArray_PSuP_MeshElem(psup_ind, psup,
                               aElemRod.data(), aElemRod.size()/3, 3, np);
  }
  JArray_Sort(psup_ind, psup);
  mats.Initialize(np, 4, true);
  mats.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
}


DFM2_INLINE void delfem2::MakeDirectorOrthogonal_RodHair(
    std::vector<CVec3d>& aS,
    const std::vector<CVec3d>& aP)
{
  for(unsigned int is=0;is<aP.size()-1;++is){
    assert( is < aS.size() );
    const unsigned int ip0 = is+0;
    const unsigned int ip1 = is+1;
    const CVec3d& p0 = aP[ip0];
    const CVec3d& p1 = aP[ip1];
    const CVec3d e01 = (p1-p0).Normalize();
    aS[is] -= (aS[is]*e01)*e01;
    aS[is].SetNormalizedVector();
  }
  /*
  for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
    unsigned int ips = aIP_HairRoot[ihair];
    unsigned int ns = aIP_HairRoot[ihair+1] - aIP_HairRoot[ihair] -1;
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0;
      const unsigned int ip1 = ips+is+1;
      const CVec3d& p0 = aP[ip0];
      const CVec3d& p1 = aP[ip1];
      const CVec3d e01 = (p1-p0).Normalize();
      aS[ip0] -= (aS[ip0]*e01)*e01;
      aS[ip0].SetNormalizedVector();
    }
  }
  */
}

double delfem2::MergeLinSys_Hair(
    std::vector<double>& vec_r,
    CMatrixSparse<double>& mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    const std::vector<unsigned int>& aIP_HairRoot,
    const std::vector<CVec3d>& aP,
    const std::vector<CVec3d>& aS,
    const std::vector<CVec3d>& aP0,
    const std::vector<CVec3d>& aS0)
{
  std::vector<int> tmp_buffer;
  double W = 0;
  for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
    unsigned int ips = aIP_HairRoot[ihair];
    unsigned int ns = aIP_HairRoot[ihair+1] - aIP_HairRoot[ihair] -1;
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0;
      const unsigned int ip1 = ips+is+1;
      const unsigned int aINoel[2] = {ip0,ip1};
      const double L0 = (aP0[ip0]-aP0[ip1]).Length();
      const CVec3d aPE[2] = { aP[ip0], aP[ip1] };
      // --------------
      CVec3d dW_dP[2];
      CMat3d ddW_ddP[2][2];
      W += WdWddW_SquareLengthLineseg3D(dW_dP, ddW_ddP,
                                        stiff_stretch, aPE, L0);
      {
        double eM[2*2*4*4]; for(auto& m : eM){ m = 0.0; }
        for(int in=0;in<2;++in){
          for(int jn=0;jn<2;++jn){
            ddW_ddP[in][jn].CopyToMat4( eM+(in*2+jn)*16 );
          }
        }
        mats.Mearge(2, aINoel, 2, aINoel, 16, eM, tmp_buffer);
      }
      for (int in=0; in<2; in++){
        const unsigned int ip = aINoel[in];
        vec_r[ip*4+0] -= dW_dP[in].x();
        vec_r[ip*4+1] -= dW_dP[in].y();
        vec_r[ip*4+2] -= dW_dP[in].z();
      }
    }
  }
  // --------------------------
  for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
    const unsigned int ips = aIP_HairRoot[ihair];
    const unsigned int nr = aIP_HairRoot[ihair+1] - aIP_HairRoot[ihair] - 2;
    for(unsigned int ir=0;ir<nr;++ir){
      const unsigned int ip0 = ips+ir+0;
      const unsigned int ip1 = ips+ir+1;
      const unsigned int ip2 = ips+ir+2;
      const unsigned int aINoel[3] = {ip0,ip1,ip2};
      const CVec3d aPE[3] = { aP[ip0], aP[ip1], aP[ip2] };
      const CVec3d aSE[2] = { aS[ip0], aS[ip1] };
      CVec3d Darboux0;
      {
        const CVec3d aPE0[3] = { aP0[ip0], aP0[ip1], aP0[ip2] };
        const CVec3d aSE0[2] = { aS0[ip0], aS0[ip1] };
        Darboux0 = Darboux_Rod(aPE0, aSE0);
      }
      // ------
      CVec3d dW_dP[3];
      double dW_dt[2];
      CMat3d ddW_ddP[3][3];
      CVec3d ddW_dtdP[2][3];
      double ddW_ddt[2][2];
      W +=  WdWddW_Rod(dW_dP,dW_dt,ddW_ddP,ddW_dtdP,ddW_ddt,
                       stiff_bendtwist,
                       aPE,aSE,Darboux0, false);
      {
        double eM[3*3*4*4]; for(auto& m: eM){ m = 0.0; }
        for(int in=0;in<3;++in){
          for(int jn=0;jn<3;++jn){
            ddW_ddP[in][jn].CopyToMat4(eM+(in*3+jn)*16);
          }
        }
        for(int in=0;in<3;++in){
          for(int jn=0;jn<2;++jn){
            eM[(in*3+jn)*16+0*4+3] = eM[(jn*3+in)*16+3*4+0] = ddW_dtdP[jn][in].x();
            eM[(in*3+jn)*16+1*4+3] = eM[(jn*3+in)*16+3*4+1] = ddW_dtdP[jn][in].y();
            eM[(in*3+jn)*16+2*4+3] = eM[(jn*3+in)*16+3*4+2] = ddW_dtdP[jn][in].z();
          }
        }
        for(int in=0;in<2;++in){
          for(int jn=0;jn<2;++jn){
            eM[(in*3+jn)*16+4*3+3] = ddW_ddt[in][jn];
          }
        }
        mats.Mearge(3, aINoel, 3, aINoel, 16, eM, tmp_buffer);
      }
      {
        for (int ino=0; ino<3; ino++){
          const unsigned int ip = aINoel[ino];
          vec_r[ip*4+0] -= dW_dP[ino].x();
          vec_r[ip*4+1] -= dW_dP[ino].y();
          vec_r[ip*4+2] -= dW_dP[ino].z();
        }
        for (int in=0; in<2; in++){
          const unsigned int in0 = aINoel[in];
          vec_r[in0*4+3] -= dW_dt[in];
        }
      }
    }
  }
  return W;
}

DFM2_INLINE void delfem2::UpdateSolutionHair(
    std::vector<CVec3d>& aP,
    std::vector<CVec3d>& aS,
    const std::vector<double>& vec_x,
    const std::vector<unsigned int>& aIP_HairRoot,
    const std::vector<int>& aBCFlag)
{
  for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
    unsigned int ips = aIP_HairRoot[ihair];
    unsigned int ns = aIP_HairRoot[ihair+1] - aIP_HairRoot[ihair] -1;
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0;
      const unsigned int ip1 = ips+is+1;
      CVec3d V01 = aP[ip1]-aP[ip0];
      CVec3d du(vec_x[ip1*4+0]-vec_x[ip0*4+0],
                vec_x[ip1*4+1]-vec_x[ip0*4+1],
                vec_x[ip1*4+2]-vec_x[ip0*4+2]);
      const double dtheta = vec_x[ip0*4+3];
      CVec3d frm[3];
      RodFrameTrans(frm,
                    aS[ip0],V01,du,dtheta);
      aS[ip0] = frm[0];
    }
  }
  for(unsigned int ip=0;ip<aP.size();++ip){
    if( aBCFlag[ip*4+0] != 0 ) continue; 
    aP[ip].p[0] += vec_x[ip*4+0];
    aP[ip].p[1] += vec_x[ip*4+1];
    aP[ip].p[2] += vec_x[ip*4+2];
  }
}

DFM2_INLINE void delfem2::Solve_RodHair(
    std::vector<CVec3d>& aP,
    std::vector<CVec3d>& aS,
    CMatrixSparse<double>& mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d>& aP0,
    const std::vector<CVec3d>& aS0,
    const std::vector<int>& aBCFlag,
    const std::vector<unsigned int>& aIP_HairRoot)
{
  assert( mats.len_col == 4 );
  assert( mats.len_row == 4 );
  const unsigned int np = aP.size();
  assert( aP0.size() == np );
  assert( aS0.size() == np );
  assert( aP.size() == np );
  assert( aS.size() == np );
  mats.SetZero();
  std::vector<double> vec_r;
  vec_r.assign(np*4, 0.0);
  double W = MergeLinSys_Hair(
      vec_r,mats,
      stiff_stretch,stiff_bendtwist,
      aIP_HairRoot,aP,aS,aP0,aS0);
  for(unsigned int ip=0;ip<aP.size();++ip){
    mats.valDia[ip*16+0*4+0] += mdtt;
    mats.valDia[ip*16+1*4+1] += mdtt;
    mats.valDia[ip*16+2*4+2] += mdtt;
  }
  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << CheckSymmetry(mats) << std::endl;
  assert( aBCFlag.size() == np*4 );
  mats.SetFixedBC(aBCFlag.data());
  setRHS_Zero(vec_r, aBCFlag,0);
  std::vector<double> vec_x;
  vec_x.assign(np*4, 0.0);
  {
    auto aConvHist = Solve_CG(vec_r.data(),vec_x.data(),
                              vec_r.size(), 1.0e-4, 300, mats);
    if( aConvHist.size() > 0 ){
      std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
    }
  }
  UpdateSolutionHair(aP,aS,
      vec_x,aIP_HairRoot,aBCFlag);
}

namespace delfem2 {
namespace femrod {
  class CMatContact {
  public:
    CMatContact(
        const CMatrixSparse<double>& mats_,
        const std::vector<CContactHair>& aContact_,
        const double stiffness_)
        : mats(mats_), aContact(aContact_), stiffness(stiffness_)
        {}
    void MatVec(double *y,
                double alpha, const double *x,
                double beta) const{
      mats.MatVec(y,
                  alpha,x,beta);
      for(unsigned int ic=0;ic<aContact.size();++ic){
        const CContactHair& ch = aContact[ic];
        const unsigned int aIP[4] = {ch.ip0, ch.ip1, ch.iq0, ch.iq1};
        const double aW[4] = {1-ch.s, ch.s, -(1-ch.t), -ch.t};
        const CVec3d& nrm = ch.norm;
//        std::cout << "lennorm: " << nrm.Length() << std::endl;
        const auto NN = Mat3_OuterProduct(nrm,nrm);
        for(int iip=0;iip<4;++iip){
          unsigned int ip0 = aIP[iip];
          for(int jjp=0;jjp<4;++jjp){
            unsigned int jp0 = aIP[jjp];
            double v0[3];
            MatVec3(v0,
                NN.mat, x+jp0*4);
            y[ip0*4+0] += alpha*stiffness*aW[iip]*aW[jjp]*v0[0];
            y[ip0*4+1] += alpha*stiffness*aW[iip]*aW[jjp]*v0[1];
            y[ip0*4+2] += alpha*stiffness*aW[iip]*aW[jjp]*v0[2];
          }
        }
      }
    }
  public:
    const CMatrixSparse<double>& mats;
    const std::vector<CContactHair>& aContact;
    const double stiffness;
  };
}
}


DFM2_INLINE void delfem2::Solve_RodHairContact(
    std::vector<CVec3d>& aP,
    std::vector<CVec3d>& aS,
    CMatrixSparse<double>& mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d>& aPt0,
    const std::vector<CVec3d>& aP0,
    const std::vector<CVec3d>& aS0,
    const std::vector<int>& aBCFlag,
    const std::vector<unsigned int>& aIP_HairRoot,
    const double clearance,
    const double stiff_contact,
    const std::vector<CContactHair>& aContact)
{
  assert( mats.len_col == 4 );
  assert( mats.len_row == 4 );
  const unsigned int np = aP.size();
  assert( aP0.size() == np );
  assert( aS0.size() == np );
  assert( aP.size() == np );
  assert( aS.size() == np );
  mats.SetZero();
  std::vector<double> vec_r;
  vec_r.assign(np*4, 0.0);
  double W = MergeLinSys_Hair(
      vec_r,mats,
      stiff_stretch,stiff_bendtwist,
      aIP_HairRoot,aP,aS,aP0,aS0);
  for(unsigned int ip=0;ip<aP.size();++ip){
    mats.valDia[ip*16+0*4+0] += mdtt;
    mats.valDia[ip*16+1*4+1] += mdtt;
    mats.valDia[ip*16+2*4+2] += mdtt;
  }
  for(unsigned int ip=0;ip<aP.size();++ip){ // this term has effect only in the second nonliear iteration
    const CVec3d dp = aPt0[ip]-aP[ip];
    vec_r[ip*4+0] += dp.x()*mdtt;
    vec_r[ip*4+1] += dp.y()*mdtt;
    vec_r[ip*4+2] += dp.z()*mdtt;
  }
//  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << CheckSymmetry(mats) << std::endl;
  for(const auto& ch : aContact){
    const unsigned int aIP[4] = {ch.ip0, ch.ip1, ch.iq0, ch.iq1};
    const double aW[4] = {1-ch.s, ch.s, -(1-ch.t), -ch.t};
    CVec3d a = aW[0]*aP[aIP[0]] + aW[1]*aP[aIP[1]] + aW[2]*aP[aIP[2]] + aW[3]*aP[aIP[3]];
    double r0 = a*ch.norm - clearance;
//    std::cout << "    contact: " << r0 << std::endl;
    for(int iip=0;iip<4;++iip){
      const unsigned int ip0 = aIP[iip];
      vec_r[ip0*4+0] -= stiff_contact*r0*aW[iip]*ch.norm.x();
      vec_r[ip0*4+1] -= stiff_contact*r0*aW[iip]*ch.norm.y();
      vec_r[ip0*4+2] -= stiff_contact*r0*aW[iip]*ch.norm.z();
    }
  }
  // --------------
  W = 0.0; // to remove warning
  assert( aBCFlag.size() == np*4 );
  mats.SetFixedBC(aBCFlag.data());
  std::vector<double> vec_x;
  vec_x.assign(np*4, 0.0);
  setRHS_Zero(vec_r, aBCFlag,0);
  femrod::CMatContact mc(mats,aContact,stiff_contact);
  {
    auto aConvHist = Solve_CG(vec_r.data(),vec_x.data(),
                              vec_r.size(), 1.0e-6, 3000, mc);
    /*
    if( aConvHist.size() > 0 ){
      std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
    }
     */
//    std::cout << "  updates: " << DotX(vec_x.data(), vec_x.data(), vec_x.size()) << std::endl;
  }
  UpdateSolutionHair(aP,aS,
                     vec_x,aIP_HairRoot,aBCFlag);
  /*
  std::cout << "hogehoge " << aContact.size() << std::endl;
  for(const auto& ch : aContact){
    const unsigned int aIP[4] = {ch.ip0, ch.ip1, ch.iq0, ch.iq1};
    const double aW[4] = {1-ch.s, ch.s, -(1-ch.t), -ch.t};
    CVec3d a = aW[0]*aP[aIP[0]] + aW[1]*aP[aIP[1]] + aW[2]*aP[aIP[2]] + aW[3]*aP[aIP[3]];
    double r0 = a*ch.norm - clearance;
    std::cout << "       contact: " << r0 << " " << clearance << std::endl;
  }
   */
}

