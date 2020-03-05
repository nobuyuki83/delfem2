#include <cstdlib>
#include <vector>
#include <set>
#include <random>
#include "delfem2/mat3.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
//
#include "delfem2/v23m3q.h"

// --------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_v23dtricad.h"
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------

/**
 * @brief energy W and its derivative dW and second derivative ddW
 * where W = a^T R(dn) b(theta)
 */
void RodFrameTrans(dfm2::CVec3d frm[3],
                   const dfm2::CVec3d& S0,
                   const dfm2::CVec3d& V01,
                   const dfm2::CVec3d& du,
                   double dtheta)
{
//  std::cout << "      "  << S0.Length() << std::endl;
  assert( fabs(S0.Length() - 1.0) < 1.0e-3 );
  assert( fabs(S0*V01) < 1.0e-3 );
  const dfm2::CVec3d U0 = V01.Normalize();
  const dfm2::CVec3d T0 = U0^S0;
  frm[2] = (V01+du).Normalize();
  dfm2::CMat3d R = dfm2::Mat3_MinimumRotation(U0, frm[2]);
  frm[0] = R*(cos(dtheta)*S0 + sin(dtheta)*T0);
  frm[1] = R*(cos(dtheta)*T0 - sin(dtheta)*S0);
}

void DiffFrameRod(
    dfm2::CMat3d dF_dv[3], // first-order derivative
    dfm2::CVec3d dF_dt[3],
    //
    double l01,
    const dfm2::CVec3d Frm[3])
{
  dF_dt[0] = +Frm[1];
  dF_dt[1] = -Frm[0];
  dF_dt[2].SetZero();
  dF_dv[0] = (-1.0/l01)*dfm2::Mat3_OuterProduct(Frm[2],Frm[0]);
  dF_dv[1] = (-1.0/l01)*dfm2::Mat3_OuterProduct(Frm[2],Frm[1]);
  dF_dv[2] = (+1.0/l01)*(dfm2::Mat3_Identity(1.0) - dfm2::Mat3_OuterProduct(Frm[2], Frm[2]) );
}

/**
 * @brief energy W and its derivative dW and second derivative ddW
 * where W = a^T R(dn) b(theta)
 */
void DifDifFrameRod(
    dfm2::CMat3d& ddW_ddv,
    dfm2::CVec3d& ddW_dvdt, // second-order derrivative
    double& ddW_dtt,
    //
    unsigned int iaxis,
    double l01,
    const dfm2::CVec3d& Q,
    const dfm2::CVec3d Frm[3])
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
    ddW_dvdt = dfm2::CVec3d(0,0,0);
  }
  {
    dfm2::CMat3d S = dfm2::Mat3_Spin(Frm[2]);
    dfm2::CMat3d A = dfm2::Mat3_Spin(Frm[iaxis])*dfm2::Mat3_Spin(Q);
    dfm2::CMat3d M0a = -S*(A*S);
    dfm2::CVec3d b0 = (-A+A.Trans())*Frm[2];
    dfm2::CMat3d M1 = dfm2::Mat3_OuterProduct(Frm[2], b0);
    dfm2::CMat3d M3 = (b0*Frm[2])*(3*dfm2::Mat3_OuterProduct(Frm[2], Frm[2])-dfm2::Mat3_Identity(1.0));
    ddW_ddv = (1.0/(l01*l01))*( M0a + M1 + M1.Trans() + M3 );
  }
}

void Check_WdWddW_RodFrameTrans(){
  dfm2::CVec3d V01;
  V01.SetRandom();
//  V01.SetNormalizedVector();
  dfm2::CVec3d Frm[3];
  {
    Frm[2] = V01;
    Frm[2].SetNormalizedVector();
    Frm[0].SetRandom();
    Frm[0] -= (Frm[0]*Frm[2])*Frm[2];
    Frm[0].SetNormalizedVector();
    Frm[1] = (Frm[2]^Frm[0]);
  }
  dfm2::CVec3d Q; Q.SetRandom();
//  Q = Frm[2];
  // --------------------------------
  double W[3] = { Q*Frm[0], Q*Frm[1], Q*Frm[2] };
  dfm2::CVec3d DW_Dv[3];
  double DW_Dt[3];
  {
    dfm2::CMat3d dF_dv[3];
    dfm2::CVec3d dF_dt[3];
    DiffFrameRod(dF_dv, dF_dt,
                 V01.Length(), Frm);
    DW_Dv[0] = Q*dF_dv[0];
    DW_Dv[1] = Q*dF_dv[1];
    DW_Dv[2] = Q*dF_dv[2];
    DW_Dt[0] = Q*dF_dt[0];
    DW_Dt[1] = Q*dF_dt[1];
    DW_Dt[2] = Q*dF_dt[2];
  }
  // ---------
  const double eps = 1.0e-6;
  dfm2::CVec3d du; du.SetRandom(); du *= eps;
  const double dtheta = (2.0*rand()/(RAND_MAX+1.0)-1.0)*eps;
  // ------
  dfm2::CVec3d frm[3];
  RodFrameTrans(frm,
                Frm[0], V01,
                du, dtheta);
  const double w[3] = { Q*frm[0], Q*frm[1], Q*frm[2] };
  dfm2::CVec3d dw_dv[3];
  double dw_dt[3];
  {
    dfm2::CMat3d df_dv[3];
    dfm2::CVec3d df_dt[3];
    DiffFrameRod(df_dv, df_dt,
                 (V01+du).Length(), frm);
    dw_dv[0] = Q*df_dv[0];
    dw_dv[1] = Q*df_dv[1];
    dw_dv[2] = Q*df_dv[2];
    dw_dt[0] = Q*df_dt[0];
    dw_dt[1] = Q*df_dt[1];
    dw_dt[2] = Q*df_dt[2];
  }
  for(int i=0;i<3;++i){
    double val0 = (w[i]-W[i])/eps;
    double val1 = (DW_Dt[i]*dtheta+DW_Dv[i]*du)/eps;
    std::cout << "dtv: " << i << " " << fabs(val0-val1) << "    ####  " << val0 << " " << val1 << std::endl;
  }
  //
  for(int i=0;i<3;++i){
    dfm2::CMat3d DDW_DDv;
    dfm2::CVec3d DDW_DvDt;
    double DDW_DDt;
    DifDifFrameRod(DDW_DDv, DDW_DvDt, DDW_DDt,
                   i,V01.Length(),Q,Frm);
    double val0 = (dw_dt[i]-DW_Dt[i])/eps;
    double val1 = (DDW_DDt*dtheta+DDW_DvDt*du)/eps;
    std::cout << "ddt: " << i << " " << fabs(val0-val1) << "    ####  " << val0 << " " << val1 << std::endl;
  }
  for(int i=0;i<3;++i){
    dfm2::CMat3d DDW_DDv;
    dfm2::CVec3d DDW_DvDt;
    double DDW_DDt;
    DifDifFrameRod(DDW_DDv, DDW_DvDt, DDW_DDt,
                   i,V01.Length(),Q,Frm);
    dfm2::CVec3d vec0 = (dw_dv[i]-DW_Dv[i])/eps;
    dfm2::CVec3d vec1 = (DDW_DvDt*dtheta+DDW_DDv*du)/eps;
    std::cout << "ddv: " << i << " " << (vec0-vec1).Length() << "    ####  " << vec0 << "   " << vec1 << std::endl;
  }
}


void AddDiff_DotFrames
 (dfm2::CVec3d dV_dP[3],
  double dV_dt[2],
  //
  double c,
  unsigned int i,
  unsigned int j,
  const dfm2::CVec3d Frm0[3],
  const dfm2::CVec3d Frm1[3],
  const dfm2::CMat3d dF0_dv[3],
  const dfm2::CVec3d dF0_dt[3],
  const dfm2::CMat3d dF1_dv[3],
  const dfm2::CVec3d dF1_dt[3])
{
  dV_dt[0] += c*Frm1[j]*dF0_dt[i];
  dV_dt[1] += c*Frm0[i]*dF1_dt[j];
  dV_dP[0] -= c*Frm1[j]*dF0_dv[i];
  dV_dP[1] += c*Frm1[j]*dF0_dv[i];
  dV_dP[1] -= c*Frm0[i]*dF1_dv[j];
  dV_dP[2] += c*Frm0[i]*dF1_dv[j];
}

void AddDiffDiff_DotFrames
(dfm2::CMat3d ddV_ddP[3][3],
 dfm2::CVec3d ddV_dtdP[2][3],
 double ddV_ddt[2][2],
 //
 double c,
 unsigned int i,
 unsigned int j,
 const dfm2::CVec3d P[3],
 const dfm2::CVec3d F0[3],
 const dfm2::CVec3d F1[3],
 const dfm2::CMat3d dF0_dv[3],
 const dfm2::CVec3d dF0_dt[3],
 const dfm2::CMat3d dF1_dv[3],
 const dfm2::CVec3d dF1_dt[3])
{
  {
    dfm2::CMat3d ddW_ddv;
    dfm2::CVec3d ddW_dvdt;
    double ddW_ddt;
    DifDifFrameRod(ddW_ddv, ddW_dvdt, ddW_ddt, i, (P[1]-P[0]).Length(), F1[j], F0);
    ddV_dtdP[0][0] += c*(-ddW_dvdt);
    ddV_dtdP[0][1] += c*(+ddW_dvdt - dF0_dt[i]*dF1_dv[j]);
    ddV_dtdP[0][2] += c*(+dF0_dt[i]*dF1_dv[j]);
    ddV_ddt[0][0] += c*ddW_ddt;
    ddV_ddt[0][1] += c*dF0_dt[i]*dF1_dt[j];
    const dfm2::CMat3d T = dF0_dv[i].Trans()*dF1_dv[j];
    ddV_ddP[0][0] += c*ddW_ddv;
    ddV_ddP[0][1] += c*(-ddW_ddv + T);
    ddV_ddP[0][2] += c*(-T);
    ddV_ddP[1][0] += c*(-ddW_ddv);
    ddV_ddP[1][1] += c*(+ddW_ddv - T);
    ddV_ddP[1][2] += c*(+T);
  }
  {
    dfm2::CMat3d ddW_ddv;
    dfm2::CVec3d ddW_dvdt;
    double ddW_ddt;
    DifDifFrameRod(ddW_ddv, ddW_dvdt, ddW_ddt, j, (P[2]-P[1]).Length(), F0[i], F1);
    ddV_dtdP[1][0] += c*-dF1_dt[j]*dF0_dv[i];
    ddV_dtdP[1][1] += c*(-ddW_dvdt + dF1_dt[j]*dF0_dv[i]);
    ddV_dtdP[1][2] += c*+ddW_dvdt;
    ddV_ddt[1][0] += c*dF0_dt[i]*dF1_dt[j];
    ddV_ddt[1][1] += c*ddW_ddt;
    const dfm2::CMat3d T = dF1_dv[j].Trans()*dF0_dv[i];
    ddV_ddP[1][0] += c*+T;
    ddV_ddP[1][1] += c*(+ddW_ddv - T);
    ddV_ddP[1][2] += c*(-ddW_ddv);
    ddV_ddP[2][0] += c*(-T);
    ddV_ddP[2][1] += c*(-ddW_ddv + T);
    ddV_ddP[2][2] += c*(+ddW_ddv);
  }
}


double WdWddW_DotFrame(
     dfm2::CVec3d dV_dP[3],
     double dV_dt[2],
     dfm2::CMat3d ddV_ddP[3][3],
     dfm2::CVec3d ddV_dtdP[2][3],
     double ddV_ddt[2][2],
     //
     const dfm2::CVec3d P[3],
     const dfm2::CVec3d S[2],
     const double off[3])
{
  assert( fabs(S[0].Length() - 1.0) < 1.0e-10 );
  assert( fabs(S[0]*(P[1]-P[0]).Normalize()) < 1.0e-10 );
  assert( fabs(S[1].Length() - 1.0) < 1.0e-10 );
  assert( fabs(S[1]*(P[2]-P[1]).Normalize()) < 1.0e-10 );
  dfm2::CVec3d Frm0[3];
  {
    Frm0[2] = (P[1]-P[0]).Normalize();
    Frm0[0] = S[0];
    Frm0[1] = dfm2::Cross(Frm0[2],Frm0[0]);
  }
  dfm2::CVec3d Frm1[3];
  {
    Frm1[2] = (P[2]-P[1]).Normalize();
    Frm1[0] = S[1];
    Frm1[1] = dfm2::Cross(Frm1[2],Frm1[0]);
  }
  // ----------
  dfm2::CMat3d dF0_dv[3];
  dfm2::CVec3d dF0_dt[3];
  DiffFrameRod(dF0_dv, dF0_dt,
               (P[1]-P[0]).Length(), Frm0);
  dfm2::CMat3d dF1_dv[3];
  dfm2::CVec3d dF1_dt[3];
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
    AddDiff_DotFrames(dV_dP, dV_dt,
                      c, i, j, Frm0, Frm1, dF0_dv, dF0_dt, dF1_dv, dF1_dt);
    AddDiffDiff_DotFrames(ddV_ddP, ddV_dtdP,ddV_ddt,
                          c, i, j, P,Frm0, Frm1, dF0_dv, dF0_dt,dF1_dv,dF1_dt);
  }
  }
  return V;
}

void Check_WdWddW_DotFrame()
{
  dfm2::CVec3d P[3];
  P[0].SetRandom();
  P[1].SetRandom();
  P[2].SetRandom();
  dfm2::CVec3d S[2];
  {
    S[0].SetRandom();
    const dfm2::CVec3d U0 = (P[1]-P[0]).Normalize();
    S[0] -= (S[0]*U0)*U0;
    S[0].SetNormalizedVector();
  }
  {
    S[1].SetRandom();
    const dfm2::CVec3d U1 = (P[2]-P[1]).Normalize();
    S[1] -= (S[1]*U1)*U1;
    S[1].SetNormalizedVector();
  }
  const double off[3] = {
    2.0*rand()/(RAND_MAX+1.0)-1.0,
    2.0*rand()/(RAND_MAX+1.0)-1.0,
    2.0*rand()/(RAND_MAX+1.0)-1.0 };
  // ------------------------
  dfm2::CVec3d dW_dP[3];
  double dW_dt[2];
  dfm2::CMat3d ddW_ddP[3][3];
  dfm2::CVec3d ddW_dtdP[2][3];
  double ddW_ddt[2][2];
  double W = WdWddW_DotFrame(dW_dP,dW_dt,
                              ddW_ddP, ddW_dtdP,ddW_ddt,
                              P, S, off);
  // -----------------------
  double eps = 1.0e-7;
  dfm2::CVec3d dP[3];
  dP[0].SetRandom(); dP[0] *= eps;
  dP[1].SetRandom(); dP[1] *= eps;
  dP[2].SetRandom(); dP[2] *= eps;
  const double dT[2] = {
    (2.0*rand()/(RAND_MAX+1.0)-1.0)*eps,
    (2.0*rand()/(RAND_MAX+1.0)-1.0)*eps };
  dfm2::CVec3d frm0[3], frm1[3];
  RodFrameTrans(frm0,
                S[0], P[1]-P[0], dP[1]-dP[0], dT[0]);
  RodFrameTrans(frm1,
                S[1], P[2]-P[1], dP[2]-dP[1], dT[1]);
  const dfm2::CVec3d p[3] = { P[0] + dP[0], P[1] + dP[1], P[2] + dP[2] };
  const dfm2::CVec3d s[2] = { frm0[0], frm1[0] };
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
    const double val0 = (w-W)/eps;
    const double val1 = (+dW_dt[0]*dT[0]
                         +dW_dt[1]*dT[1]
                         +dW_dP[0]*dP[0]
                         +dW_dP[1]*dP[1]
                         +dW_dP[2]*dP[2])/eps;
    std::cout << "diff_W : " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const double val0 = (dw_dt[0]-dW_dt[0])/eps;
    const double val1 = (+ddW_ddt[ 0][0]*dT[0]
                         +ddW_ddt[ 0][1]*dT[1]
                         +ddW_dtdP[0][0]*dP[0]
                         +ddW_dtdP[0][1]*dP[1]
                         +ddW_dtdP[0][2]*dP[2])/eps;
    std::cout << "diff_T0: " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const double val0 = (dw_dt[1]-dW_dt[1])/eps;
    const double val1 = (+ddW_ddt[ 1][0]*dT[0]
                         +ddW_ddt[ 1][1]*dT[1]
                         +ddW_dtdP[1][0]*dP[0]
                         +ddW_dtdP[1][1]*dP[1]
                         +ddW_dtdP[1][2]*dP[2])/eps;
    std::cout << "diff_T1: " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[0]-dW_dP[0])/eps;
    const dfm2::CVec3d val1 = (+ddW_dtdP[0][0]*dT[0]
                               +ddW_dtdP[1][0]*dT[1]
                               +ddW_ddP[0][0]*dP[0]
                               +ddW_ddP[0][1]*dP[1]
                               +ddW_ddP[0][2]*dP[2])/eps;
    std::cout << "diff_V0: " << (val0-val1).Length() << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[1]-dW_dP[1])/eps;
    const dfm2::CVec3d val1 = (+ddW_dtdP[0][1]*dT[0]
                               +ddW_dtdP[1][1]*dT[1]
                               +ddW_ddP[1][0]*dP[0]
                               +ddW_ddP[1][1]*dP[1]
                               +ddW_ddP[1][2]*dP[2])/eps;
    std::cout << "diff_V2: " << (val0-val1).Length() << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[2]-dW_dP[2])/eps;
    const dfm2::CVec3d val1 = (+ddW_dtdP[0][2]*dT[0]
                               +ddW_dtdP[1][2]*dT[1]
                               +ddW_ddP[2][0]*dP[0]
                               +ddW_ddP[2][1]*dP[1]
                               +ddW_ddP[2][2]*dP[2])/eps;
    std::cout << "diff_V2: " << (val0-val1).Length() << " ---> " << val0 << " " << val1 << std::endl;
  }
  std::cout << std::endl;
}


void AddOuterProduct(
    dfm2::CMat3d ddV_ddP[3][3],
    dfm2::CVec3d ddV_dtdP[2][3],
    double ddV_ddt[2][2],
    //
    double c,
    const dfm2::CVec3d dA_dP[3],
    const double dA_dt[2],
    const dfm2::CVec3d dB_dP[3],
    const double dB_dt[2])
{
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      ddV_ddP[i][j] += c*dfm2::Mat3_OuterProduct(dA_dP[i], dB_dP[j]);
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

double WdWddW_Rod(
    dfm2::CVec3d dV_dP[3],
    double dV_dt[2],
    dfm2::CMat3d ddV_ddP[3][3],
    dfm2::CVec3d ddV_dtdP[2][3],
    double ddV_ddt[2][2],
    //
    const dfm2::CVec3d P[3],
    const dfm2::CVec3d S[2],
    const double off[3],
    bool is_exact)
{
//  std::cout << S[0].Length() << std::endl;
  assert( fabs(S[0].Length() - 1.0) < 1.0e-5 );
  assert( fabs(S[0]*(P[1]-P[0]).Normalize()) < 1.0e-5 );
  assert( fabs(S[1].Length() - 1.0) < 1.0e-5 );
  assert( fabs(S[1]*(P[2]-P[1]).Normalize()) < 1.0e-5 );
  dfm2::CVec3d F0[3];
  {
    F0[2] = (P[1]-P[0]).Normalize();
    F0[0] = S[0];
    F0[1] = dfm2::Cross(F0[2],F0[0]);
  }
  dfm2::CVec3d F1[3];
  {
    F1[2] = (P[2]-P[1]).Normalize();
    F1[0] = S[1];
    F1[1] = dfm2::Cross(F1[2],F1[0]);
  }
  // ----------
  dfm2::CMat3d dF0_dv[3];
  dfm2::CVec3d dF0_dt[3];
  DiffFrameRod(dF0_dv, dF0_dt,
               (P[1]-P[0]).Length(), F0);
  dfm2::CMat3d dF1_dv[3];
  dfm2::CVec3d dF1_dt[3];
  DiffFrameRod(dF1_dv, dF1_dt,
               (P[2]-P[1]).Length(), F1);
  for(int i=0;i<3;++i){ dV_dP[i].SetZero(); }
  for(int i=0;i<2;++i){ dV_dt[i] = 0.0; }
  for(int i=0;i<4;++i){ (&ddV_ddt[0][0])[i] = 0.0; }
  for(int i=0;i<6;++i){ (&ddV_dtdP[0][0])[i].SetZero(); }
  for(int i=0;i<9;++i){ (&ddV_ddP[0][0])[i].SetZero(); }
  const double Y = 1 + F0[0]*F1[0] + F0[1]*F1[1] + F0[2]*F1[2];
  dfm2::CVec3d dY_dp[3]; double dY_dt[2];
  {
    dY_dp[0].SetZero();  dY_dp[1].SetZero();  dY_dp[2].SetZero();
    dY_dt[0] = 0.0;  dY_dt[1] = 0.0;
    AddDiff_DotFrames(dY_dp, dY_dt, +1, 0, 0, F0, F1, dF0_dv, dF0_dt, dF1_dv, dF1_dt);
    AddDiff_DotFrames(dY_dp, dY_dt, +1, 1, 1, F0, F1, dF0_dv, dF0_dt, dF1_dv, dF1_dt);
    AddDiff_DotFrames(dY_dp, dY_dt, +1, 2, 2, F0, F1, dF0_dv, dF0_dt, dF1_dv, dF1_dt);
  }
  // ---------------------
  const double X[3] = {
    F0[1]*F1[2] - F0[2]*F1[1],
    F0[2]*F1[0] - F0[0]*F1[2],
    F0[0]*F1[1] - F0[1]*F1[0] };
  const double R[3] = {
    X[0]/Y-off[0],
    X[1]/Y-off[1],
    X[2]/Y-off[2] };
  for(unsigned int iaxis=0;iaxis<3;++iaxis){
    const unsigned int jaxis = (iaxis+1)%3;
    const unsigned int kaxis = (iaxis+2)%3;
    dfm2::CVec3d dX_dP[3]; double dX_dt[2];
    {
      dX_dP[0].SetZero();  dX_dP[1].SetZero();  dX_dP[2].SetZero();
      dX_dt[0] = 0.0;  dX_dt[1] = 0.0;
      AddDiff_DotFrames(dX_dP, dX_dt, +1, jaxis, kaxis, F0, F1, dF0_dv, dF0_dt, dF1_dv, dF1_dt);
      AddDiff_DotFrames(dX_dP, dX_dt, -1, kaxis, jaxis, F0, F1, dF0_dv, dF0_dt, dF1_dv, dF1_dt);
    }
    dfm2::CVec3d dR0_P[3]; double dR0_t[2];
    {
      const double t0 = 1.0/Y;
      const double t1 = - X[iaxis]/(Y*Y);
      dR0_P[0] = t0*dX_dP[0] + t1*dY_dp[0];
      dR0_P[1] = t0*dX_dP[1] + t1*dY_dp[1];
      dR0_P[2] = t0*dX_dP[2] + t1*dY_dp[2];
      dR0_t[0] = t0*dX_dt[0] + t1*dY_dt[0];
      dR0_t[1] = t0*dX_dt[1] + t1*dY_dt[1];
    }
    {
      dV_dP[0] += R[iaxis]*dR0_P[0];
      dV_dP[1] += R[iaxis]*dR0_P[1];
      dV_dP[2] += R[iaxis]*dR0_P[2];
      dV_dt[0] += R[iaxis]*dR0_t[0];
      dV_dt[1] += R[iaxis]*dR0_t[1];
    }
    AddOuterProduct(ddV_ddP, ddV_dtdP, ddV_ddt,
                    1.0, dR0_P, dR0_t, dR0_P, dR0_t);
    if( is_exact ){
      double t0 = R[iaxis]/Y;
      AddDiffDiff_DotFrames(ddV_ddP, ddV_dtdP,ddV_ddt,
                            +t0, jaxis, kaxis, P,F0, F1, dF0_dv, dF0_dt,dF1_dv,dF1_dt);
      AddDiffDiff_DotFrames(ddV_ddP, ddV_dtdP,ddV_ddt,
                            -t0, kaxis, jaxis, P,F0, F1, dF0_dv, dF0_dt,dF1_dv,dF1_dt);
    }
    {
      double t0 = -R[iaxis]/(Y*Y);
      AddOuterProduct(ddV_ddP, ddV_dtdP, ddV_ddt,
                      t0, dX_dP, dX_dt, dY_dp, dY_dt);
      AddOuterProduct(ddV_ddP, ddV_dtdP, ddV_ddt,
                      t0, dY_dp, dY_dt, dX_dP, dX_dt);
    }
  }
  // ---------------
  if( is_exact ){
    double t0 = -(R[0]*X[0]+R[1]*X[1]+R[2]*X[2])/(Y*Y);
    AddDiffDiff_DotFrames(ddV_ddP, ddV_dtdP,ddV_ddt,
                          t0, 0, 0, P,F0, F1, dF0_dv, dF0_dt,dF1_dv,dF1_dt);
    AddDiffDiff_DotFrames(ddV_ddP, ddV_dtdP,ddV_ddt,
                          t0, 1, 1, P,F0, F1, dF0_dv, dF0_dt,dF1_dv,dF1_dt);
    AddDiffDiff_DotFrames(ddV_ddP, ddV_dtdP,ddV_ddt,
                          t0, 2, 2, P,F0, F1, dF0_dv, dF0_dt,dF1_dv,dF1_dt);
  }
  {
    double t0 = +(R[0]*X[0]+R[1]*X[1]+R[2]*X[2])*2.0/(Y*Y*Y);
    AddOuterProduct(ddV_ddP, ddV_dtdP, ddV_ddt,
                    t0, dY_dp, dY_dt, dY_dp, dY_dt);
  }
  return 0.5*(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
}


void Check_WdWddW_Rod()
{
  dfm2::CVec3d P[3];
  P[0].SetRandom();
  P[1].SetRandom();
  P[2].SetRandom();
  dfm2::CVec3d S[2];
  {
    S[0].SetRandom();
    const dfm2::CVec3d U0 = (P[1]-P[0]).Normalize();
    S[0] -= (S[0]*U0)*U0;
    S[0].SetNormalizedVector();
  }
  {
    S[1].SetRandom();
    const dfm2::CVec3d U1 = (P[2]-P[1]).Normalize();
    S[1] -= (S[1]*U1)*U1;
    S[1].SetNormalizedVector();
  }
  const double off[3] = {
    2.0*rand()/(RAND_MAX+1.0)-1.0,
    2.0*rand()/(RAND_MAX+1.0)-1.0,
    2.0*rand()/(RAND_MAX+1.0)-1.0 };
  // ------------------------
  dfm2::CVec3d dW_dP[3];
  double dW_dt[2];
  dfm2::CMat3d ddW_ddP[3][3];
  dfm2::CVec3d ddW_dtdP[2][3];
  double ddW_ddt[2][2];
  double W = WdWddW_Rod(dW_dP,dW_dt,
                             ddW_ddP, ddW_dtdP,ddW_ddt,
                             P, S, off, true);
  // -----------------------
  double eps = 1.0e-7;
  dfm2::CVec3d dP[3];
  dP[0].SetRandom(); dP[0] *= eps;
  dP[1].SetRandom(); dP[1] *= eps;
  dP[2].SetRandom(); dP[2] *= eps;
  const double dT[2] = {
    (2.0*rand()/(RAND_MAX+1.0)-1.0)*eps,
    (2.0*rand()/(RAND_MAX+1.0)-1.0)*eps };
  dfm2::CVec3d frm0[3], frm1[3];
  RodFrameTrans(frm0,
                S[0], P[1]-P[0], dP[1]-dP[0], dT[0]);
  RodFrameTrans(frm1,
                S[1], P[2]-P[1], dP[2]-dP[1], dT[1]);
  const dfm2::CVec3d p[3] = { P[0] + dP[0], P[1] + dP[1], P[2] + dP[2] };
  const dfm2::CVec3d s[2] = { frm0[0], frm1[0] };
  dfm2::CVec3d dw_dP[3];
  double dw_dt[2];
  double w = 0;
  {
    dfm2::CMat3d ddw_ddP[3][3];
    dfm2::CVec3d ddw_dtdP[2][3];
    double ddw_ddt[2][2];
    w = WdWddW_Rod(dw_dP, dw_dt,
                        ddw_ddP, ddw_dtdP, ddw_ddt,
                        p, s, off, true);
  }
  {
    const double val0 = (w-W)/eps;
    const double val1 = (+dW_dt[0]*dT[0]
                         +dW_dt[1]*dT[1]
                         +dW_dP[0]*dP[0]
                         +dW_dP[1]*dP[1]
                         +dW_dP[2]*dP[2])/eps;
    std::cout << "diff_W : " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const double val0 = (dw_dt[0]-dW_dt[0])/eps;
    const double val1 = (+ddW_ddt[ 0][0]*dT[0]
                         +ddW_ddt[ 0][1]*dT[1]
                         +ddW_dtdP[0][0]*dP[0]
                         +ddW_dtdP[0][1]*dP[1]
                         +ddW_dtdP[0][2]*dP[2])/eps;
    std::cout << "diff_T0: " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const double val0 = (dw_dt[1]-dW_dt[1])/eps;
    const double val1 = (+ddW_ddt[ 1][0]*dT[0]
                         +ddW_ddt[ 1][1]*dT[1]
                         +ddW_dtdP[1][0]*dP[0]
                         +ddW_dtdP[1][1]*dP[1]
                         +ddW_dtdP[1][2]*dP[2])/eps;
    std::cout << "diff_T1: " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[0]-dW_dP[0])/eps;
    const dfm2::CVec3d val1 = (+ddW_dtdP[0][0]*dT[0]
                               +ddW_dtdP[1][0]*dT[1]
                               +ddW_ddP[0][0]*dP[0]
                               +ddW_ddP[0][1]*dP[1]
                               +ddW_ddP[0][2]*dP[2])/eps;
    std::cout << "diff_V0: " << (val0-val1).Length() << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[1]-dW_dP[1])/eps;
    const dfm2::CVec3d val1 = (+ddW_dtdP[0][1]*dT[0]
                               +ddW_dtdP[1][1]*dT[1]
                               +ddW_ddP[1][0]*dP[0]
                               +ddW_ddP[1][1]*dP[1]
                               +ddW_ddP[1][2]*dP[2])/eps;
    std::cout << "diff_V2: " << (val0-val1).Length() << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[2]-dW_dP[2])/eps;
    const dfm2::CVec3d val1 = (+ddW_dtdP[0][2]*dT[0]
                               +ddW_dtdP[1][2]*dT[1]
                               +ddW_ddP[2][0]*dP[0]
                               +ddW_ddP[2][1]*dP[1]
                               +ddW_ddP[2][2]*dP[2])/eps;
    std::cout << "diff_V2: " << (val0-val1).Length() << " ---> " << val0 << " " << val1 << std::endl;
  }
  std::cout << std::endl;
}


double WdWddW_SquareLengthLineseg3D
(dfm2::CVec3d dW_dP[2],
 dfm2::CMat3d ddW_ddP[2][2],
 //
 const dfm2::CVec3d P[2],
 double L0)
{
  double l  = sqrt(+ (P[0][0]-P[1][0])*(P[0][0]-P[1][0])
                    + (P[0][1]-P[1][1])*(P[0][1]-P[1][1])
                    + (P[0][2]-P[1][2])*(P[0][2]-P[1][2]));
  dfm2::CVec3d v = P[0]-P[1];
  double R = L0-l;
  dW_dP[0] = (-R/l)*v;
  dW_dP[1] = (+R/l)*v;
  dfm2::CMat3d m = L0/(l*l*l)*dfm2::Mat3_OuterProduct(v, v) + (l-L0)/l*dfm2::Mat3_Identity(1.0);
  ddW_ddP[0][0] = m;
  ddW_ddP[0][1] = -m;
  ddW_ddP[1][0] = -m;
  ddW_ddP[1][1] = m;
  return 0.5*R*R;
}

void Check_WdWddW_SquareLengthLineseg3D()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(-1,+1);
  dfm2::CVec3d P[2];
  P[0].SetRandom();
  P[1].SetRandom();
  dfm2::CVec3d dW_dP[2];
  dfm2::CMat3d ddW_ddP[2][2];
  const double L0 = 1.0;
  double W = WdWddW_SquareLengthLineseg3D(dW_dP, ddW_ddP,
                                          P, L0);
  // -----
  double eps = 1.0e-5;
  dfm2::CVec3d dP[2];
  dP[0].SetRandom(); dP[0] *= eps;
  dP[1].SetRandom(); dP[1] *= eps;
  const dfm2::CVec3d p[2] = { P[0]+dP[0], P[1]+dP[1] };
  double w;
  dfm2::CVec3d dw_dP[2];
  {
    dfm2::CMat3d ddw_ddP[2][2];
    const double L0 = 1.0;
    w = WdWddW_SquareLengthLineseg3D(dw_dP, ddw_ddP,
                                     p, L0);
  }
  {
    const double val0 = (w-W)/eps;
    const double val1 = (+dW_dP[0]*dP[0]+dW_dP[1]*dP[1])/eps;
    std::cout << "diff_W : " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[0]-dW_dP[0])/eps;
    const dfm2::CVec3d val1 = (+ddW_ddP[0][0]*dP[0]+ddW_ddP[0][1]*dP[1])/eps;
    std::cout << "diff_P0 : " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
  {
    const dfm2::CVec3d val0 = (dw_dP[1]-dW_dP[1])/eps;
    const dfm2::CVec3d val1 = (+ddW_ddP[1][0]*dP[0]+ddW_ddP[1][1]*dP[1])/eps;
    std::cout << "diff_P0 : " << val0-val1 << " ---> " << val0 << " " << val1 << std::endl;
  }
}




void Check(){
  std::cout << "#######" << std::endl;
  std::cout << "##Frame" << std::endl;
  Check_WdWddW_RodFrameTrans();
  std::cout << "#" << std::endl;
  Check_WdWddW_RodFrameTrans();
  std::cout << "#" << std::endl;
  Check_WdWddW_RodFrameTrans();
  std::cout << "##########" << std::endl;
  std::cout << "##FrameDot" << std::endl;
  Check_WdWddW_DotFrame();
  std::cout << "#" << std::endl;
  Check_WdWddW_DotFrame();
  std::cout << "#" << std::endl;
  Check_WdWddW_DotFrame();
  std::cout << "#####" << std::endl;
  std::cout << "##Rod" << std::endl;
  Check_WdWddW_Rod();
  std::cout << "#" << std::endl;
  Check_WdWddW_Rod();
  std::cout << "#" << std::endl;
  Check_WdWddW_Rod();
  std::cout << "#" << std::endl;
  Check_WdWddW_Rod();
  std::cout << "#########" << std::endl;
  std::cout << "##LineSeg" << std::endl;
  Check_WdWddW_SquareLengthLineseg3D();
  std::cout << "#" << std::endl;
  Check_WdWddW_SquareLengthLineseg3D();
  std::cout << "#" << std::endl;
  Check_WdWddW_SquareLengthLineseg3D();
}

// -------------------------------------

void myGlutDisplay
(const std::vector<dfm2::CVec3d>& aP,
 const std::vector<dfm2::CVec3d>& aS,
 const std::vector<unsigned int>& aElemSeg)
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for(int ip=0;ip<aP.size();++ip){
    ::glVertex3d(aP[ip].x(), aP[ip].y(), aP[ip].z());
  }
  ::glEnd();
  // ------------
  ::glColor3d(0,0,0);
  ::glLineWidth(3);
  ::glBegin(GL_LINES);
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    unsigned int i0 = aElemSeg[iseg*2+0]; assert( i0 < aP.size() );
    unsigned int i1 = aElemSeg[iseg*2+1]; assert( i1 < aP.size() );
    ::glVertex3d(aP[i0].x(), aP[i0].y(), aP[i0].z());
    ::glVertex3d(aP[i1].x(), aP[i1].y(), aP[i1].z());
  }
  ::glEnd();
  // --------------
  ::glBegin(GL_LINES);
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    unsigned int i0 = aElemSeg[iseg*2+0]; assert( i0 < aP.size() );
    unsigned int i1 = aElemSeg[iseg*2+1]; assert( i1 < aP.size() );
    dfm2::CVec3d p01 = 0.5*(aP[i0]+aP[i1]);
    double l01 = (aP[i0]-aP[i1]).Length();
    dfm2::opengl::myGlVertex(p01);
    dfm2::opengl::myGlVertex(p01+(l01*0.5)*aS[iseg]);
  }
  ::glEnd();
}

void Solve
(std::vector<dfm2::CVec3d>& aP,
 std::vector<dfm2::CVec3d>& aS,
 dfm2::CMatrixSparse<double>& mats,
 const std::vector<dfm2::CVec3d>& aP0,
 const std::vector<unsigned int>& aElemSeg,
 const std::vector<unsigned int>& aElemRod,
 const std::vector<int>& aBCFlag)
{
  const unsigned int nNode = aBCFlag.size()/3;
  mats.SetZero();
  std::vector<double> vec_r;
  vec_r.assign(nNode*3, 0.0);
  std::vector<int> tmp_buffer;
  double W = 0;
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    const unsigned int i0 = aElemSeg[iseg*2+0];
    const unsigned int i1 = aElemSeg[iseg*2+1];
    const unsigned int* aINoel = aElemSeg.data()+iseg*2;
    const double L0 = (aP0[i0]-aP0[i1]).Length();
    const dfm2::CVec3d aPE[2] = { aP[i0], aP[i1] };
    // --------------
    dfm2::CVec3d dW_dP[2];
    dfm2::CMat3d ddW_ddP[2][2];
    W += WdWddW_SquareLengthLineseg3D(dW_dP, ddW_ddP,
                                      aPE, L0);
    {
      double eM[2*2*3*3];
      for(int in=0;in<2;++in){
        for(int jn=0;jn<2;++jn){
          ddW_ddP[in][jn].CopyValueToPtr(eM+(in*2+jn)*9);
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
  for(int irod=0;irod<aElemRod.size()/5;++irod){
    const unsigned int* aINoel = aElemRod.data()+irod*5;
    const unsigned int nP = aP.size();
    const dfm2::CVec3d aPE[3] = { aP[aINoel[0]], aP[aINoel[1]], aP[aINoel[2]] };
    const dfm2::CVec3d aSE[2] = {
      aS[aINoel[3]-nP],
      aS[aINoel[4]-nP] };
    const double off[3] = {0,0,0};
    // ------
    dfm2::CVec3d dW_dP[3];
    double dW_dt[2];
    dfm2::CMat3d ddW_ddP[3][3];
    dfm2::CVec3d ddW_dtdP[2][3];
    double ddW_ddt[2][2];
    W +=  WdWddW_Rod(dW_dP,dW_dt,ddW_ddP,ddW_dtdP,ddW_ddt,
                     aPE,aSE,off, false);
    {
      double eM[5][5][3][3];
      for(int i=0;i<5*5*3*3;++i){ (&eM[0][0][0][0])[i] = 0.0; }
      for(int in=0;in<3;++in){
        for(int jn=0;jn<3;++jn){
//          dfm2::CMat3d t = ddW_ddP[in][jn].Trans();
//          t.CopyValueToPtr(&eM[in][jn][0][0]);
          ddW_ddP[in][jn].CopyValueToPtr(&eM[in][jn][0][0]);
        }
      }
      for(int in=0;in<2;++in){
        for(int jn=0;jn<3;++jn){
          eM[3+in][jn][0][0] = ddW_dtdP[in][jn].x();
          eM[3+in][jn][0][1] = ddW_dtdP[in][jn].y();
          eM[3+in][jn][0][2] = ddW_dtdP[in][jn].z();
        }
      }
      for(int in=0;in<3;++in){
        for(int jn=0;jn<2;++jn){
          eM[in][jn+3][0][0] = ddW_dtdP[jn][in].x();
          eM[in][jn+3][1][0] = ddW_dtdP[jn][in].y();
          eM[in][jn+3][2][0] = ddW_dtdP[jn][in].z();
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
//  std::cout << dfm2::CheckSymmetry(mats) << std::endl;
//  mats.AddDia(0.1);
  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << dfm2::CheckSymmetry(mats) << std::endl;
  mats.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_r, aBCFlag,0);
  std::vector<double> vec_x;
  vec_x.assign(nNode*3, 0.0);
 {
   auto aConvHist = Solve_CG(vec_r.data(),vec_x.data(),
   vec_r.size(), 1.0e-4, 100, mats);
   if( aConvHist.size() > 0 ){
     std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
   }
 }
  /*
  {
    auto aConvHist = dfm2::Solve_BiCGStab(vec_r,vec_x,
                                          1.0e-4, 300, mats);
    if( aConvHist.size() > 0 ){
      std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
    }
  }
   */
  //    for(int i=0;i<vec_x.size();++i){
  //      std::cout << i << " " << vec_x[i] << std::endl;
  //    }
  assert( aS.size() == aElemSeg.size()/2 );
  for(int is=0;is<aS.size();++is){
    unsigned int i0 = aElemSeg[is*2+0];
    unsigned int i1 = aElemSeg[is*2+1];
    dfm2::CVec3d V01 = aP[i1]-aP[i0];
    dfm2::CVec3d du(vec_x[i1*3+0]-vec_x[i0*3+0],
                    vec_x[i1*3+1]-vec_x[i0*3+1],
                    vec_x[i1*3+2]-vec_x[i0*3+2]);
    const unsigned int np = aP.size();
    const double dtheta = vec_x[ np*3 + is*3 ];
    dfm2::CVec3d frm[3];
    RodFrameTrans(frm,
                  aS[is],V01,du,dtheta);
    aS[is] = frm[0];
  }
  for(int ip=0;ip<aP.size();++ip){
    aP[ip].p[0] += vec_x[ip*3+0];
    aP[ip].p[1] += vec_x[ip*3+1];
    aP[ip].p[2] += vec_x[ip*3+2];
  }
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    const unsigned int i0 = aElemSeg[iseg*2+0];
    const unsigned int i1 = aElemSeg[iseg*2+1];
    const dfm2::CVec3d& p0 = aP[i0];
    const dfm2::CVec3d& p1 = aP[i1];
    const dfm2::CVec3d e01 = (p1-p0).Normalize();
    aS[iseg] -= (aS[iseg]*e01)*e01;
    aS[iseg].SetNormalizedVector();
  }
}


void MakeProblemSetting
(std::vector<dfm2::CVec3d>& aP0,
 std::vector<unsigned int>& aElemSeg,
 std::vector<unsigned int>& aElemRod,
 std::vector<int>& aBCFlag, // if value this is not 0, it is fixed boundary condition
 unsigned int np)
{
  double elen = 2.0/np;
  aP0.resize(np);
  for(int ip=0;ip<np;++ip){
    aP0[ip] = dfm2::CVec3d(-1.0+ip*elen, 0.0, 0.0);
  };
  unsigned int ns = np-1;
  aElemSeg.resize(ns*2);
  for(int is=0;is<ns;++is){
    aElemSeg[is*2+0] = is+0;
    aElemSeg[is*2+1] = is+1;
  }
  unsigned int nr = ns-1;
  aElemRod.resize(nr*5);
  for(int ir=0;ir<nr;++ir){
    aElemRod[ir*5+0] = ir+0;
    aElemRod[ir*5+1] = ir+1;
    aElemRod[ir*5+2] = ir+2;
    aElemRod[ir*5+3] = np+ir+0;
    aElemRod[ir*5+4] = np+ir+1;
  };
  const unsigned int nNode = np+ns;
  aBCFlag.assign(nNode*3, 0);
  {
    unsigned int np = aP0.size();
    aBCFlag[0*3+0] = 1; aBCFlag[0*3+1] = 1; aBCFlag[0*3+2] = 1;
    aBCFlag[1*3+0] = 1; aBCFlag[1*3+1] = 1; aBCFlag[1*3+2] = 1;
    aBCFlag[(np+0)*3+0] = 1;
    for(int is=0;is<ns;++is){
      aBCFlag[(np+is)*3+1] = 1;
      aBCFlag[(np+is)*3+2] = 1;
    }
  }
}

int main(int argc,char* argv[])
{
  Check();
  // -----
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.5;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  // -----
  std::vector<dfm2::CVec3d> aP0;
  std::vector<unsigned int> aElemSeg, aElemRod;
  std::vector<int> aBCFlag;
  MakeProblemSetting(aP0,aElemSeg,aElemRod, aBCFlag,
                     30);
  dfm2::CMatrixSparse<double> mats;
  {
    unsigned int nNode = aElemSeg.size()/2 + aP0.size();
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aElemRod.data(), aElemRod.size()/5, 5, nNode);
    /*
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aElemSeg.data(), aElemRod.size()/2, 2, nNode);
     */
    dfm2::JArray_Sort(psup_ind, psup);
    mats.Initialize(nNode, 3, true);
    mats.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
  }
  // -----------------
  std::vector<dfm2::CVec3d> aP = aP0;
  std::vector<dfm2::CVec3d> aS(aElemSeg.size()/2);
  assert( aP.size() == aP0.size() );
  assert( aS.size() == aElemSeg.size()/2);
  // ------
  int iframe = 0;
  while (true)
  {
    if( iframe % 200 == 0 ){
      for(int ip=0;ip<aP.size();++ip){
        aP[ip] = aP0[ip];
        auto rnd = dfm2::CVec3d::Random()*3.0;
        if( aBCFlag[ip*3+0] == 0 ){ aP[ip].p[0] += rnd.x(); }
        if( aBCFlag[ip*3+1] == 0 ){ aP[ip].p[1] += rnd.y(); }
        if( aBCFlag[ip*3+2] == 0 ){ aP[ip].p[2] += rnd.z(); }
      }
      unsigned int ns = aS.size();
      for(int is=0;is<ns;++is){
        aS[is] = dfm2::CVec3d(0,1,0);
        auto rnd = dfm2::CVec3d::Random()*3.0;
        const unsigned int np = aP.size();
        if( aBCFlag[(np+is)*3+0] == 0 ){ aS[is].p[0] += rnd.x(); }
        if( aBCFlag[(np+is)*3+1] == 0 ){ aS[is].p[1] += rnd.y(); }
        if( aBCFlag[(np+is)*3+2] == 0 ){ aS[is].p[2] += rnd.z(); }
      }
      for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
        const unsigned int i0 = aElemSeg[iseg*2+0];
        const unsigned int i1 = aElemSeg[iseg*2+1];
        const dfm2::CVec3d& p0 = aP[i0];
        const dfm2::CVec3d& p1 = aP[i1];
        const dfm2::CVec3d e01 = (p1-p0).Normalize();
        assert( iseg < aS.size() );
        aS[iseg] -= (aS[iseg]*e01)*e01;
        aS[iseg].SetNormalizedVector();
      }
    }
    {
      Solve(aP, aS, mats,
            aP0, aElemSeg, aElemRod, aBCFlag);
    }
    iframe = iframe+1;
    // -------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aP,aS,aElemSeg);
    viewer.DrawEnd_oldGL();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
