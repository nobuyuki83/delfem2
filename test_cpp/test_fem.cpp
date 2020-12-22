/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h" // need to be defiend in the beginning
//
#include "delfem2/femrod.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/objf_geo3.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"
#include "delfem2/jagarray.h"
#include "delfem2/vec2.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/emat.h"
#include "delfem2/mats.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/primitive.h"
#include <random>


namespace dfm2 = delfem2;

// --------------------------------------

TEST(objfunc_v23, Check_CdC_TriStrain){
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_m5p5(-5, 5);
  // -----
  for(int itr=0;itr<200;++itr){
    const double P[3][2] = {
        {dist_m5p5(randomEng), dist_m5p5(randomEng) },
        {dist_m5p5(randomEng), dist_m5p5(randomEng) },
        {dist_m5p5(randomEng), dist_m5p5(randomEng) } };
    if( dfm2::Distance2(P[0],P[1]) < 0.1 ){ continue; }
    if( dfm2::Distance2(P[1],P[2]) < 0.1 ){ continue; }
    if( dfm2::Distance2(P[0],P[2]) < 0.1 ){ continue; }
    double a0 = dfm2::Area_Tri2(P[0], P[1], P[2]);
    if( fabs(a0) < 0.2 ) continue;
    const double p[3][3] = {
        {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng) },
        {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng) },
        {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng) } };
    const double eps = 1.0e-5;
    double C[3], dCdp[3][9];
    dfm2::PBD_CdC_TriStrain2D3D(C, dCdp, P, p);
    for(int ine=0;ine<3;++ine){
      for(int idim=0;idim<3;++idim){
        double p1[3][3]; for(int i=0;i<9;++i){ (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1[3], dCdp1[3][9];
        dfm2::PBD_CdC_TriStrain2D3D(C1, dCdp1, P, p1);
        EXPECT_NEAR( (C1[0]-C[0])/eps, dCdp[0][ine*3+idim], 1.0e-2);
        EXPECT_NEAR( (C1[1]-C[1])/eps, dCdp[1][ine*3+idim], 1.0e-2);
        EXPECT_NEAR( (C1[2]-C[2])/eps, dCdp[2][ine*3+idim], 1.0e-2);
      }
    }
  }
}

TEST(objfunc_v23, Bend)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0,1);
  const double eps = 1.0e-5;
  for(unsigned int itr=0;itr<10000;++itr) {
    const dfm2::CVec3d p[4] = {
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng)};
    if (dfm2::Distance(p[0], p[1]) < 0.1) { continue; }
    if (dfm2::Distance(p[0], p[2]) < 0.1) { continue; }
    if (dfm2::Distance(p[0], p[3]) < 0.1) { continue; }
    if (dfm2::Distance(p[1], p[2]) < 0.1) { continue; }
    if (dfm2::Distance(p[1], p[3]) < 0.1) { continue; }
    if (dfm2::Distance(p[2], p[3]) < 0.1) { continue; }
    if (dfm2::Area_Tri(p[0],p[2],p[3]) < 0.01){ continue; }
    if (dfm2::Area_Tri(p[1],p[2],p[3]) < 0.01){ continue; }
    double C;
    dfm2::CVec3d dC[4];
    GetConstConstDiff_Bend(C, dC, p[0], p[1], p[2], p[3]);
    for (int ino = 0; ino < 4; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        dfm2::CVec3d p1[4] = {p[0], p[1], p[2], p[3]};
        p1[ino][idim] += eps;
        double C1;
        dfm2::CVec3d dC1[4];
        GetConstConstDiff_Bend(C1, dC1, p1[0], p1[1], p1[2], p1[3]);
        const double val0 = (C1 - C) / eps;
        const double val1 = dC[ino][idim];
        EXPECT_NEAR(val0, val1, (1.0+fabs(val1)) * 1.5e-2);
      }
    }
  }
}


TEST(objfunc_v23, MIPS)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0, 1);
  const double eps = 1.0e-5;
  //
  for(unsigned int itr=0;itr<10000;++itr) {
    const double P[3][3] = {
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)}};
    if (dfm2::Distance3(P[0], P[1]) < 0.1) { continue; }
    if (dfm2::Distance3(P[0], P[2]) < 0.1) { continue; }
    if (dfm2::Distance3(P[1], P[2]) < 0.1) { continue; }
    if (dfm2::Area_Tri3(P[0],P[1],P[2])<0.01){ continue; }
    double p[3][3];
    {
      dfm2::CMat3d m;
      m.SetRotMatrix_Cartesian(0.3, 1.0, 0.5);
      for (int ino = 0; ino < 3; ++ino) {
        m.MatVec(P[ino], p[ino]);
      }
    }
    if (dfm2::Area_Tri3(p[0],p[1],p[2])<0.01){ continue; }
    double E, dE[3][3], ddE[3][3][3][3];
    dfm2::WdWddW_MIPS(E, dE, ddE,
        p, P);
    for (int ino = 0; ino < 3; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        double c1[3][3] = {
            {p[0][0], p[0][1], p[0][2]},
            {p[1][0], p[1][1], p[1][2]},
            {p[2][0], p[2][1], p[2][2]}};
        c1[ino][idim] += eps;
        double E1, dE1[3][3], ddE1[3][3][3][3];
        dfm2::WdWddW_MIPS(E1, dE1, ddE1,
                          c1, P);
        {
          const double val0 = (E1 - E) / eps;
          const double val1 = dE[ino][idim];
          EXPECT_NEAR(val0, val1, 1.0e-2*(1+fabs(val1)));
        }
        for (int jno = 0; jno < 3; ++jno) {
          for (int jdim = 0; jdim < 3; ++jdim) {
            const double val0 = (dE1[jno][jdim] - dE[jno][jdim]) / eps;
            const double val1 = ddE[jno][ino][jdim][idim];
            EXPECT_NEAR( val0, val1, 3.e-2*(1+fabs(val1)));
          }
        }
      }
    }
  }
}

TEST(objfunc_v23, distancetri2d3d)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0,1);
  //
  for(unsigned int itr=0;itr<10000;++itr) {
    const double P[3][2] = {  // undeformed triangle vertex positions
        { dist_01(randomEng), dist_01(randomEng)},
        { dist_01(randomEng), dist_01(randomEng)},
        { dist_01(randomEng), dist_01(randomEng)} };
    if (dfm2::Distance3(P[0], P[1]) < 0.1) { continue; }
    if (dfm2::Distance3(P[0], P[2]) < 0.1) { continue; }
    if (dfm2::Distance3(P[1], P[2]) < 0.1) { continue; }
    const double p[3][3] = { //  deformed triangle vertex positions)
        { dist_01(randomEng), dist_01(randomEng), dist_01(randomEng) },
        { dist_01(randomEng), dist_01(randomEng), dist_01(randomEng) },
        { dist_01(randomEng), dist_01(randomEng), dist_01(randomEng) } };
    // --------------
    double C[3], dCdp[3][9];
    dfm2::PBD_ConstraintProjection_DistanceTri2D3D(C, dCdp, P, p);
    for (int ine = 0; ine < 3; ++ine) {
      for (int idim = 0; idim < 3; ++idim) {
        double eps = 1.0e-6;
        double p1[3][3];
        for (int i = 0; i < 9; ++i) { (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1[3], dCdp1[3][9];
        dfm2::PBD_ConstraintProjection_DistanceTri2D3D(C1, dCdp1, P, p1);
        for (int jdim = 0; jdim < 3; ++jdim) {
          const double val0 = (C1[jdim] - C[jdim]) / eps;
          const double val1 = dCdp[jdim][ine * 3 + idim];
          EXPECT_NEAR(val0, val1, 1.0e-2*(fabs(val1)+1.0));
        }
      }
    }
  }
}


TEST(objfunc_v23, pbd_energy_stvk)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0,1);
  std::uniform_real_distribution<double> dist_12(1,2);

  //
  for(int itr=0;itr<1000;++itr){
    const double P[3][2] = {  // undeformed triangle vertex positions
        { dist_01(randomEng), dist_01(randomEng) },
        { dist_01(randomEng), dist_01(randomEng) },
        { dist_01(randomEng), dist_01(randomEng) } };
    const double A0 = dfm2::Area_Tri2(P[0], P[1], P[2]);
    if( A0 < 0.1 ) continue;
    const double p[3][3] = { //  deformed triangle vertex positions)
        { dist_01(randomEng),dist_01(randomEng), dist_01(randomEng) },
        { dist_01(randomEng),dist_01(randomEng), dist_01(randomEng) },
        { dist_01(randomEng),dist_01(randomEng), dist_01(randomEng) } };
    const double lambda = dist_12(randomEng);
    const double myu = dist_12(randomEng);
    // ---------------------
    double C, dCdp[9];
    dfm2::PBD_ConstraintProjection_EnergyStVK(C, dCdp, P, p, lambda, myu);
    for(int ine=0;ine<3;++ine){
      for(int idim=0;idim<3;++idim){
        double eps = 1.0e-8;
        double p1[3][3]; for(int i=0;i<9;++i){ (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1, dCdp1[9];
        dfm2::PBD_ConstraintProjection_EnergyStVK(C1, dCdp1, P, p1, lambda, myu);
        EXPECT_NEAR( (C1-C)/eps, dCdp[ine*3+idim], 1.0e-3*(1+fabs(dCdp[ine*3+idim])) );
      }
    }
  }
}

TEST(objfunc_v23, dWddW_RodFrameTrans)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0.0, 1.0);
  std::uniform_real_distribution<double> dist_m1p1(-1.0, +1.0);
  const double eps = 1.0e-6;
  //
  for(int itr=0;itr<10000;++itr){
    const dfm2::CVec3d V01 = dfm2::CVec3d::Random(dist_01, randomEng);
    if( V01.Length() < 0.1 ){ continue; }
    dfm2::CVec3d Frm[3];
    {
      Frm[2] = V01;
      Frm[2].SetNormalizedVector();
      Frm[0].SetRandom(dist_01, randomEng);
      Frm[0] -= (Frm[0]*Frm[2])*Frm[2];
      Frm[0].SetNormalizedVector();
      Frm[1] = (Frm[2]^Frm[0]);
    }
    dfm2::CVec3d Q; Q.SetRandom(dist_01, randomEng);
    //  Q = Frm[2];
    // --------------------------------
    double W[3] = { Q*Frm[0], Q*Frm[1], Q*Frm[2] };
    dfm2::CVec3d DW_Dv[3];
    double DW_Dt[3];
    {
      dfm2::CMat3d dF_dv[3];
      dfm2::CVec3d dF_dt[3];
      dfm2::DiffFrameRod(dF_dv, dF_dt,
                         V01.Length(), Frm);
      DW_Dv[0] = Q*dF_dv[0];
      DW_Dv[1] = Q*dF_dv[1];
      DW_Dv[2] = Q*dF_dv[2];
      DW_Dt[0] = Q*dF_dt[0];
      DW_Dt[1] = Q*dF_dt[1];
      DW_Dt[2] = Q*dF_dt[2];
    }
    // ---------
    const dfm2::CVec3d du = dfm2::CVec3d::Random(dist_01, randomEng)*eps;
    const double dtheta = dist_m1p1(randomEng)*eps;
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
      EXPECT_NEAR(val0, val1, 1.0e-3);
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
      EXPECT_NEAR(val0, val1, 1.0e-3);
    }
    for(int i=0;i<3;++i){
      dfm2::CMat3d DDW_DDv;
      dfm2::CVec3d DDW_DvDt;
      double DDW_DDt;
      DifDifFrameRod(DDW_DDv, DDW_DvDt, DDW_DDt,
          i,V01.Length(),Q,Frm);
      const dfm2::CVec3d vec0 = (dw_dv[i]-DW_Dv[i])/eps;
      const dfm2::CVec3d vec1 = (DDW_DvDt*dtheta+DDW_DDv*du)/eps;
      EXPECT_LT( (vec0-vec1).Length(), 3.5e-3);
    }
  }
}

TEST(objfunc_v23, WdWddW_DotFrame)
{
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0,1);
  std::uniform_real_distribution<double> dist_m1p1(-1,+1);
  const double eps = 1.0e-5;
  //
  for(int itr=0;itr<10000;++itr){
    dfm2::CVec3d P[3] = {
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng) };
    if( (P[1]-P[0]).Length() < 0.01 ){ continue; }
    if( (P[2]-P[1]).Length() < 0.01 ){ continue; }
    dfm2::CVec3d S[2];
    {
      S[0].SetRandom(dist_01,randomEng);
      const dfm2::CVec3d U0 = (P[1]-P[0]).Normalize();
      S[0] -= (S[0]*U0)*U0;
      S[0].SetNormalizedVector();
    }
    {
      S[1].SetRandom(dist_01,randomEng);
      const dfm2::CVec3d U1 = (P[2]-P[1]).Normalize();
      S[1] -= (S[1]*U1)*U1;
      S[1].SetNormalizedVector();
    }
    const double off[3] = {
        dist_m1p1(randomEng),
        dist_m1p1(randomEng),
        dist_m1p1(randomEng) };
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
    const dfm2::CVec3d dP[3] = {
        dfm2::CVec3d::Random(dist_01,randomEng)*eps,
        dfm2::CVec3d::Random(dist_01,randomEng)*eps,
        dfm2::CVec3d::Random(dist_01,randomEng)*eps };
    const double dT[2] = {
        dist_m1p1(randomEng)*eps,
        dist_m1p1(randomEng)*eps };
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
      EXPECT_NEAR(val0, val1, 3.5e-2*(1.0+fabs(val1)));
    }
    {
      const double val0 = (dw_dt[0]-dW_dt[0])/eps;
      const double val1 = (+ddW_ddt[ 0][0]*dT[0]
                           +ddW_ddt[ 0][1]*dT[1]
                           +ddW_dtdP[0][0]*dP[0]
                           +ddW_dtdP[0][1]*dP[1]
                           +ddW_dtdP[0][2]*dP[2])/eps;
      EXPECT_NEAR(val0, val1, 3.0e-2*(1.0+fabs(val1)));
    }
    {
      const double val0 = (dw_dt[1]-dW_dt[1])/eps;
      const double val1 = (+ddW_ddt[ 1][0]*dT[0]
                           +ddW_ddt[ 1][1]*dT[1]
                           +ddW_dtdP[1][0]*dP[0]
                           +ddW_dtdP[1][1]*dP[1]
                           +ddW_dtdP[1][2]*dP[2])/eps;
      EXPECT_NEAR(val0, val1, 3.0e-2*(1.0+fabs(val1)));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[0]-dW_dP[0])/eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][0]*dT[0]
                                 +ddW_dtdP[1][0]*dT[1]
                                 +ddW_ddP[0][0]*dP[0]
                                 +ddW_ddP[0][1]*dP[1]
                                 +ddW_ddP[0][2]*dP[2])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-2*(1+val1.Length()) );
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[1]-dW_dP[1])/eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][1]*dT[0]
                                 +ddW_dtdP[1][1]*dT[1]
                                 +ddW_ddP[1][0]*dP[0]
                                 +ddW_ddP[1][1]*dP[1]
                                 +ddW_ddP[1][2]*dP[2])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-2*(1+val1.Length()) );
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[2]-dW_dP[2])/eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][2]*dT[0]
                                 +ddW_dtdP[1][2]*dT[1]
                                 +ddW_ddP[2][0]*dP[0]
                                 +ddW_ddP[2][1]*dP[1]
                                 +ddW_ddP[2][2]*dP[2])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-2*(1+val1.Length()) );
    }
  }
}


TEST(objfunc_v23, WdWddW_Rod)
{
  std::random_device rd;
  std::mt19937 randomEng(rd());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  std::uniform_real_distribution<double> dist_12(+1, +2);
  std::uniform_real_distribution<double> dist_01(+0, +1);
  double eps = 1.0e-7;
  //
  for(int itr=0;itr<100;++itr){
    const double stiff_bendtwist[3] = {
        dist_12(randomEng),
        dist_12(randomEng),
        dist_12(randomEng) };
    const dfm2::CVec3d P[3] = {
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng) };
    if( dfm2::Distance(P[0],P[1]) < 0.1 ){ continue; }
    if( dfm2::Distance(P[1],P[2]) < 0.1 ){ continue; }

    dfm2::CVec3d S[2];
    {
      S[0].SetRandom(dist_01,randomEng);
      const dfm2::CVec3d U0 = (P[1]-P[0]).Normalize();
      S[0] -= (S[0]*U0)*U0;
      S[0].SetNormalizedVector();
    }
    {
      S[1].SetRandom(dist_01,randomEng);
      const dfm2::CVec3d U1 = (P[2]-P[1]).Normalize();
      S[1] -= (S[1]*U1)*U1;
      S[1].SetNormalizedVector();
    }
    const double off[3] = {
        dist_m1p1(randomEng),
        dist_m1p1(randomEng),
        dist_m1p1(randomEng) };
    // ------------------------
    dfm2::CVec3d dW_dP[3];
    double dW_dt[2];
    dfm2::CMat3d ddW_ddP[3][3];
    dfm2::CVec3d ddW_dtdP[2][3];
    double ddW_ddt[2][2];
    double W = WdWddW_Rod(dW_dP,dW_dt,
                          ddW_ddP, ddW_dtdP,ddW_ddt,
                          stiff_bendtwist, P, S, off, true);
    // -----------------------
    const dfm2::CVec3d dP[3] = {
        dfm2::CVec3d::Random(dist_01,randomEng)*eps,
        dfm2::CVec3d::Random(dist_01,randomEng)*eps,
        dfm2::CVec3d::Random(dist_01,randomEng)*eps };
    const double dT[2] = {
        dist_m1p1(randomEng)*eps,
        dist_m1p1(randomEng)*eps };
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
                     stiff_bendtwist, p, s, off, true);
    }
    {
      const double val0 = (w-W)/eps;
      const double val1 = (+dW_dt[0]*dT[0]
                           +dW_dt[1]*dT[1]
                           +dW_dP[0]*dP[0]
                           +dW_dP[1]*dP[1]
                           +dW_dP[2]*dP[2])/eps;
      EXPECT_NEAR(val0, val1, 5.0e-3*(1+fabs(val1)) );
    }
    {
      const double val0 = (dw_dt[0]-dW_dt[0])/eps;
      const double val1 = (+ddW_ddt[ 0][0]*dT[0]
                           +ddW_ddt[ 0][1]*dT[1]
                           +ddW_dtdP[0][0]*dP[0]
                           +ddW_dtdP[0][1]*dP[1]
                           +ddW_dtdP[0][2]*dP[2])/eps;
      EXPECT_NEAR(val0, val1, 2.5e-3*(1+fabs(val1)) );
    }
    {
      const double val0 = (dw_dt[1]-dW_dt[1])/eps;
      const double val1 = (+ddW_ddt[ 1][0]*dT[0]
                           +ddW_ddt[ 1][1]*dT[1]
                           +ddW_dtdP[1][0]*dP[0]
                           +ddW_dtdP[1][1]*dP[1]
                           +ddW_dtdP[1][2]*dP[2])/eps;
      EXPECT_NEAR(val0, val1, 1.0e-3*(1+fabs(val1)) );
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[0]-dW_dP[0])/eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][0]*dT[0]
                                 +ddW_dtdP[1][0]*dT[1]
                                 +ddW_ddP[0][0]*dP[0]
                                 +ddW_ddP[0][1]*dP[1]
                                 +ddW_ddP[0][2]*dP[2])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-3*(1+val1.Length()) );
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[1]-dW_dP[1])/eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][1]*dT[0]
                                 +ddW_dtdP[1][1]*dT[1]
                                 +ddW_ddP[1][0]*dP[0]
                                 +ddW_ddP[1][1]*dP[1]
                                 +ddW_ddP[1][2]*dP[2])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-3*(1+val1.Length()) );
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[2]-dW_dP[2])/eps;
      const dfm2::CVec3d val1 = (+ddW_dtdP[0][2]*dT[0]
                                 +ddW_dtdP[1][2]*dT[1]
                                 +ddW_ddP[2][0]*dP[0]
                                 +ddW_ddP[2][1]*dP[1]
                                 +ddW_ddP[2][2]*dP[2])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-3*(1+val1.Length()) );
    }
  }
}

TEST(objfunc_v23, WdWddW_SquareLengthLineseg3D)
{
  std::random_device rd;
  std::mt19937 rndeng(rd());
  std::uniform_real_distribution<double> dist_m1p1(-1,+1);
  std::uniform_real_distribution<double> dist_12(+1,+2);
  std::uniform_real_distribution<double> dist_01(+0,+1);
  double eps = 1.0e-5;
  //
  for(int itr=0;itr<10000;++itr){
    const double stiff_stretch = dist_12(rndeng);
    const dfm2::CVec3d P[2] = {
        dfm2::CVec3d::Random(dist_01,rndeng),
        dfm2::CVec3d::Random(dist_01,rndeng) };
    if( (P[0]-P[1]).Length() < 0.1 ){ continue; }
    dfm2::CVec3d dW_dP[2];
    dfm2::CMat3d ddW_ddP[2][2];
    const double L0 = 1.0;
    double W = WdWddW_SquareLengthLineseg3D(dW_dP, ddW_ddP,
                                            stiff_stretch, P, L0);
    // -----
    const dfm2::CVec3d dP[2] = {
        dfm2::CVec3d::Random(dist_01,rndeng)*eps,
        dfm2::CVec3d::Random(dist_01,rndeng)*eps };
    const dfm2::CVec3d p[2] = { P[0]+dP[0], P[1]+dP[1] };
    double w;
    dfm2::CVec3d dw_dP[2];
    {
      dfm2::CMat3d ddw_ddP[2][2];
      w = WdWddW_SquareLengthLineseg3D(dw_dP, ddw_ddP,
          stiff_stretch, p, L0);
    }
    {
      const double val0 = (w-W)/eps;
      const double val1 = (+dW_dP[0]*dP[0]+dW_dP[1]*dP[1])/eps;
      EXPECT_NEAR( val0, val1, 1.0e-3*(1+fabs(val1)) );
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[0]-dW_dP[0])/eps;
      const dfm2::CVec3d val1 = (+ddW_ddP[0][0]*dP[0]+ddW_ddP[0][1]*dP[1])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-3*(1+val1.Length()) );
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[1]-dW_dP[1])/eps;
      const dfm2::CVec3d val1 = (+ddW_ddP[1][0]*dP[0]+ddW_ddP[1][1]*dP[1])/eps;
      EXPECT_LT( (val0-val1).Length(), 1.0e-3*(1+val1.Length()) );
    }
  }
}

TEST(objfunc_v23, arap)
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3D_Cube(aXYZ0, aTri, 10);
  const size_t np = aXYZ0.size()/3;
  
  std::vector<unsigned int> psup_ind, psup;
  {
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aTri.data(), aTri.size()/3, 3,
                               aXYZ0.size()/3);
    dfm2::JArray_Sort(psup_ind, psup);
  }
  
  // base mesh
  // ------------------
  // deformed mesh
  
  std::vector<double> aXYZ1 = aXYZ0;
  dfm2::Translate_Points3(aXYZ1, 0.1, 0.2, 0.3);
  dfm2::Rotate_Points3(aXYZ1, 0.1, 0.2, 0.3);
  
  std::vector<double> aQuat1;
  { // initialize rotation
    aQuat1.resize(np*4);
    for(unsigned int ip=0;ip<np;++ip){
      dfm2::Quat_Identity(aQuat1.data()+4*ip);
    }
    for(int itr=0;itr<40;++itr){
      dfm2::UpdateRotationsByMatchingCluster_Linear(aQuat1,
                                                    aXYZ0,aXYZ1,psup_ind,psup);
    }
  }
  
  // ------------------

  std::random_device rd;
  std::mt19937 rndeng(rd());
  std::uniform_real_distribution<double> dist_m1p1(-1,+1);
  std::vector<double> dXYZ12(np*3);
  for(unsigned int i=0;i<np*3;++i){
    dXYZ12[i] = dist_m1p1(rndeng)*0.01;
  }
  assert( aXYZ1.size() == np*3 );
  assert( aQuat1.size() == np*4 );
  double w1 = dfm2::W_ArapEnergy(aXYZ0, aXYZ1, aQuat1, psup_ind, psup);
  
  const double eps = 1.0e-5;
  std::vector<double> aXYZ2 = aXYZ1;
  for(int i=0;i<aXYZ2.size();++i){ aXYZ2[i] += eps*dXYZ12[i]; }
  
  std::vector<double> aQuat2 = aQuat1;
  dfm2::UpdateRotationsByMatchingCluster_Linear(
      aQuat2,
      aXYZ0, aXYZ2, psup_ind, psup);
  
  double w2 = dfm2::W_ArapEnergy(
      aXYZ0, aXYZ2, aQuat2, psup_ind, psup);
  std::vector<double> aRes2;
  dfm2::dW_ArapEnergy(
      aRes2,
      aXYZ0, aXYZ2, aQuat2, psup_ind, psup);
  
  // ---------------------------------
  
  std::vector<double> aRes1;
  dfm2::dW_ArapEnergy(aRes1,
                      aXYZ0, aXYZ1, aQuat1, psup_ind, psup);
  
  {
    double dw = dfm2::DotX(aRes1.data(),dXYZ12.data(),aRes1.size());
    double val2 = (w2-w1)/eps;
    EXPECT_NEAR(dw, val2, 1.0e-5*(fabs(dw)+1.0) );
  }
  
  // ---------------------------------

  dfm2::CMatrixSparse<double> Mat;
  {
    Mat.Initialize(np, 3, true);
    std::vector<unsigned int> psup_ind1, psup1;
    dfm2::JArray_Extend(psup_ind1, psup1,
                        psup_ind.data(), psup_ind.size(), psup.data());
    dfm2::JArray_Sort(psup_ind1, psup1);
    assert( psup_ind1.size() == np+1 );
    Mat.SetPattern(psup_ind1.data(), psup_ind1.size(), psup1.data(), psup1.size());
    Mat.SetZero();
    std::vector<int> tmp_buffer;
    for(unsigned int ip=0;ip<np;++ip){
      std::vector<unsigned int> aIP;
      for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        const unsigned int jp = psup[ipsup];
        aIP.push_back(jp);
      }
      aIP.push_back(ip);
      std::vector<double> eM;
      dfm2::ddW_ArapEnergy(eM,
                           aIP,aXYZ0,aQuat1);
      Mat.Mearge(aIP.size(), aIP.data(), aIP.size(), aIP.data(), 9, eM.data(), tmp_buffer);
    }
  }
  EXPECT_LT( dfm2::CheckSymmetry(Mat), 1.0e-10);
  
  {
    std::vector<double> aRes12(np*3);
    Mat.MatVec(aRes12.data(),
               1.0, dXYZ12.data(), 0.0);
    
    for(int i=0;i<np*3;++i){
      double val0 = (aRes2[i]-aRes1[i])/eps;
      double val1 = aRes12[i];
//      std::cout << i << " " << val0 << " " << val1 << std::endl;
      EXPECT_NEAR( val0, val1, 1.0e-5*(1+fabs(val1)) );
    }
  }
    
}

// ------------------------------------------------------------

TEST(fem,plate_bending_mitc3_emat)
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<> dist0(-0.5, +0.5);
  std::uniform_real_distribution<> dist1(+1.0e-10, +1.0);
  for(int itr=0;itr<200;++itr){
    double C[3][2];
    for(int i=0;i<6;++i){
      (&C[0][0])[i] = 10.0*dist0(mt);
    }
    double a0 = dfm2::Area_Tri2(C[0], C[1], C[2]);
    if( a0 < 0.1 ) continue;
    double u0[3][3];
    for(int i=0;i<9;++i){
      (&u0[0][0])[i] = 1.0*dist0(mt);
    }
    double thickness1 = dist1(mt);
    double lambda1 = dist1(mt);
    double myu1 = dist1(mt);
    const double eps = 1.0e-5;
    // -------------------
    double W0, dW0[3][3], ddW0[3][3][3][3];
    W0 = 0.0;
    for(int i=0;i<9;++i){ (&dW0[0][0])[i] = 0.0; }
    for(int i=0;i<81;++i){ (&ddW0[0][0][0][0])[i] = 0.0; }
    dfm2::WdWddW_PlateBendingMITC3(W0,dW0,ddW0,
                                   C,u0,
                                   thickness1,lambda1,myu1);
    for(int ino=0;ino<3;++ino){
      for(int idof=0;idof<3;++idof){
        double u1[3][3]; for(int i=0;i<9;++i){ (&u1[0][0])[i] = (&u0[0][0])[i]; }
        u1[ino][idof] += eps;
        double W1, dW1[3][3], ddW1[3][3][3][3];
        W1 = 0.0;
        for(int i=0;i<9;++i){ (&dW1[0][0])[i] = 0.0; }
        for(int i=0;i<81;++i){ (&ddW1[0][0][0][0])[i] = 0.0; }
        dfm2::WdWddW_PlateBendingMITC3(W1,dW1,ddW1,
                                       C,u1,
                                       thickness1,lambda1,myu1);
        EXPECT_NEAR( (W1-W0)/eps,
                    dW0[ino][idof],
                    3.0e-3*(1.0+fabs(dW0[ino][idof])) );
        for(int jno=0;jno<3;++jno){
          for(int jdof=0;jdof<3;++jdof){
            EXPECT_NEAR( (dW1[jno][jdof]-dW0[jno][jdof])/eps,
                        ddW0[ino][jno][idof][jdof],
                        1.0e-3*(1.0+fabs(ddW0[ino][jno][idof][jdof])) );
          }
        }
      }
    }
  }
}


TEST(fem,plate_bending_mitc3_cantilever)
{
  std::random_device rd;
  std::mt19937 rndeng(rd());
  std::uniform_real_distribution<double> dist_12(1,+2);
  //
  const double lambda = 0.0;
  const double lenx0 = 1.0;
  const double leny0 = 0.2;
  const double thickness0 = 0.05;
  const double myu0 = 10000.0;
  const double rho0 = 1.0;
  const double gravity_z0 = -10.0;
  const double elen0 = 0.03;
  for(int itr = 0;itr<10;++itr ){
    const double lenx = lenx0*dist_12(rndeng);
    const double leny = leny0*dist_12(rndeng);
    const double thickness = thickness0*dist_12(rndeng);
    const double myu = myu0*dist_12(rndeng);
    const double rho = rho0*dist_12(rndeng);
    const double gravity_z = gravity_z0*dist_12(rndeng);
    const double elen = elen0*dist_12(rndeng);
    std::vector<unsigned int> aTri;
    std::vector<double> aXY0;
    {
      std::vector< std::vector<double> > aaXY;
      {
        aaXY.resize(1);
        aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(-leny*0.5);
        aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(-leny*0.5);
        aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(+leny*0.5);
        aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(+leny*0.5);
      }
      // ---------------------
      std::vector<dfm2::CDynPntSur> aPo2D;
      std::vector<dfm2::CDynTri> aETri;
      std::vector<dfm2::CVec2d> aVec2;
      GenMesh(aPo2D, aETri, aVec2,
              aaXY, elen, elen);
      MeshTri2D_Export(aXY0,aTri,
                       aVec2,aETri);
    }
    std::vector<int> aBCFlag; // boundary condition flag
    {
      const int np = (int)aXY0.size()/2;
      aBCFlag.assign(np*3, 0);
      for(int ip=0;ip<np;++ip){
        const double px = aXY0[ip*2+0];
        //    const double py = aXY0[ip*2+1];
        if( fabs(px-(-lenx*0.5)) < 0.0001 ){
          aBCFlag[ip*3+0] = 1;
          aBCFlag[ip*3+1] = 1;
          aBCFlag[ip*3+2] = 1;
        }
      }
    }
    dfm2::CMatrixSparse<double> mat_A;
    dfm2::CPreconditionerILU<double> ilu_A;
    {
      std::vector<unsigned int> psup_ind, psup;
      dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                                 aTri.data(), aTri.size()/3, 3,
                                 (int)aXY0.size()/2);
      dfm2::JArray_Sort(psup_ind, psup);
      //
      const int np = (int)aXY0.size()/2;
      mat_A.Initialize(np, 3, true);
      mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
      ilu_A.Initialize_ILU0(mat_A);
    }
    std::vector<double> aVal;
    aVal.assign(aXY0.size()/2*3, 0.0);
    std::vector<double> vec_b;
    {
      const int np = (int)aXY0.size()/2;
      const int nDoF = np*3;
      // -------------------
      mat_A.SetZero();
      vec_b.assign(nDoF, 0.0);
      dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(
          mat_A,vec_b.data(),
          thickness,lambda,myu,
          rho,gravity_z,
          aXY0.data(), aXY0.size()/2,
          aTri.data(), aTri.size()/3,
          aVal.data());
      mat_A.SetFixedBC(aBCFlag.data());
      dfm2::setRHS_Zero(vec_b, aBCFlag,0);
      // --------------------------
      std::vector<double> vec_x;
      {
        ilu_A.SetValueILU(mat_A);
        ilu_A.DoILUDecomp();
        vec_x.resize(vec_b.size());
        std::vector<double> conv = Solve_PCG(
            vec_b.data(), vec_x.data(),
            vec_b.size(),
            1.0e-5, 1000,
            mat_A, ilu_A);
//        std::cout << "convergence   nitr:" << conv.size() << "    res:" << conv[conv.size()-1] << std::endl;
        EXPECT_LT( conv.size(), 1000 );
        EXPECT_LT( conv[conv.size()-1], 1.0e-5);
      }
      // -------------------------
      dfm2::XPlusAY(aVal,nDoF,aBCFlag,
                    1.0,vec_x);
    }
    {
      assert( fabs(lambda)<1.0e-10 );
      const double E = myu*2.0;
      const double I = thickness*thickness*thickness*leny/12.0;
      const double W = thickness*lenx*leny*rho*gravity_z;
      const double w = W/lenx;
      const double disp = w*(lenx*lenx*lenx*lenx)/(8.0*E*I);
//      std::cout << "disp:" << disp << std::endl;
      for(int ip=0;ip<aXY0.size()/2;++ip){
        const double px = aXY0[ip*2+0];
        if( fabs(px-(+lenx*0.5)) > 0.0001 ){ continue; }
        EXPECT_LE( fabs(aVal[ip*3+0] - disp), 0.002*fabs(disp) );
//        std::cout << aVal[ip*3+0] << " " << disp << "  " << fabs(aVal[ip*3+0] - disp)/disp << std::endl;
      }
    }
  }
}
