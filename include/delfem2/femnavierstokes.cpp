/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femnavierstokes.h"

#include <cmath>
#include <cassert>

DFM2_INLINE void delfem2::MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1(
    double rho,
    double myu,
    double g_x,
    double g_y,
    double dt,
    double gamma,
    const double coords[][2],
    const double velo[][2],
    const double press[],
    const double acc[][2],
    const double apress[],
    double emat_uu[][3][2][2],
    double emat_up[][3][2],
    double emat_pu[][3][2],
    double emat_pp[][3],
    double eres_u[3][2],
    double eres_p[3])
{
  const int nno = 3;
  const int ndim = 2;
  
  const double area = femutil::TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[3][2], const_term[3];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  const double velo_ave[2] = {
    (velo[0][0]+velo[1][0]+velo[2][0])/3.0,
    (velo[0][1]+velo[1][1]+velo[2][1])/3.0 };
  
  double tau;
  {  // Calc Stabilization Parameter
    
    const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
    const double h = sqrt( area / 3.14 )*2;
    const double tau_c = h*0.5/norm_v;
    const double cou_c = norm_v*dt/h;
    if( norm_v*h*rho*1.0e-30 > myu ){ // Re = \infty
      const double dtmp1 = 1/(cou_c*cou_c)+1;
      tau = tau_c / sqrt(dtmp1);
    }
    else if( norm_v*h*rho < myu*1.0e-30 ){ // Re = 0
      tau = h*h*rho*0.5/myu;
    }
    else{
      const double re_c = 0.5*norm_v*h*rho/myu;  // 0.5*norm_v*h*rho/myu;
      const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
      tau = tau_c / sqrt(dtmp1);
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  double eCmat_uu[3][3][2][2];
  // viscousity
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_uu[ino][jno][0][0] = area*myu*dldx[ino][0]*dldx[jno][0];
      eCmat_uu[ino][jno][0][1] = area*myu*dldx[ino][1]*dldx[jno][0];
      eCmat_uu[ino][jno][1][0] = area*myu*dldx[ino][0]*dldx[jno][1];
      eCmat_uu[ino][jno][1][1] = area*myu*dldx[ino][1]*dldx[jno][1];
      const double dtmp1 = area*myu*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
      eCmat_uu[ino][jno][0][0] += dtmp1;
      eCmat_uu[ino][jno][1][1] += dtmp1;
    }
  }
  {  // advection
    const double dtmp0[2] = { velo[0][0]+velo[1][0]+velo[2][0], velo[0][1]+velo[1][1]+velo[2][1] };
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = (dldx[jno][0]*dtmp0[0]+dldx[jno][1]*dtmp0[1]);
      for(int ino=0;ino<nno;ino++){
        double dtmp2 = dtmp1 + (dldx[jno][0]*velo[ino][0]+dldx[jno][1]*velo[ino][1]);
        dtmp2 *= area*rho*0.083333333333333;
        eCmat_uu[ino][jno][0][0] += dtmp2;
        eCmat_uu[ino][jno][1][1] += dtmp2;
      }
    }
  }
  {  // SUPG for advection term
    double tmp_mat[ndim][ndim] = { {0,0}, {0,0} };
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        tmp_mat[0][0] += velo[ino][0]*velo[jno][0];
        tmp_mat[0][1] += velo[ino][0]*velo[jno][1];
        tmp_mat[1][0] += velo[ino][1]*velo[jno][0];
        tmp_mat[1][1] += velo[ino][1]*velo[jno][1];
      }
      tmp_mat[0][0] += velo[ino][0]*velo[ino][0];
      tmp_mat[1][1] += velo[ino][1]*velo[ino][1];
    }
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        double dtmp1 = 0.0;
        dtmp1 += dldx[ino][0]*dldx[jno][0]*tmp_mat[0][0]
        +dldx[ino][0]*dldx[jno][1]*tmp_mat[0][1]
        +dldx[ino][1]*dldx[jno][0]*tmp_mat[1][0]
        +dldx[ino][1]*dldx[jno][1]*tmp_mat[1][1];
        dtmp1 *= tau*rho*area*0.083333333333333333333;
        eCmat_uu[ino][jno][0][0] += dtmp1;
        eCmat_uu[ino][jno][1][1] += dtmp1;
      }
    }
  }
  
  double eCmat_up[3][3][2], eCmat_pu[3][3][2];
  {   // add pressure gradient, non compression term
    const double dtmp1 = area*0.33333333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eCmat_up[ino][jno][0] = -dtmp1*dldx[ino][0];
        eCmat_up[ino][jno][1] = -dtmp1*dldx[ino][1];
        eCmat_pu[ino][jno][0] =  dtmp1*dldx[jno][0];
        eCmat_pu[ino][jno][1] =  dtmp1*dldx[jno][1];
      }
    }
  }
  // SUPG for pressure gradient
  for(int ino=0;ino<nno;ino++){
    double dtmp1 = (dldx[ino][0]*velo_ave[0]+dldx[ino][1]*velo_ave[1])*tau*area;
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] += dtmp1*dldx[jno][0];
      eCmat_up[ino][jno][1] += dtmp1*dldx[jno][1];
    }
  }
  // PSPG for advection term
  for(int jno=0;jno<nno;jno++){
    const double dtmp1 = (dldx[jno][0]*velo_ave[0]+dldx[jno][1]*velo_ave[1])*tau*area;
    for(int ino=0;ino<nno;ino++){
      eCmat_pu[ino][jno][0] += dtmp1*dldx[ino][0];
      eCmat_pu[ino][jno][1] += dtmp1*dldx[ino][1];
    }
  }
  
  double eCmat_pp[3][3];
  // PSPG for pressure term
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pp[ino][jno] = area*tau/rho*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
    }
  }
  
  // ------------------------------
  double eMmat_uu[3][3][2][2];
  {  // add inertia term
    const double dtmp1 = area*rho*0.0833333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat_uu[ino][jno][0][0] = dtmp1;
        eMmat_uu[ino][jno][0][1] = 0.0;
        eMmat_uu[ino][jno][1][0] = 0.0;
        eMmat_uu[ino][jno][1][1] = dtmp1;
      }
      eMmat_uu[ino][ino][0][0] += dtmp1;
      eMmat_uu[ino][ino][1][1] += dtmp1;
    }
  }
  // supg for inertia term
  for(int jno=0;jno<nno;jno++){
    double tmp_vec[ndim] = { 0.0, 0.0 };
    for(unsigned int kno=0;kno<nno;kno++){
      tmp_vec[0] += velo[kno][0];
      tmp_vec[1] += velo[kno][1];
    }
    tmp_vec[0] += velo[jno][0];
    tmp_vec[1] += velo[jno][1];
    for(int ino=0;ino<nno;ino++){
      const double dtmp1 = (dldx[ino][0]*tmp_vec[0]+dldx[ino][1]*tmp_vec[1])*rho*tau*area*0.083333333333333;
      eMmat_uu[ino][jno][0][0] += dtmp1;
      eMmat_uu[ino][jno][1][1] += dtmp1;
    }
  }
  
  // PSPG for inertia term
  double eMmat_pu[3][3][2];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eMmat_pu[ino][jno][0] = -tau*area*0.33333333333333333*dldx[ino][0];
      eMmat_pu[ino][jno][1] = -tau*area*0.33333333333333333*dldx[ino][1];
    }
  }
  
  // eternal force term
  for(int ino=0;ino<nno;ino++){
    eres_u[ino][0] = area*rho*g_x*0.3333333333333333333;
    eres_u[ino][1] = area*rho*g_y*0.3333333333333333333;
  }
  for(int ino=0;ino<nno;ino++){
    eres_p[ino] = 0.0;
  }
  
  /*
   // SUPG external force
   for(unsigned int ino=0;ino<nno;ino++){
   const double dtmp1 = area*tau*rho*(ave_velo[0]*dldx[ino][0]+ave_velo[1]*dldx[ino][1]);
   eres_u[ino][0] += dtmp1*g_x;
   eres_u[ino][1] += dtmp1*g_y;
   }
   */
  
  /*
   // PSPG external force
   for(unsigned int ino=0;ino<nno;ino++){
   eres_p[ino] += area*tau*(dldx[ino][0]*g_x+dldx[ino][1]*g_y);
   }
   */
  //////////////////////////////////////////////////////////////////////
  
  {
    double dtmp1 = gamma*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
    }
    for(int i=0;i<nno*nno;i++){
      (&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
      +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
      +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
      eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
      +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
      +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
      +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
      +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
    }
  }
}

DFM2_INLINE void delfem2::EMat_NavierStokes2D_Dynamic_P1(
    double myu, double rho, double g_x, double g_y,
    const double dt_timestep, const double gamma_newmark,
    const double coords[][2],
    const double velo_press[3][3], const double acc_apress[3][3],
    double emat[3][3][3][3],
    double eres[3][3])
{
  const int nno = 3;
  
  const double velo[3][2] = {
    {velo_press[0][0],velo_press[0][1]},
    {velo_press[1][0],velo_press[1][1]},
    {velo_press[2][0],velo_press[2][1]} };
  const double press[3] = { velo_press[0][2], velo_press[1][2], velo_press[2][2] };
  const double acc[3][2] = {
    {acc_apress[0][0],acc_apress[0][1]},
    {acc_apress[1][0],acc_apress[1][1]},
    {acc_apress[2][0],acc_apress[2][1]} };
  const double apress[3] = { acc_apress[0][2], acc_apress[1][2], acc_apress[2][2] };
  //
  double emat_uu[3][3][2][2], emat_up[3][3][2], emat_pu[3][3][2], emat_pp[3][3];
  double eres_u[3][2], eres_p[3];
  MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1(rho, myu, g_x, g_y,
                                              dt_timestep, gamma_newmark,
                                              coords,
                                              velo, press, acc, apress,
                                              emat_uu, emat_up, emat_pu, emat_pp,
                                              eres_u, eres_p);
  //
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
      emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
      emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
      emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
      emat[ino][jno][0][2] = emat_up[ino][jno][0];
      emat[ino][jno][1][2] = emat_up[ino][jno][1];
      emat[ino][jno][2][0] = emat_pu[ino][jno][0];
      emat[ino][jno][2][1] = emat_pu[ino][jno][1];
      emat[ino][jno][2][2] = emat_pp[ino][jno];
    }
  }
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = eres_u[ino][0];
    eres[ino][1] = eres_u[ino][1];
    eres[ino][2] = eres_p[ino];
  }
}


DFM2_INLINE void delfem2::MakeMat_NavierStokes3D_Dynamic_Newmark_P1P1(
  double rho, double myu, double g_x, double g_y, double g_z,
  double dt, double gamma,
  const double coords[4][3],
  const double velo[4][3], const double press[4], const double acc[4][3], const double apress[4],
  double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
   double eres_u[4][3], double eres_p[4])
{
  const int nno = 4;
  const int ndim = 3;
  //
  const double vol = delfem2::femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);

  const double velo_ave[3] = {
    (velo[0][0]+velo[1][0]+velo[2][0]+velo[3][0])*0.25,
    (velo[0][1]+velo[1][1]+velo[2][1]+velo[3][1])*0.25,
    (velo[0][2]+velo[1][2]+velo[2][2]+velo[3][2])*0.25 };

  ////
  double tau;
  {  // Calc Stabilization Parameter

    const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]+velo_ave[2]*velo_ave[2]);
    double tmp = 0.0;
    for(int ino=0;ino<nno;ino++){
      double tmp1 = dldx[ino][0]*velo_ave[0]+dldx[ino][1]*velo_ave[1]+dldx[ino][2]*velo_ave[2];
      tmp += fabs(tmp1);
    }
    double h_ugn = 0;
    if( fabs(tmp) > 1.0e-10 ){
      h_ugn = 2.0 * norm_v / tmp;
    }
    double tau_sugn1 = 0.0;
    if( fabs(norm_v) > 1.0e-10 ){
      tau_sugn1 = h_ugn * 0.5 / norm_v;
    }
    double tau_sugn2 = dt * 0.50;
    double tau_sugn3 = h_ugn * h_ugn * 0.25 / myu;
    double inv1 = 0.0;
    if( fabs(tau_sugn1) > 1.0e-10 ){
      inv1 = 1.0 / tau_sugn1;
    }
    double inv2 = 0.0;
    if( fabs(tau_sugn2) > 1.0e-10 ){
      inv2 = 1.0 / tau_sugn2;
    }
    double inv3 = 0.0;
    if( fabs(tau_sugn3) > 1.0e-10 ){
      inv3 = 1.0 / tau_sugn3;
    }
    double tau_supg = sqrt( inv1*inv1 + inv2*inv2 + inv3*inv3 );
    tau_supg = 1.0 / tau_supg;
//    double tau_pspg = tau_supg;
    tau = tau_supg*0.25;
//    std::cout << tau << "   " << tau_sugn1 << " " << tau_sugn2 << " " << tau_sugn3 << " " << inv1 << " " << inv2 << " " << inv3 << std::endl;
    /*
  ////
    const double h = pow(vol/3.14, 0.333333333333)*2;
    const double tau_c = h*0.5/norm_v;
    const double cou_c = norm_v*dt/h;
    if( norm_v*h*rho*1.0e-30 > myu ){ // Re = \infty
      const double dtmp1 = 1/(cou_c*cou_c)+1;
      tau = tau_c / sqrt(dtmp1);
    }
    else if( norm_v*h*rho < myu*1.0e-30 ){ // Re = 0
      tau = h*h*rho*0.5/myu;
    }
    else{
      const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
      const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
      tau = tau_c / sqrt(dtmp1);
    }
     */
  }

  /*
  double ave_velo[3];
  { // average velocity
    ave_velo[0] = 0.0;
    ave_velo[1] = 0.0;
    ave_velo[2] = 0.0;
    for(int ino=0;ino<nno;ino++){
      ave_velo[0] += velo[ino][0];
      ave_velo[1] += velo[ino][1];
      ave_velo[2] += velo[ino][2];
    }
    ave_velo[0] *= 0.25;
    ave_velo[1] *= 0.25;
    ave_velo[2] *= 0.25;
  }
  */
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  double eCmat_uu[4][4][3][3];
  // viscousity
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_uu[ino][jno][0][0] = vol*myu*dldx[ino][0]*dldx[jno][0];
      eCmat_uu[ino][jno][0][1] = vol*myu*dldx[ino][1]*dldx[jno][0];
      eCmat_uu[ino][jno][0][2] = vol*myu*dldx[ino][2]*dldx[jno][0];
      ///
      eCmat_uu[ino][jno][1][0] = vol*myu*dldx[ino][0]*dldx[jno][1];
      eCmat_uu[ino][jno][1][1] = vol*myu*dldx[ino][1]*dldx[jno][1];
      eCmat_uu[ino][jno][1][2] = vol*myu*dldx[ino][2]*dldx[jno][1];
      ///
      eCmat_uu[ino][jno][2][0] = vol*myu*dldx[ino][0]*dldx[jno][2];
      eCmat_uu[ino][jno][2][1] = vol*myu*dldx[ino][1]*dldx[jno][2];
      eCmat_uu[ino][jno][2][2] = vol*myu*dldx[ino][2]*dldx[jno][2];
      ////
      const double dtmp1 = vol*myu*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      eCmat_uu[ino][jno][0][0] += dtmp1;
      eCmat_uu[ino][jno][1][1] += dtmp1;
      eCmat_uu[ino][jno][2][2] += dtmp1;
    }
  }
  {	// advection
    const double dtmp0[3] = {
      velo[0][0]+velo[1][0]+velo[2][0]+velo[3][0],
      velo[0][1]+velo[1][1]+velo[2][1]+velo[3][1],
      velo[0][2]+velo[1][2]+velo[2][2]+velo[3][2] };
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = (dldx[jno][0]*dtmp0[0]+dldx[jno][1]*dtmp0[1]+dldx[jno][2]*dtmp0[2]);
      for(int ino=0;ino<nno;ino++){
        double dtmp2 = dtmp1 + (dldx[jno][0]*velo[ino][0]+dldx[jno][1]*velo[ino][1]+dldx[jno][2]*velo[ino][2]);
        dtmp2 *= vol*rho*0.05;
        eCmat_uu[ino][jno][0][0] += dtmp2;
        eCmat_uu[ino][jno][1][1] += dtmp2;
        eCmat_uu[ino][jno][2][2] += dtmp2;
      }
    }
  }
  {	// SUPG for advection term
    double tmp_mat[ndim][ndim] = { {0,0,0}, {0,0,0}, {0,0,0} };
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        tmp_mat[0][0] += velo[ino][0]*velo[jno][0];
        tmp_mat[0][1] += velo[ino][0]*velo[jno][1];
        tmp_mat[0][2] += velo[ino][0]*velo[jno][2];
        ////
        tmp_mat[1][0] += velo[ino][1]*velo[jno][0];
        tmp_mat[1][1] += velo[ino][1]*velo[jno][1];
        tmp_mat[1][2] += velo[ino][1]*velo[jno][2];
        ////
        tmp_mat[2][0] += velo[ino][2]*velo[jno][0];
        tmp_mat[2][1] += velo[ino][2]*velo[jno][1];
        tmp_mat[2][2] += velo[ino][2]*velo[jno][2];
      }
      tmp_mat[0][0] += velo[ino][0]*velo[ino][0];
      tmp_mat[1][1] += velo[ino][1]*velo[ino][1];
      tmp_mat[2][2] += velo[ino][2]*velo[ino][2];
    }
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        double dtmp1 = 0.0;
        dtmp1 +=
         dldx[ino][0]*dldx[jno][0]*tmp_mat[0][0]
        +dldx[ino][0]*dldx[jno][1]*tmp_mat[0][1]
        +dldx[ino][0]*dldx[jno][2]*tmp_mat[0][2]
        +dldx[ino][1]*dldx[jno][0]*tmp_mat[1][0]
        +dldx[ino][1]*dldx[jno][1]*tmp_mat[1][1]
        +dldx[ino][1]*dldx[jno][2]*tmp_mat[1][2]
        +dldx[ino][2]*dldx[jno][0]*tmp_mat[2][0]
        +dldx[ino][2]*dldx[jno][1]*tmp_mat[2][1]
        +dldx[ino][2]*dldx[jno][2]*tmp_mat[2][2];
        dtmp1 *= tau*rho*vol*0.05;
        eCmat_uu[ino][jno][0][0] += dtmp1;
        eCmat_uu[ino][jno][1][1] += dtmp1;
        eCmat_uu[ino][jno][2][2] += dtmp1;
      }
    }
  }
  
  double eCmat_up[4][4][3], eCmat_pu[4][4][3];
  {   // add pressure gradient, non compression term
    const double dtmp1 = vol*0.25;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eCmat_up[ino][jno][0] = -dtmp1*dldx[ino][0];
        eCmat_up[ino][jno][1] = -dtmp1*dldx[ino][1];
        eCmat_up[ino][jno][2] = -dtmp1*dldx[ino][2];
        eCmat_pu[ino][jno][0] = +dtmp1*dldx[jno][0];
        eCmat_pu[ino][jno][1] = +dtmp1*dldx[jno][1];
        eCmat_pu[ino][jno][2] = +dtmp1*dldx[jno][2];
      }
    }
  }
  // SUPG for pressure gradient
  for(int ino=0;ino<nno;ino++){
    const double dtmp1 = (dldx[ino][0]*velo_ave[0]+dldx[ino][1]*velo_ave[1]+dldx[ino][2]*velo_ave[2])*tau*vol;
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] += dtmp1*dldx[jno][0];
      eCmat_up[ino][jno][1] += dtmp1*dldx[jno][1];
      eCmat_up[ino][jno][2] += dtmp1*dldx[jno][2];
    }
  }
  // PSPG for advection term
  for(int jno=0;jno<nno;jno++){
    const double dtmp1 = (dldx[jno][0]*velo_ave[0]+dldx[jno][1]*velo_ave[1]+dldx[jno][2]*velo_ave[2])*tau*vol;
    for(int ino=0;ino<nno;ino++){
      eCmat_pu[ino][jno][0] += dtmp1*dldx[ino][0];
      eCmat_pu[ino][jno][1] += dtmp1*dldx[ino][1];
      eCmat_pu[ino][jno][2] += dtmp1*dldx[ino][2];
    }
  }
  
  double eCmat_pp[4][4];
  // PSPG for pressure term
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pp[ino][jno] = vol*tau/rho*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  ////////////////
  double eMmat_uu[4][4][3][3];
  {	// add inertia term
    const double dtmp1 = vol*rho*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat_uu[ino][jno][0][0] = dtmp1;
        eMmat_uu[ino][jno][0][1] = 0.0;
        eMmat_uu[ino][jno][0][2] = 0.0;
        ////
        eMmat_uu[ino][jno][1][0] = 0.0;
        eMmat_uu[ino][jno][1][1] = dtmp1;
        eMmat_uu[ino][jno][1][2] = 0.0;
        ////
        eMmat_uu[ino][jno][2][0] = 0.0;
        eMmat_uu[ino][jno][2][1] = 0.0;
        eMmat_uu[ino][jno][2][2] = dtmp1;
      }
      eMmat_uu[ino][ino][0][0] += dtmp1;
      eMmat_uu[ino][ino][1][1] += dtmp1;
      eMmat_uu[ino][ino][2][2] += dtmp1;
    }
  }
  // supg for inertia term
  for(int jno=0;jno<nno;jno++){
    double tmp_vec[ndim] = { 0.0, 0.0, 0.0 };
    for(unsigned int kno=0;kno<nno;kno++){
      tmp_vec[0] += velo[kno][0];
      tmp_vec[1] += velo[kno][1];
      tmp_vec[2] += velo[kno][2];
    }
    tmp_vec[0] += velo[jno][0];
    tmp_vec[1] += velo[jno][1];
    tmp_vec[2] += velo[jno][2];
    for(int ino=0;ino<nno;ino++){
      const double dtmp1 = (dldx[ino][0]*tmp_vec[0]+dldx[ino][1]*tmp_vec[1]+dldx[ino][2]*tmp_vec[2])*rho*tau*vol*0.05;
      eMmat_uu[ino][jno][0][0] += dtmp1;
      eMmat_uu[ino][jno][1][1] += dtmp1;
      eMmat_uu[ino][jno][2][2] += dtmp1;
    }
  }
  
  // PSPG for inertia term
  double eMmat_pu[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eMmat_pu[ino][jno][0] = -tau*vol*0.25*dldx[ino][0];
      eMmat_pu[ino][jno][1] = -tau*vol*0.25*dldx[ino][1];
      eMmat_pu[ino][jno][2] = -tau*vol*0.25*dldx[ino][2];
    }
  }
  
  // eternal force term
  for(int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*rho*g_x*0.25;
    eres_u[ino][1] = vol*rho*g_y*0.25;
    eres_u[ino][2] = vol*rho*g_z*0.25;
  }
  for(int ino=0;ino<nno;ino++){
    eres_p[ino] = 0.0;
  }
  
   // SUPG external force
//   for(unsigned int ino=0;ino<nno;ino++){
//			const double dtmp1 = area*tau*rho*(ave_velo[0]*dldx[ino][0]+ave_velo[1]*dldx[ino][1]);
//			eres_u[ino][0] += dtmp1*g_x;
//			eres_u[ino][1] += dtmp1*g_y;
//   }
  
  
   // PSPG external force
//   for(unsigned int ino=0;ino<nno;ino++){
//			eres_p[ino] += area*tau*(dldx[ino][0]*g_x+dldx[ino][1]*g_y);
//   }
  
  // ----------------------------------
  {
    double dtmp1 = gamma*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
    }
    for(int i=0;i<nno*nno;i++){
      (&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_u[ino][0]
      -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_uu[ino][jno][0][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_uu[ino][jno][0][0]*acc[jno][0]
      +  eMmat_uu[ino][jno][0][1]*acc[jno][1]
      +  eMmat_uu[ino][jno][0][2]*acc[jno][2]
      +  eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
      eres_u[ino][1]
      -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_uu[ino][jno][1][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_uu[ino][jno][1][0]*acc[jno][0]
      +  eMmat_uu[ino][jno][1][1]*acc[jno][1]
      +  eMmat_uu[ino][jno][1][2]*acc[jno][2]
      +  eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
      eres_u[ino][2]
      -= eCmat_uu[ino][jno][2][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_uu[ino][jno][2][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_uu[ino][jno][2][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_uu[ino][jno][2][0]*acc[jno][0]
      +  eMmat_uu[ino][jno][2][1]*acc[jno][1]
      +  eMmat_uu[ino][jno][2][2]*acc[jno][2]
      +  eCmat_up[ino][jno][2]*(press[jno]+dt*apress[jno]);
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_p[ino]
      -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_pu[ino][jno][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_pu[ino][jno][0]*acc[jno][0]
      +  eMmat_pu[ino][jno][1]*acc[jno][1]
      +  eMmat_pu[ino][jno][2]*acc[jno][2]
      +  eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
    }
  }
}


DFM2_INLINE void delfem2::MakeMat_NavierStokes3D_Dynamic_P1
(double myu, double rho, double g_x, double g_y, double g_z,
 const double dt_timestep, const double gamma_newmark,
 const double coords[4][3],
 const double velo_press[4][4], const double acc_apress[4][4],
 double emat[4][4][4][4],
 double eres[4][4])
{
  const int nno = 4;
  const double velo[4][3] = {
    {velo_press[0][0],velo_press[0][1],velo_press[0][2]},
    {velo_press[1][0],velo_press[1][1],velo_press[1][2]},
    {velo_press[2][0],velo_press[2][1],velo_press[2][2]},
    {velo_press[3][0],velo_press[3][1],velo_press[3][2]} };
  const double press[4] = { velo_press[0][3], velo_press[1][3], velo_press[2][3], velo_press[3][3] };
  const double acc[4][3] = {
    {acc_apress[0][0],acc_apress[0][1],acc_apress[0][2]},
    {acc_apress[1][0],acc_apress[1][1],acc_apress[1][2]},
    {acc_apress[2][0],acc_apress[2][1],acc_apress[2][2]},
    {acc_apress[3][0],acc_apress[3][1],acc_apress[3][2]} };
  const double apress[4] = { acc_apress[0][3], acc_apress[1][3], acc_apress[2][3], acc_apress[3][3] };
  ////
  double emat_uu[4][4][3][3], emat_up[4][4][3], emat_pu[4][4][3], emat_pp[4][4];
  double eres_u[4][3], eres_p[4];
  MakeMat_NavierStokes3D_Dynamic_Newmark_P1P1(rho, myu, g_x, g_y,g_z,
                                              dt_timestep, gamma_newmark,
                                              coords,
                                              velo, press, acc, apress,
                                              emat_uu, emat_up, emat_pu, emat_pp,
                                              eres_u, eres_p);
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
      emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
      emat[ino][jno][0][2] = emat_uu[ino][jno][0][2];
      //
      emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
      emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
      emat[ino][jno][1][2] = emat_uu[ino][jno][1][2];
      //
      emat[ino][jno][2][0] = emat_uu[ino][jno][2][0];
      emat[ino][jno][2][1] = emat_uu[ino][jno][2][1];
      emat[ino][jno][2][2] = emat_uu[ino][jno][2][2];
      //
      emat[ino][jno][0][3] = emat_up[ino][jno][0];
      emat[ino][jno][1][3] = emat_up[ino][jno][1];
      emat[ino][jno][2][3] = emat_up[ino][jno][2];
      //
      emat[ino][jno][3][0] = emat_pu[ino][jno][0];
      emat[ino][jno][3][1] = emat_pu[ino][jno][1];
      emat[ino][jno][3][2] = emat_pu[ino][jno][2];
      //
      emat[ino][jno][3][3] = emat_pp[ino][jno];
    }
  }
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = eres_u[ino][0];
    eres[ino][1] = eres_u[ino][1];
    eres[ino][2] = eres_u[ino][2];
    eres[ino][3] = eres_p[ino];
  }
}