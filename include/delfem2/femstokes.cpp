/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femstokes.h"

#include <cmath>

// ------------------------
// below tri

DFM2_INLINE void delfem2::MakeMat_Stokes2D_Static_P1P1(
    double alpha,
    double g_x,
    double g_y,
    const double coords[][2],
    const double velo[3][2],
    const double press[3],
    double emat_uu[][3][2][2],
    double emat_up[][3][2],
    double emat_pu[][3][2],
    double emat_pp[][3],
    double eres_u[][2],
    double eres_p[3])
{
  const unsigned int nno = 3;
  const unsigned int ndim = 2;
  
  const double area = femutil::TriArea2D(coords[0],coords[1],coords[2]);
  
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  
  //-------------------------
  
  for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat_uu[0][0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      const double dtmp1 = area*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
      emat_uu[ino][jno][0][0] = dtmp1;
      emat_uu[ino][jno][1][1] = dtmp1;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_up[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_up[ino][jno][0] += area*dldx[ino][0]*0.333333333333333333333333;
      emat_up[ino][jno][1] += area*dldx[ino][1]*0.333333333333333333333333;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_pu[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pu[ino][jno][0] += area*dldx[jno][0]*0.333333333333333333333333;
      emat_pu[ino][jno][1] += area*dldx[jno][1]*0.333333333333333333333333;
    }
  }
  for(unsigned int i=0;i<nno*nno;i++){ *(&emat_pp[0][0]+i) = 0.0; }
  double tau; // relaxation parameter
  {
    const double h = sqrt( area / 3.14 )*2;
    tau = -h*h/alpha*0.1;
//    tau = 0.0;
  }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pp[ino][jno] = area*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
    }
  }
  
  for(unsigned int ino=0;ino<nno;ino++){
    eres_u[ino][0] = area*g_x*0.3333333333333333333;
    eres_u[ino][1] = area*g_y*0.3333333333333333333;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_uu[ino][jno][0][0]*velo[jno][0]+emat_uu[ino][jno][0][1]*velo[jno][1];
      eres_u[ino][1] -= emat_uu[ino][jno][1][0]*velo[jno][0]+emat_uu[ino][jno][1][1]*velo[jno][1];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_up[ino][jno][0]*press[jno];
      eres_u[ino][1] -= emat_up[ino][jno][1]*press[jno];
    }
  }
  for(unsigned int ino=0;ino<nno;ino++){
    eres_p[ino] = 0.0;
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pu[ino][jno][0]*velo[jno][0]+emat_pu[ino][jno][1]*velo[jno][1];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pp[ino][jno]*press[jno];
    }
  }
}

DFM2_INLINE void delfem2::EMat_Stokes2D_Static_P1(
    double alpha, double g_x, double g_y,
    const double coords[][2],
    const double velo_press[3][3],
    double emat[3][3][3][3],
    double eres[3][3])
{
  const int nno = 3;
  
  const double velo[3][2] = {
    {velo_press[0][0],velo_press[0][1]},
    {velo_press[1][0],velo_press[1][1]},
    {velo_press[2][0],velo_press[2][1]} };
  const double press[3] = { velo_press[0][2], velo_press[1][2], velo_press[2][2] };
  ////
  double emat_uu[3][3][2][2], emat_up[3][3][2], emat_pu[3][3][2], emat_pp[3][3];
  double eres_u[3][2], eres_p[3];
  MakeMat_Stokes2D_Static_P1P1(
      alpha, g_x, g_y,
      coords, velo, press,
      emat_uu, emat_up, emat_pu, emat_pp,
      eres_u, eres_p);
  ////
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


DFM2_INLINE void delfem2::EMat_Stokes2D_Dynamic_P1(
    double alpha,
    double rho,
    double g_x,
    double g_y,
    const double dt_timestep,
    const double gamma_newmark,
    const double coords[][2],
    const double velo_press[3][3],
    const double acc_apress[3][3],
    double emat[3][3][3][3],
    double eres[3][3])
{
  const int nno = 3;
  //
  const double area = femutil::TriArea2D(coords[0],coords[1],coords[2]);
  //
  double dldx[nno][2], const_term[nno];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  //
  double tau;
  {
    const double h = sqrt( area / 3.14 )*2;
    tau = -h*h/alpha*0.1;
  }
  
  double eCmat[3][3][3][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = area*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
      eCmat[ino][jno][0][0] = dtmp1;
      eCmat[ino][jno][1][1] = dtmp1;
      eCmat[ino][jno][0][1] = 0.0;
      eCmat[ino][jno][1][0] = 0.0;
      eCmat[ino][jno][2][0] = area*dldx[jno][0]*0.333333333333333333333;
      eCmat[ino][jno][2][1] = area*dldx[jno][1]*0.333333333333333333333;
      eCmat[ino][jno][0][2] = area*dldx[ino][0]*0.333333333333333333333;
      eCmat[ino][jno][1][2] = area*dldx[ino][1]*0.333333333333333333333;
      eCmat[ino][jno][2][2] = area*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
    }
  }
  
  //
  double eMmat[3][3][3][3];
  EmatConsistentMassTri2<3>(
      eMmat,
      area*rho, false);
  
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = area*g_x*0.33333333333333333333;
    eres[ino][1] = area*g_y*0.33333333333333333333;
    eres[ino][2] = 0.0;
  }
  
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<3*3*3*3;++i){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eCmat[0][0][0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<3;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<3;jdim++){
          eres[ino][idim]
          -= eCmat[ino][jno][idim][jdim]*(velo_press[jno][jdim]+dt_timestep*acc_apress[jno][jdim])
          + eMmat[ino][jno][idim][jdim]*acc_apress[jno][jdim];
        }
      }
    }
  }
}

// above: tri
// ----------------------------------------------
// below: tet

DFM2_INLINE void delfem2::MakeMat_Stokes3D_Static_P1P1(
  double alpha,
  double g_x,
  double g_y,
  double g_z,
  const double coords[4][3],
  const double velo[4][3],
  const double press[4],
  double emat_uu[4][4][3][3],
  double emat_up[4][4][3],
  double emat_pu[4][4][3],
  double emat_pp[4][4],
  double eres_u[4][3],
  double eres_p[4])
{
  const unsigned int nno = 4;
  const unsigned int ndim = 3;
  //
  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  //
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);
  //
  SetEMatLaplaceTet(
      emat_uu, vol*alpha, dldx);
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_up[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_up[ino][jno][0] += vol*dldx[ino][0]*0.25;
      emat_up[ino][jno][1] += vol*dldx[ino][1]*0.25;
      emat_up[ino][jno][2] += vol*dldx[ino][2]*0.25;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_pu[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pu[ino][jno][0] += vol*dldx[jno][0]*0.25;
      emat_pu[ino][jno][1] += vol*dldx[jno][1]*0.25;
      emat_pu[ino][jno][2] += vol*dldx[jno][2]*0.25;
    }
  }
  {
    const double h = pow(vol/3.14, 0.3333333333)*2;
    const double tau = -h*h/alpha*0.1;
    SetEMatLaplaceTet(
        emat_pp,
        vol*tau, dldx);
  }

  // ---------
  for(unsigned int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*g_x*0.25;
    eres_u[ino][1] = vol*g_y*0.25;
    eres_u[ino][2] = vol*g_z*0.25;
  }
  AddEmatEvecScale3<4>(eres_u,emat_uu,velo,-1.);
  for(unsigned int ino=0;ino<nno;ino++) {
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_up[ino][jno][0]*press[jno];
      eres_u[ino][1] -= emat_up[ino][jno][1]*press[jno];
      eres_u[ino][2] -= emat_up[ino][jno][2]*press[jno];
    }
  }
  // ---------
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  eres_p[3] = 0;
  for(unsigned int ino=0;ino<nno;ino++) {
    eres_p[ino] = 0.0;
    for (unsigned int jno = 0; jno < nno; jno++) {
      eres_p[ino] -= emat_pu[ino][jno][0] * velo[jno][0] + emat_pu[ino][jno][1] * velo[jno][1] +
                     emat_pu[ino][jno][2] * velo[jno][2];
    }
  }
  for(unsigned int ino=0;ino<nno;ino++) {
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pp[ino][jno]*press[jno];
    }
  }
}



DFM2_INLINE void delfem2::MakeMat_Stokes3D_Static_P1(
    double alpha,
    double g_x,
    double g_y,
    double g_z,
    const double coords[4][3],
    const double velo_press[4][4],
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
  //
  double emat_uu[4][4][3][3], emat_up[4][4][3], emat_pu[4][4][3], emat_pp[4][4];
  double eres_u[4][3], eres_p[4];
  MakeMat_Stokes3D_Static_P1P1(
      alpha, g_x, g_y, g_z,
      coords, velo, press,
      emat_uu, emat_up, emat_pu, emat_pp,
      eres_u, eres_p);
  //
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
      ////
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



DFM2_INLINE void delfem2::MakeMat_Stokes3D_Dynamic_Newmark_P1P1(
    double alpha,
    double rho,
    double g_x,
    double g_y,
    double g_z,
    const double dt_timestep,
    const double gamma_newmark,
    const double coords[4][3],
    const double velo[4][3],
    const double press[4],
    const double acc[4][3],
    const double apress[4],
    double emat_uu[4][4][3][3],
    double emat_up[4][4][3],
    double emat_pu[4][4][3],
    double emat_pp[4][4],
    double eres_u[4][3],
    double eres_p[4])
{
  constexpr int nno = 4;
  constexpr int ndim = 3;
  //
  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0],coords[1],coords[2],coords[3]);
  
  // ------------------------------
  double eCmat_uu[4][4][3][3];
  SetEMatLaplaceTet(eCmat_uu, vol*alpha, dldx);
  
  double eCmat_up[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] = vol*dldx[ino][0]*0.25;
      eCmat_up[ino][jno][1] = vol*dldx[ino][1]*0.25;
      eCmat_up[ino][jno][2] = vol*dldx[ino][2]*0.25;
    }
  }
  
  double eCmat_pu[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pu[ino][jno][0] = vol*dldx[jno][0]*0.25;
      eCmat_pu[ino][jno][1] = vol*dldx[jno][1]*0.25;
      eCmat_pu[ino][jno][2] = vol*dldx[jno][2]*0.25;
    }
  }
  
  double tau;
  {
    const double h = pow( vol / 3.14, 0.3333333 )*2;
    tau = -h*h/alpha*0.1;
  }
  
  double eCmat_pp[4][4];
  SetEMatLaplaceTet(eCmat_pp,vol*tau,dldx);

  // -------------------------
  double eMmat_uu[4][4][3][3];
  SetEmatConsistentMassTet(
      eMmat_uu,vol*rho);
  
  for(int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*g_x*0.25;
    eres_u[ino][1] = vol*g_y*0.25;
    eres_u[ino][2] = vol*g_z*0.25;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  eres_p[3] = 0;
  
  // --------------------------------
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          for(int jdim=0;jdim<ndim;jdim++){
            emat_uu[ino][jno][idim][jdim] = eMmat_uu[ino][jno][idim][jdim]+dtmp1*eCmat_uu[ino][jno][idim][jdim];
          }
        }
      }
    }
    
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          emat_up[ino][jno][idim] = dtmp1*eCmat_up[ino][jno][idim];
          emat_pu[ino][jno][idim] = dtmp1*eCmat_pu[ino][jno][idim];
        }
      }
    }
    
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        emat_pp[ino][jno] = dtmp1*eCmat_pp[ino][jno];
      }
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<ndim;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<ndim;jdim++){
          eres_u[ino][idim] -= eCmat_uu[ino][jno][idim][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim])
          + eMmat_uu[ino][jno][idim][jdim]*acc[jno][jdim];
        }
      }
      for(int jno=0;jno<nno;jno++){
        eres_u[ino][idim] -= eCmat_up[ino][jno][idim]*(press[jno]+dt_timestep*apress[jno]);
      }
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      for(int jdim=0;jdim<ndim;jdim++){
        eres_p[ino] -= eCmat_pu[ino][jno][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim]);
      }
    }
    for(int jno=0;jno<nno;jno++){
      eres_p[ino] -= eCmat_pp[ino][jno]*(press[jno]+dt_timestep*apress[jno]);
    }
  }
}


DFM2_INLINE void delfem2::MakeMat_Stokes3D_Dynamic_P1(
    double alpha,
    double rho,
    double g_x,
    double g_y,
    double g_z,
    const double dt_timestep,
    const double gamma_newmark,
    const double coords[4][3],
    const double velo_press[4][4],
    const double acc_apress[4][4],
    double emat[4][4][4][4],
    double eres[4][4])
{
  constexpr int nno = 4;
  constexpr int ndim = 4;
  //
  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  //
  double dldx[nno][3], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);
  //
  double tau;
  {
    const double h = pow( vol / 3.14, 0.333333333 )*2;
    tau = -h*h/alpha*0.1;
  }
  // -------------------------
  double eCmat[4][4][4][4];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = vol*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      eCmat[ino][jno][0][0] = dtmp1;
      eCmat[ino][jno][1][1] = dtmp1;
      eCmat[ino][jno][2][2] = dtmp1;
      eCmat[ino][jno][0][1] = 0.0;
      eCmat[ino][jno][0][2] = 0.0;
      eCmat[ino][jno][1][0] = 0.0;
      eCmat[ino][jno][1][2] = 0.0;
      eCmat[ino][jno][2][0] = 0.0;
      eCmat[ino][jno][2][1] = 0.0;
      eCmat[ino][jno][3][0] = vol*dldx[jno][0]*0.25;
      eCmat[ino][jno][3][1] = vol*dldx[jno][1]*0.25;
      eCmat[ino][jno][3][2] = vol*dldx[jno][2]*0.25;
      eCmat[ino][jno][0][3] = vol*dldx[ino][0]*0.25;
      eCmat[ino][jno][1][3] = vol*dldx[ino][1]*0.25;
      eCmat[ino][jno][2][3] = vol*dldx[ino][2]*0.25;
      eCmat[ino][jno][3][3] = vol*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }

  double eMmat[4][4][4][4];
  for(int i=0;i<nno*nno*ndim*ndim;++i){ (&eMmat[0][0][0][0])[i] = 0.0; }
  AddEmatConsistentMassTet<4>(
      eMmat,
      rho*vol);
  
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = vol*g_x*0.25;
    eres[ino][1] = vol*g_y*0.25;
    eres[ino][2] = vol*g_z*0.25;
    eres[ino][3] = 0.0;
  }
  
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno*ndim*ndim;++i){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eCmat[0][0][0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<ndim;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<ndim;jdim++){
          eres[ino][idim]
          -= eCmat[ino][jno][idim][jdim]*(velo_press[jno][jdim]+dt_timestep*acc_apress[jno][jdim])
          + eMmat[ino][jno][idim][jdim]*acc_apress[jno][jdim];
        }
      }
    }
  }
}
