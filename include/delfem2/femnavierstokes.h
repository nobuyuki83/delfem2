/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMNAVIERSTOKES_H
#define DFM2_FEMNAVIERSTOKES_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/femutil.h"

#ifdef DFM2_STATIC_LIBRARY
// Merge use explicitly use the template so for static library we need to include the template itself.
#  include "delfem2/lsmats.h"
#endif

namespace delfem2 {

DFM2_INLINE void MakeMat_Stokes2D_Dynamic_Newmark_P1P1(
    double alpha,
    double rho,
    double g_x,
    double g_y,
    const double dt_timestep,
    const double gamma_newmark,
    const double coords[3][2],
    const double velo[3][2],
    const double press[3],
    const double acc[3][2],
    const double apress[3],
    double emat_uu[3][3][2][2],
    double emat_up[][3][2],
    double emat_pu[][3][2],
    double emat_pp[][3],
    double eres_u[3][2],
    double eres_p[3]);

DFM2_INLINE void MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1(
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
    double eres_p[3]);

void EMat_NavierStokes2D_NonStatic_Newmark_P1P1(
    double dt,
    double gamma,
    double rho,
    double myu,
    double g_x,
    double g_y,
    const double coords[][2],
    const double velo[][2],
    const double acc[][2],
    const double press[],
    const double apress[],
    double eres_u[3][2],
    double eres_p[3],
    double eCmat_uu[][3][2][2],
    double eCmat_up[][3][2],
    double eCmat_pu[][3][2],
    double eCmat_pp[][3],
    double eMmat_uu[][3][2][2],
    double eMmat_pu[][3][2]);

void EMat_NavierStokes2D_Dynamic_P1(
    double myu,
    double rho,
    double g_x,
    double g_y,
    const double dt_timestep,
    const double gamma_newmark,
    const double coords[][2],
    const double velo_press[3][3],
    const double acc_apress[3][3],
    double emat[3][3][3][3],
    double eres[3][3]);    

void MakeMat_NavierStokes3D_Dynamic_P1(
    double myu, double rho, double g_x, double g_y, double g_z,
    const double dt_timestep, const double gamma_newmark,
    const double coords[4][3],
    const double velo_press[4][4], const double acc_apress[4][4],
    double emat[4][4][4][4],
    double eres[4][4]);

void MakeMat_NavierStokes3D_Dynamic_Newmark_P1P1(
    double rho, double myu, double g_x, double g_y, double g_z,
    double dt, double gamma,
    const double coords[4][3],
    const double velo[4][3], const double press[4], const double acc[4][3], const double apress[4],
    double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
    double eres_u[4][3], double eres_p[4]);

template <class MAT>
void MergeLinSys_NavierStokes2D(
    MAT& mat_A,
    double* vec_b,
    const double myu,
    const double rho,
    const double g_x,
    const double g_y,
    const double dt_timestep,
    const double gamma_newmark,
    const double* aXY1,
    size_t nXY,
    const unsigned int* aTri1,
    size_t nTri,
    const double* aVal, // vx,vy,press
    const double* aDtVal) // ax,ay,apress
{
  const size_t np = nXY;
  std::vector<unsigned int> tmp_buffer(np, UINT_MAX);
  for (unsigned int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData<3,2>(coords, aIP,aXY1);
    double velo[3][3]; FetchData<3,3>(velo, aIP,aVal);
    double acc[3][3]; FetchData<3,3>(acc, aIP,aDtVal);
    //
    double eres[3][3], emat[3][3][3][3];
    EMat_NavierStokes2D_Dynamic_P1(
        myu, rho,  g_x, g_y,
        dt_timestep, gamma_newmark,
        coords, velo, acc,
        emat, eres);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
//    mat_A.Mearge(3, aIP, 3, aIP, 9, &emat[0][0][0][0], tmp_buffer);
    Merge<3,3,3,3,double>(mat_A,aIP,aIP,emat,tmp_buffer);
  }
}

template <class MAT>
void MergeLinSys_NavierStokes3D_Dynamic(
    MAT& mat_A,
    std::vector<double>& vec_b,
    const double myu,
    const double rho,
    const double g_x,
    const double g_y,
    const double g_z,
    const double dt_timestep,
    const double gamma_newmark,
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aVal,
    const std::vector<double>& aVelo)
{
  const unsigned int np = static_cast<unsigned int>(aXYZ.size()/3);
  const unsigned int nDoF = np*4;
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  std::vector<unsigned int> tmp_buffer(np, UINT_MAX);
  for (unsigned int iel = 0; iel<aTet.size()/4; ++iel){
    const unsigned int aIP[4] = {
      aTet[iel*4+0],
      aTet[iel*4+1],
      aTet[iel*4+2],
      aTet[iel*4+3] };
    double coords[4][3]; FetchData<4,3>(coords, aIP,aXYZ.data());
    double velo_press[4][4]; FetchData<4,4>(velo_press, aIP,aVal.data());
    double acc_apress[4][4]; FetchData<4,4>(acc_apress, aIP,aVelo.data());
    double eres[4][4], emat[4][4][4][4];
    MakeMat_NavierStokes3D_Dynamic_P1(
        myu, rho,  g_x, g_y,g_z,
        dt_timestep, gamma_newmark,
        coords, velo_press, acc_apress,
        emat, eres);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*4+0] += eres[ino][0];
      vec_b[ip*4+1] += eres[ino][1];
      vec_b[ip*4+2] += eres[ino][2];
      vec_b[ip*4+3] += eres[ino][3];
    }
//    mat_A.Mearge(4, aIP, 4, aIP,16, &emat[0][0][0][0], tmp_buffer);
    Merge<4,4,4,4,double>(mat_A,aIP,aIP,emat,tmp_buffer);
  }
}



} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/femnavierstokes.cpp"
#endif
  
#endif /* DFM2_FEMNAVIERSTOKES_H */
