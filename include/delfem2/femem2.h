/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMEM2_H
#define DFM2_FEMEM2_H

#include "delfem2/dfm2_inline.h"
#include <complex>

namespace delfem2 {

/**
 * @brief derivative of a shape function of a triangle and constant compornent
 */
void TriDlDx(
    double dldx[3][2],
    double const_term[3],
    const double p0[2],
    const double p1[2],
    const double p2[2]);

// -------------------------------------------------------------------
// below: fem element matrix for 2D mesh

// -[\alpha]\nabla^2[value] = [source]
void EMat_Poisson_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha, const double source,
    const double coords[3][2],
    const double value[3]);

void EMat_Helmholtz_Tri2D(
    std::complex<double> eres[3],
    std::complex<double> emat[][3],
    const double wave_length,
    const double coords[3][2],
    const std::complex<double> value[3]);

void EMat_SommerfeltRadiationBC_Line2D(
    std::complex<double> eres[2],
    std::complex<double> emat[2][2],
    double wave_length,
    const double P[2][2],
    const std::complex<double> val[2]);

// [\rho][velo] - [\alpha]\nabla^2[value] = [source]
void EMat_Diffusion_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha,
    const double source,
    const double dt,
    const double gamma,
    const double rho,
    const double coords[3][2],
    const double value[3],
    const double velo[3]);

// ------------------------
// below: linear solid

void EMat_SolidStaticLinear_Tri2D(
    double eres[3][2],
    double emat[][3][2][2],
    const double myu,
    const double lambda,
    const double rho,
    const double g_x,
    const double g_y,
    const double disp[3][2],
    const double coords[3][2]);

void EMat_SolidDynamicLinear_Tri2D(
    double eres[3][2],
    double emat[3][3][2][2],
    const double myu,
    const double lambda,
    const double rho,
    const double g_x,
    const double g_y,
    const double dt_timestep,
    const double gamma_newmark,
    const double beta_newmark,
    const double disp[3][2],
    const double velo[3][2],
    const double acc[3][2],
    const double coords[3][2],
    bool is_initial);

void EMat_Stokes2D_Static_P1(
    double alpha, double g_x, double g_y,
    const double coords[][2],
    const double velo_press[3][3],
    double emat[3][3][3][3],
    double eres[3][3]);

void MakeMat_Stokes2D_Static_P1P1(
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
    double eres_p[3]);

void EMat_Stokes2D_Dynamic_P1(
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
    double eres[3][3]);

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

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femem2.cpp"
#endif
  
#endif /* fem_ematrix_h */
