/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMEM2_H
#define DFM2_FEMEM2_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/femutil.h"
#include <complex>

namespace delfem2 {



// -------------------------------------------------------------------
// below: fem element matrix for 2D mesh

// -[\alpha]\nabla^2[value] = [source]
void EMat_Poisson_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha, const double source,
    const double coords[3][2],
    const double value[3]);

void EMat_Poission2_QuadOrth(
    double emat[4][4],
    double lx,
    double ly);

/**
 *
 * @param emat
 * @param lx
 * @param ly
 * @param[in] ngauss ngauss=1 is enough for analytically exact integration
 */
void EMat_Poisson2_QuadOrth_GaussInt(
    double emat[4][4],
    double lx,
    double ly,
    unsigned int ngauss);

void EMat_SolidLinear2_QuadOrth_GaussInt(
    double emat[4][4][2][2],
    double lx,
    double ly,
    double myu,
    double lambda,
    unsigned int ngauss);

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


} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femem2.cpp"
#endif
  
#endif /* fem_ematrix_h */
