/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMEM3_H
#define DFM2_FEMEM3_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/femutil.h"
#include <complex>

namespace delfem2 {

void EMat_Poisson_Tet3D(double eres[4],
    double emat[][4],
    const double alpha, const double source,
    const double coords[4][3],
    const double value[4]);

void EMat_Diffusion_Newmark_Tet3D(double eres[4],
    double emat[4][4],
    const double alpha, const double source,
    const double dt_timestep, const double gamma_newmark, const double rho,
    const double coords[4][3],
    const double value[4], const double velo[4]);

DFM2_INLINE void MakeCurvetureDKT(
    double B1[][3],
    double B2[][3][2],
    const double coord0[],
    const double coord1[],
    const double coord2[],
    const double l1,
    const double l2 );

DFM2_INLINE void MakeMat_PlateBendingDKT
 (double emat_ww[3][3],
  double emat_wr[3][3][2],
  double emat_rw[3][3][2],
  double emat_rr[3][3][2][2],
  double eres_w[3],
  double eres_r[3][2],
  const double young, const double poisson, const double thickness,
  const double coord[][2], const double w[], const double rot[][2]);

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femem3.cpp"
#endif
  
#endif /* fem_ematrix_h */
