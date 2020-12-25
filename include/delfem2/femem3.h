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

void ddW_MassConsistentVal3D_Tet3D(
    double* eMmat,
    double rho, double vol,
    bool is_add,
    unsigned int nstride = 3);

/**
 * @brief compute energy and its 1st and 2nd derivative for cloth bending
 * @param[out] W strain energy
 * @param[out] dW 1st derivative of energy
 * @param[out] ddW 2nd derivative of energy
 * @param[in] C undeformed quad vertex positions
 * @param[in] c deformed quad vertex positions
 */
void WdWddW_Bend(
    double& W,
    double dW[4][3],
    double ddW[4][4][3][3],
    //
    const double C[4][3],
    const double c[4][3],
    double stiff);

void WdWddW_CST(
    double& W, // (out) energy
    double dW[3][3], // (out) 1st derivative of energy
    double ddW[3][3][3][3], // (out) 2nd derivative of energy
    //
    const double C[3][3], // (in) undeformed triangle vertex positions
    const double c[3][3], // (in) deformed triangle vertex positions
    const double lambda, // (in) Lame's 1st parameter
    const double myu);   // (in) Lame's 2nd parameter

/**
 * @brief compute energy and its 1st and 2nd derivative for contact against object
 */
class CInput_Contact
{
public:
  virtual double penetrationNormal(
      double& nx, double& ny, double& nz,
      double px, double py, double pz) const = 0;
};

void WdWddW_Contact(
    double& W,  // (out) energy
    double dW[3], // (out) 1st derivative of energy
    double ddW[3][3], // (out) 2nd derivative of energy
    //
    const double c[3], // (in) deformed triangle vertex positions
    double stiff_contact,
    double contact_clearance,
    const CInput_Contact& input);

// above: 2D
// --------------------------------------------------------------------------
// below: 3D

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
