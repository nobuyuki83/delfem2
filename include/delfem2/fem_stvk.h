//
// Created by Nobuyuki Umetani on 2021/11/20.
//

#ifndef DFM2_FEM_STVK_H
#define DFM2_FEM_STVK_H

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 *
 * @param[out] C
 * @param[out] dCdp
 * @param[in] P undeformed triangle vertex positions
 * @param[in] p deformed triangle vertex positions
 */
DFM2_INLINE void CdC_StVK(
    double C[3],
    double dCdp[3][9],
    const double P[3][2],
    const double p[3][3]);

/**
  *
  * @param[out] C energy
  * @param[out] dCdp 1st derivative of energy
  * @param[in] P undeformed triangle vertex positions
  * @param[in] p deformed triangle vertex positions
  * @param[in] lambda Lame's 1st parameter
  * @param[in] myu Lame's 2nd parameter
  */
DFM2_INLINE void CdC_EnergyStVK(
    double &C,
    double dCdp[9],
    const double P[3][2],
    const double p[3][3],
    const double lambda,
    const double myu);

void WdWddW_CST(
    double& W, // (out) energy
    double dW[3][3], // (out) 1st derivative of energy
    double ddW[3][3][3][3], // (out) 2nd derivative of energy
    //
    const double C[3][3], // (in) undeformed triangle vertex positions
    const double c[3][3], // (in) deformed triangle vertex positions
    const double lambda, // (in) Lame's 1st parameter
    const double myu);   // (in) Lame's 2nd parameter

void WdWddW_CST_Sensitivity(
  double Kmat[3][3][3][3], 
  double Res[3][3], 
  double dRdC[3][3][3][3],
  double lambda, 
  double myu,
  const double C[3][3], 
  const double c[3][3]);


}

#ifndef DFM2_STATIC_LIBRARY

#include "delfem2/fem_stvk.cpp"

#endif

#endif // DFM2_FEM_STVK_H
