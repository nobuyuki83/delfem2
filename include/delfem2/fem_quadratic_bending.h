//
// Created by Nobuyuki Umetani on 2021-11-20.
//

#ifndef FEM_QUADRATIC_BENDING_H_
#define FEM_QUADRATIC_BENDING_H_

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE void CdC_QuadBend(
    double C[3],
    double dCdp[3][4][3],
    const double P[4][3],
    const double p[4][3]);

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

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fem_quadratic_bending.cpp"
#endif

#endif //FEM_QUADRATIC_BENDING_H_
