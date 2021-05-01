#ifndef DFM2_FEMHYPER_H
#define DFM2_FEMHYPER_H

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 * Compute elastic potential and its grandient (residual vector) and hessian (stiffness matrix)
 * for the 2nd order Mooney-Rivlin hyper elastic material with reduced invariants
 * for hex element
 * @param[out] W elastic potential
 * @param[out] dW gradient of W
 * @param[out] ddW hessian of W
 * @param[out] vol volume
 * @param[in] c1 first parameter of the 2nd order Mooney-Rivlin material
 * @param[in] c2 second parameter of the 2nd order Mooney-Rivlin material
 * @param[in] aP0 coordinates of vertices of the hex element
 * @param[in] aU displacement of vertices of the hex element
 * @param[in] iGauss degree of Gaussian quadrature
 */
void AddWdWddW_Solid3HyperMooneyrivlin2Reduced_Hex(
    double &W,
    double dW[8][3],
    double ddW[8][8][3][3],
    double &vol,
    double c1,
    double c2,
    const double aP0[8][3],
    const double aU[8][3],
    unsigned int iGauss);

void AddWdWddW_Solid3Compression_Hex(
    double& W,
    double dW[8][3],
    double ddW[8][8][3][3],
    double& vol,
    //
    double stiff_comp,
    const double aP0[8][3],
    const double aU[8][3],
    unsigned int iGauss);

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femsolidhyper.cpp"
#endif

#endif
