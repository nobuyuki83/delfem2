#ifndef DFM2_FEMSOLIDHYPER_H
#define DFM2_FEMSOLIDHYPER_H

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
 * @param[in] stiff_c1 first parameter of the 2nd order Mooney-Rivlin material
 * @param[in] stiff_c2 second parameter of the 2nd order Mooney-Rivlin material
 * @param[in] aP0 coordinates of vertices of the hex element
 * @param[in] aU displacement of vertices of the hex element
 * @param[in] iGauss degree of Gaussian quadrature
 */
void AddWdWddW_Solid3HyperMooneyrivlin2Reduced_Hex(
    double &W,
    double dW[8][3],
    double ddW[8][8][3][3],
    double &vol,
    double stiff_c1,
    double stiff_c2,
    const double aP0[8][3],
    const double aU[8][3],
    unsigned int iGauss);

/**
 *
 * @param[out] W energy density
 * @param[out] dW derivertive of energy density w.r.t. right cauchy green tensor
 * @param[out] ddW hessian of energy density w.r.t. right cauchy green tensor
 * @param[out] vol volume
 * @param[in] stiff_comp stiffness for compression
 * @param[in] aP0 eight corner vertex positions of a hex element rest shape
 * @param[in] aU eight displacement
 * @param[in] iGauss degree of gaussian quadrature
 */
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

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/femsolidhyper.cpp"
#endif

#endif
