#ifndef DFM2_FEMHYPER_H
#define DFM2_FEMHYPER_H

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

void AddWdWddW_SolidHyper3Hex(
    double &W,
    double dW[8][3],
    double ddW[8][8][3][3],
    double &vol,
    double c1,
    double c2,
    const double aP0[8][3],
    const double aU[8][3],
    unsigned int iGauss);

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femsolidhyper.cpp"
#endif

#endif