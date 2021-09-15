/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @discussion splitting this file into "vecx.h" and "itersol.h" in the future
 */

#ifndef DFM2_LS_MASTERSLAVE_H
#define DFM2_LS_MASTERSLAVE_H

#include <vector>
#include <cassert>
#include <complex>
#include <iostream>

#include "delfem2/dfm2_inline.h"
#include "delfem2/lsmats.h"

namespace delfem2 {

DFM2_INLINE void setRHS_MasterSlave(
    double *vec_b,
    unsigned int nDoF,
    const unsigned int *aMSFlag);

DFM2_INLINE void JArray_AddMasterSlavePattern(
    std::vector<unsigned int> &index,
    std::vector<unsigned int> &array,
    const unsigned int* aMSFlag,
    unsigned int ndim,
    const unsigned int *psup_ind0,
    size_t npsup_ind0,
    const unsigned int *psup0);

DFM2_INLINE void SetMasterSlave(
    delfem2::CMatrixSparse<double> &mat,
    const unsigned int *aMSFlag);

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/ls_masterslave.cpp"
#endif
  
#endif // MATDIA_CRS_H
