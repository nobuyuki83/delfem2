/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details 
 */

#ifndef DFM2_CLUSTERPOINTS_H
#define DFM2_CLUSTERPOINTS_H
#include "delfem2/dfm2_inline.h"

namespace delfem2{

DFM2_INLINE unsigned int BinaryClustering_Points2d(
    unsigned int* map01,
    //
    const unsigned int np0,
    const double* aArea0,
    const unsigned int* psup_ind0,
    const unsigned int* psup0);

DFM2_INLINE void BinaryClustering_Points3d(
    std::vector<double>& aXYZ1,
    std::vector<double>& aArea1,
    std::vector<double>& aNorm1,
    std::vector<unsigned int>& map01,
    //
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aArea0,
    const std::vector<double>& aNorm0,
    const std::vector<unsigned int>& psup_ind0,
    const std::vector<unsigned int>& psup0);

DFM2_INLINE void BinaryClusteringPoints_FindConnection(
    std::vector<unsigned int>& psup_ind1,
    std::vector<unsigned int>& psup1,
    //
    unsigned int np1,
    unsigned int np0,
    const unsigned int* map01,
    const unsigned int* psup_ind0,
    const unsigned int* psup0);


}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/clusterpoints.cpp"
#endif


#endif /* DFM2_CLUSTERPOINTS */





