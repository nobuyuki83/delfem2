/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/bv.h"
#include "delfem2/bvh.h"

namespace delfem2 {
namespace cuda{

void cuda_Min3Max3_Points3F(
    float *p_min3,
    float *p_max3,
    const float *p_XYZ,
    unsigned int np);

// -------------------------------------------------------------

void cuda_CentsMaxRad_MeshTri3F(
    float* pXYZ_cent,
    float* max_rad,
    const float *aXYZ,
    const unsigned int nXYZ,
    const unsigned int *aTri,
    const unsigned int nTri);

void cuda_MortonCode_Points3FSorted(
    unsigned int *aSortedId,
    std::uint32_t *aSortedMc,
    const float *aXYZ,
    const unsigned int nXYZ,
    const float* hMinXYZ,
    const float* hMaxXYZ);

void cuda_MortonCode_BVHTopology(
    CNodeBVH2* aNodeBVH,
    const unsigned int* aSortedId,
    const std::uint32_t* aSortedMc,
    unsigned int N);

void cuda_BVHGeometry(
    float* aAABB,
    const CNodeBVH2* aNodeBVH,
    const float* aXYZ,
    unsigned int nXYZ,
    const unsigned int* aTri,
    unsigned int nTri);


} // cuda
} // delfem2