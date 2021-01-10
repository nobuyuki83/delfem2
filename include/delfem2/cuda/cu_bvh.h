/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/srchbv3aabb.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"

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

void cuda_BVHGeometry_AABB3f(
    CBV3_AABB<float>* aAABB,
    const CNodeBVH2* aNodeBVH,
    const float* aXYZ,
    unsigned int nXYZ,
    const unsigned int* aTri,
    unsigned int nTri);

template <typename REAL>
void cuda_BVHGeometry_Sphere(
    CBV3_Sphere<REAL>* aSphere,
    const CNodeBVH2* aNodeBVH,
    const REAL* aXYZ,
    unsigned int nXYZ,
    const unsigned int* aTri,
    unsigned int nTri);

template <typename REAL>
void cuda_BVH_NearestPoint(
    unsigned int* hInd,
    //
    const REAL* hXYZ1,
    unsigned int nXYZ1,
    const CNodeBVH2* hNodeBVH0,
    unsigned int nNodeBVH0,
    const CBV3_Sphere<REAL>* hBVSphere0);

} // cuda
} // delfem2