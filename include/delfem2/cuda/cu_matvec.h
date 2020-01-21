/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

namespace delfem2 {
namespace cuda{

void cuda_VecScale(
    float *hOut,
    const float *hIn,
    float scale,
    const unsigned int n);

// -------------------------------------------------------------

float cuda_Dot(
    const float *h_A,
    const float *h_B,
    unsigned int n);

// -------------------------------------------------------------

void cuda_MatMat(
    float *h_C_gpu,
    const float *h_A,
    const float *h_B,
    unsigned int WIDTH);

// -------------------------------------------------------------

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

void cuda_MortonCode_Points3F(
    unsigned int *aSortedId,
    std::uint32_t *aSortedMc,
    const float *aXYZ,
    const unsigned int nXYZ);


} // cuda
} // delfem2