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

} // cuda
} // delfem2