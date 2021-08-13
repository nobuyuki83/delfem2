/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMMITC3_H
#define DFM2_FEMMITC3_H

#include <vector>

#include "delfem2/femutil.h"
#include "delfem2/dfm2_inline.h"

#ifdef DFM2_STATIC_LIBRARY
// Merge use explicitly use the template so for static library we need to include the template itself.
#  include "delfem2/lsmats.h"
#endif

namespace delfem2 {

void WdWddW_PlateBendingMITC3(double &W,
                              double dW[3][3],
                              double ddW[3][3][3][3],
                              const double C[3][2], // initial XY position
                              const double u[3][3], // z displacement + xy axis rotation
                              double thk,
                              double lambda,
                              double myu);

template<class MAT>
void MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(
    MAT &mat_A,
    double *vec_b,
    const double thick,
    const double lambda,
    const double myu,
    const double rho,
    const double gravity_z,
    const double *aXY1,
    size_t nXY,
    const unsigned int *aTri1,
    size_t nTri,
    const double *aVal) {
  const size_t np = nXY;
  std::vector<unsigned int> tmp_buffer(np, UINT_MAX);
  for (unsigned int iel = 0; iel < nTri; ++iel) {
    const unsigned int i0 = aTri1[iel * 3 + 0];
    const unsigned int i1 = aTri1[iel * 3 + 1];
    const unsigned int i2 = aTri1[iel * 3 + 2];
    const unsigned int aIP[3] = {i0, i1, i2};
    double P[3][2];
    delfem2::FetchData<3, 2>(P, aIP, aXY1);
    double u[3][3];
    delfem2::FetchData<3, 3>(u, aIP, aVal);
    //
    double W = 0.0, dW[3][3], ddW[3][3][3][3];
    for (int i = 0; i < 9; ++i) { (&dW[0][0])[i] = 0.0; }
    for (int i = 0; i < 81; ++i) { (&ddW[0][0][0][0])[i] = 0.0; }
    WdWddW_PlateBendingMITC3(
        W, dW, ddW,
        P, u,
        thick, lambda, myu);
    {
      const double A = delfem2::femutil::TriArea2D(P[0], P[1], P[2]);
      dW[0][0] = rho * A * thick / 3.0 * gravity_z;
      dW[1][0] = rho * A * thick / 3.0 * gravity_z;
      dW[2][0] = rho * A * thick / 3.0 * gravity_z;
    }
    for (unsigned int ino = 0; ino < 3; ino++) {
      const unsigned int ip = aIP[ino];
      vec_b[ip * 3 + 0] += dW[ino][0];
      vec_b[ip * 3 + 1] += dW[ino][1];
      vec_b[ip * 3 + 2] += dW[ino][2];
    }
    // marge dde
//    mat_A.Mearge(3, aIP, 3, aIP, 9, &ddW[0][0][0][0], tmp_buffer);
    Merge<3, 3, 3, 3, double>(mat_A, aIP, aIP, ddW, tmp_buffer);
  }
}

void MassLumped_ShellPlateBendingMITC3(
    double *aM,
    double rho, double thick,
    const double *aXY, unsigned int nXY,
    const unsigned int *aTri, unsigned int nTri);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/femmitc3.cpp"
#endif

#endif /* fem_ematrix_h */
