/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/fem_rod2.h"

#include "delfem2/geo3_v23m34q.h"

DFM2_INLINE void delfem2::WdWddW_Rod2(
    double &W,
    double dW[3][2],
    double ddW[3][3][2][2],
    const double ap[3][2],
    const double aL[2],
    double stiff_stretch01,
    double stiff_stretch12,
    double stiff1_bend) {
  constexpr unsigned nP = 3;
  const double lF = Distance2(ap[0], ap[1]);
  const double lG = Distance2(ap[1], ap[2]);
  const double lFi = 1.0 / lF;
  const double lGi = 1.0 / lG;
  const double F[2] = {
      (ap[1][0] - ap[0][0]) * lFi,
      (ap[1][1] - ap[0][1]) * lFi};
  const double G[2] = {
      (ap[2][0] - ap[1][0]) * lGi,
      (ap[2][1] - ap[1][1]) * lGi};
  const double cFG = F[0] * G[0] + F[1] * G[1];
  const double sFG = F[0] * G[1] - F[1] * G[0];
  const double dF[2][nP][2] = {
      {
          {-(1. - F[0] * F[0]) * lFi, +F[0] * F[1] * lFi},
          {+(1. - F[0] * F[0]) * lFi, -F[0] * F[1] * lFi},
          {0., 0.},
      },
      {
          {+F[1] * F[0] * lFi, -(1. - F[1] * F[1]) * lFi},
          {-F[1] * F[0] * lFi, +(1. - F[1] * F[1]) * lFi},
          {0., 0.},
      }
  };
  const double dG[2][nP][2] = {
      {
          {0., 0.},
          {-(1. - G[0] * G[0]) * lGi, +G[0] * G[1] * lGi},
          {+(1. - G[0] * G[0]) * lGi, -G[0] * G[1] * lGi},
      },
      {
          {0., 0.},
          {+G[0] * G[1] * lGi, -(1. - G[1] * G[1]) * lGi},
          {-G[0] * G[1] * lGi, +(1. - G[1] * G[1]) * lGi}
      }
  };
  //
  constexpr unsigned int nT = 3;
  const double aT[nT] = {
      lF - aL[0],
      lG - aL[1],
      sFG / (1 + cFG)};
  const double aWT[nT] = {
      stiff_stretch01,
      stiff_stretch12,
      stiff1_bend};
  double dT[nT][nP][2];
  std::fill(&dT[0][0][0], &dT[0][0][0] + nT * nP * 2, 0.);
  dT[0][0][0] = -F[0];
  dT[0][0][1] = -F[1];
  dT[0][1][0] = +F[0];
  dT[0][1][1] = +F[1];
  dT[1][1][0] = -G[0];
  dT[1][1][1] = -G[1];
  dT[1][2][0] = +G[0];
  dT[1][2][1] = +G[1];
  for (int ip = 0; ip < 3; ++ip) {
    for (int idim = 0; idim < 2; ++idim) {
      dT[2][ip][idim] =
          1.0 / (1. + cFG) *
          ( +
          dF[0][ip][idim] * G[1] +
          F[0] * dG[1][ip][idim] -
          dF[1][ip][idim] * G[0] -
          F[1] * dG[0][ip][idim]
          ) -
          sFG / ((1. + cFG) * (1. + cFG)) *
          ( +
          dF[0][ip][idim] * G[0] +
          F[0] * dG[0][ip][idim] +
          dF[1][ip][idim] * G[1] +
          F[1] * dG[1][ip][idim]
          );
    }
  }
  /*
  W = aT[2];
  for(int ip=0;ip<nP;++ip) {
    for (int idim = 0; idim < 2; ++idim) {
      dW[ip][idim] = dT[2][ip][idim];
    }
  }
  return;
   */
  W = 0.0;
  for (unsigned int iT = 0; iT < nT; ++iT) {
    W += 0.5 * aWT[iT] * aT[iT] * aT[iT];
  }
  for (unsigned int ip = 0; ip < nP; ++ip) {
    for (int idim = 0; idim < 2; ++idim) {
      dW[ip][idim] = 0.0;
      for (unsigned int iT = 0; iT < nT; ++iT) {
        dW[ip][idim] += aWT[iT] * dT[iT][ip][idim] * aT[iT];
      }
    }
  }
  for (unsigned int ip = 0; ip < nP; ++ip) {
    for (unsigned int jp = 0; jp < nP; ++jp) {
      for (int idim = 0; idim < 2; ++idim) {
        for (int jdim = 0; jdim < 2; ++jdim) {
          ddW[ip][jp][idim][jdim] = 0.0;
          for (unsigned int iT = 0; iT < nT; ++iT) {
            ddW[ip][jp][idim][jdim] += aWT[iT] * dT[iT][ip][idim] * dT[iT][jp][jdim];
          }
        }
      }
    }
  }
}