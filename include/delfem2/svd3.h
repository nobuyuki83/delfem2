//
// Created by Nobuyuki Umetani on 2022/01/30.
//
/**
 * Singular Value Decomposition using Jacobi method
 */

#ifndef SVD3_H_
#define SVD3_H_

#include <tuple>
#include <array>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 * @func compute eigen value & vector for symmmetric matrix
 * @details
 * sm[6] = (M_00,M_11,M_22,M_12,M_20,M_01)
 * M = ULU^T
 * u[9] = (U_00,U_01,U_02, U_10,U_11,U_12, U_20,U_21,U_22)
 */
DFM2_INLINE bool EigenSym3(
    double u[9],
    double l[3],
    const double sm[6],
    int nitr);

// m = UGV^T
DFM2_INLINE void Svd3(
    double U[9],
    double G[3],
    double V[9],
    const double m[9],
    int nitr);

template <typename MAT>
std::tuple<MAT,MAT,MAT> Svd3(
    const MAT& m,
    int nitr){
  const double m_[9] = {
      m(0,0),m(0,1),m(0,2),
      m(1,0),m(1,1),m(1,2),
      m(2,0),m(2,1),m(2,2)};
  double U[9], G[3], V[9];
  Svd3(
      U,G,V,
      m_, nitr);
  return {
      {U[0],U[1],U[2], U[3],U[4],U[5], U[6],U[7],U[8]},
      {G[0],0,0,  0,G[1],0, 0,0,G[2]},
      {V[0],V[1],V[2], V[3],V[4],V[5], V[6],V[7],V[8]} };
}

DFM2_INLINE void GetRotPolarDecomp(
    double R[9],
    //
    const double am[9],
    int nitr);

/**
 * this is put in the header because MAT can be Eigen::Matrix3
 * @tparam MAT
 * @tparam T
 * @param diff
 * @param U
 * @param S
 * @param V
 */
template<class MAT, typename T = typename MAT::Scalar>
DFM2_INLINE void Svd3Differential(
    T diff[3][3][9],
    const MAT &U,
    const MAT &S,
    const MAT &V) {
  auto invmat2 = [](T a0, T a1) -> std::array<T, 2> {
    if (fabs(a0 - a1) < 1.0e-6) { a0 += 1.0e-6; }
    T detinv = 1 / (a0 * a0 - a1 * a1);
    return {+detinv * a0, -detinv * a1};
  };

  const std::array<T, 2> Ai0 = invmat2(S(1, 1), S(2, 2));
  const std::array<T, 2> Ai1 = invmat2(S(2, 2), S(0, 0));
  const std::array<T, 2> Ai2 = invmat2(S(0, 0), S(1, 1));
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      // dSdu
      diff[i][j][3] = U(i, 0) * V(j, 0);
      diff[i][j][4] = U(i, 1) * V(j, 1);
      diff[i][j][5] = U(i, 2) * V(j, 2);
      {
        const T b0[2] = {-U(i, 2) * V(j, 1), U(i, 1) * V(j, 2)};
        diff[i][j][0] = +Ai0[0] * b0[0] + Ai0[1] * b0[1];
        diff[i][j][6] = -Ai0[1] * b0[0] - Ai0[0] * b0[1];
      }
      {
        const T b1[2] = {-U(i, 0) * V(j, 2), U(i, 2) * V(j, 0)};
        diff[i][j][1] = +Ai1[0] * b1[0] + Ai1[1] * b1[1];
        diff[i][j][7] = -Ai1[1] * b1[0] - Ai1[0] * b1[1];
      }
      {
        const T b2[2] = {-U(i, 1) * V(j, 0), U(i, 0) * V(j, 1)};
        diff[i][j][2] = +Ai2[0] * b2[0] + Ai2[1] * b2[1];
        diff[i][j][8] = -Ai2[1] * b2[0] - Ai2[0] * b2[1];
      }
    }
  }
}

}  // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/svd3.cpp"
#endif

#endif //SVD3_H_
