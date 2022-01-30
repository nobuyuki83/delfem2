//
// Created by Nobuyuki Umetani on 2022/01/30.
//
/**
 * Singular Value Decomposition using Jacobi method
 */

#ifndef SVD3_H_
#define SVD3_H_

#include <tuple>

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

}  // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/svd3.cpp"
#endif

#endif //SVD3_H_
