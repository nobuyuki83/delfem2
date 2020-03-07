/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFEM2_MATN_H
#define DFEM2_MATN_H

/**
 * @file template class/function for vector/matrix where size is determined at compiling
 * @detail this header's extension is .hpp because it is purely template header
 */

namespace delfem2 {

template <typename REAL, unsigned int N>
void MatMat
 (REAL* C,
  const REAL* A, const REAL* B)
{
  for(unsigned int i=0;i<N;++i){
    for(unsigned int j=0;j<N;++j){
      C[i*N+j] = 0;
      for(unsigned int k=0;k<N;++k){
        C[i*N+j] += A[i*N+k]*B[k*N+j];
      }
    }
  }
}


template <typename REAL, unsigned int NCOL, unsigned int NROW>
void MatVec
(REAL* y,
 const REAL* A, const REAL* x)
{
  for(unsigned int i=0;i<NCOL;++i){
    y[i] = 0;
    for(unsigned int j=0;j<NROW;++j){
      y[i] += A[i*NROW+j]*x[j];
    }
  }
}


}

#endif /* DFM2_VM_H */
