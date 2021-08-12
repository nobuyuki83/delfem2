/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_LSILU_MATS_H
#define DFM2_LSILU_MATS_H

#include <iostream>

#include "delfem2/lsmats.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 * @brief ILU decomposision preconditioner class
 */
template<typename T>
class CPreconditionerILU {
 public:
  CPreconditionerILU() noexcept {}
  CPreconditionerILU(const CPreconditionerILU &); // copy
  ~CPreconditionerILU() { m_diaInd.clear(); }
  void Clear() {
    colInd.clear();
    rowPtr.clear();
    valCrs.clear();
    valDia.clear();
    m_diaInd.clear();
  }
  void SetPattern0(const CMatrixSparse<T> &m);
  void Initialize_ILUk(const CMatrixSparse<T> &m, int fill_level);

  void CopyValue(const CMatrixSparse<T> &m);
  void SolvePrecond(T *vec) const {
    this->ForwardSubstitution(vec);
    this->BackwardSubstitution(vec);
  }
  bool Decompose();
  //
  void ForwardSubstitution(T *vec) const;
  void BackwardSubstitution(T *vec) const;

  // treat 1x1 block space matrix as N*N block sparse matrix where the block matrix is diagonal
  void ForwardSubstitutionDegenerate(T *vec, unsigned int N) const;

  // treat 1x1 block space matrix as N*N block sparse matrix where the block matrix is diagonal
  void BackwardSubstitutionDegenerate(T *vec, unsigned int N) const;
 public:
  unsigned int nblk;
  unsigned int ndim;
  std::vector<unsigned int> colInd;
  std::vector<unsigned int> rowPtr;
  std::vector<unsigned int> m_diaInd;
  std::vector<T> valCrs;
  std::vector<T> valDia;
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/lsilu_mats.cpp"
#endif

#endif 
