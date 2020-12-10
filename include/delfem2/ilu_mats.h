/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_ILU_MATS_H
#define DFM2_ILU_MATS_H

#include "mats.h"
#include "delfem2/dfm2_inline.h"
#include <iostream>

namespace delfem2 {

/**
 * @brief ILU decomposision preconditioner class
 */
template <typename T>
class CPreconditionerILU
{
public:
  CPreconditionerILU(){}
  CPreconditionerILU(const CPreconditionerILU&); // copy
  ~CPreconditionerILU(){ m_diaInd.clear(); }
  void Clear(){
    mat.Clear();
    m_diaInd.clear();
  }
  void Initialize_ILU0(const CMatrixSparse<T>& m);
  void Initialize_ILUk(const CMatrixSparse<T>& m, int fill_level);
  void SetValueILU(const CMatrixSparse<T>& m);
  void SolvePrecond(T* vec) const{
		this->ForwardSubstitution(vec);
		this->BackwardSubstitution(vec);
  }
  bool DoILUDecomp();
  //
  void ForwardSubstitution(  T* vec ) const;
  void BackwardSubstitution( T* vec ) const;
  
  // treat 1x1 block space matrix as N*N block sparse matrix where the block matrix is diagonal
  void ForwardSubstitutionDegenerate(  T* vec, unsigned int N) const;
  
  // treat 1x1 block space matrix as N*N block sparse matrix where the block matrix is diagonal
  void BackwardSubstitutionDegenerate( T* vec, unsigned int N ) const;
public:
  CMatrixSparse<T> mat;
  std::vector<unsigned int> m_diaInd;
};
   
} // end namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/ilu_mats.cpp"
#endif

#endif 
