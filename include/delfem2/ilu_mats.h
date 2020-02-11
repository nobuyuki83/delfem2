/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_ILU_MATS
#define DFM2_ILU_MATS

#include <iostream>

#include "mats.h"

namespace delfem2 {

template <typename T>
class CPreconditionerILU
{
public:
  CPreconditionerILU(){}
  CPreconditionerILU(const CPreconditionerILU&); // copy
  ~CPreconditionerILU(){ m_diaInd.clear(); }
  void Initialize_ILU0(const CMatrixSparse<T>& m);
  void Initialize_ILUk(const CMatrixSparse<T>& m, int fill_level);
  void SetValueILU(const CMatrixSparse<T>& m);
  void Solve(std::vector<T>& vec) const{
		this->ForwardSubstitution(vec);
		this->BackwardSubstitution(vec);
  }
  bool DoILUDecomp();
private:
  void ForwardSubstitution(  std::vector<T>& vec ) const;
  void BackwardSubstitution( std::vector<T>& vec ) const;
public:
  CMatrixSparse<T> mat;
  std::vector<unsigned int> m_diaInd;
};
 
  
} // end namespace delfem2

#endif 
