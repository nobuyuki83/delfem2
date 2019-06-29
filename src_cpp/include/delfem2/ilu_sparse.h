/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef ILU_SPARSE
#define ILU_SPARSE

#include <iostream>

#include "matrix_sparse.h"

class CPreconditionerILU
{
public:
  CPreconditionerILU();
  CPreconditionerILU(const CPreconditionerILU&); // copy
  ~CPreconditionerILU();
  void Initialize_ILU0(const CMatrixSparse& m);
  void Initialize_ILUk(const CMatrixSparse& m, int fill_level);
  void SetValueILU(const CMatrixSparse& m);
  void Solve(std::vector<double>& vec) const{
		this->ForwardSubstitution(vec);    
		this->BackwardSubstitution(vec);
  }
  bool DoILUDecomp();
private:
  void ForwardSubstitution(  std::vector<double>& vec ) const;
  void BackwardSubstitution( std::vector<double>& vec ) const;
public:
  CMatrixSparse mat;
  std::vector<unsigned int> m_diaInd;
};


std::vector<double> Solve_PCG(double* r_vec,
                              double* u_vec,
                              double conv_ratio,
                              int iteration,
                              const CMatrixSparse& mat,
                              const CPreconditionerILU& ilu);

std::vector<double> Solve_PBiCGStab(double* r_vec,
                                    double* x_vec,
                                    double conv_ratio,
                                    int num_iter,
                                    const CMatrixSparse& mat,
                                    const CPreconditionerILU& ilu);

#endif 
