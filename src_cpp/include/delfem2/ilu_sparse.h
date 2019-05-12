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
  void Initialize_ILU0(const CMatrixSquareSparse& m);
  void Initialize_ILUk(const CMatrixSquareSparse& m, int fill_level);
  void SetValueILU(const CMatrixSquareSparse& m);
  void Solve(std::vector<double>& vec) const{
		this->ForwardSubstitution(vec);    
		this->BackwardSubstitution(vec);
  }
  bool DoILUDecomp();
private:
  void ForwardSubstitution(  std::vector<double>& vec ) const;
  void BackwardSubstitution( std::vector<double>& vec ) const;
public:
  CMatrixSquareSparse mat;
  int* m_diaInd;
};


std::vector<double> Solve_PCG(double* r_vec,
                              double* u_vec,
                              double conv_ratio,
                              int iteration,
                              const CMatrixSquareSparse& mat,
                              const CPreconditionerILU& ilu);

std::vector<double> Solve_PBiCGStab(double* r_vec,
                                    double* x_vec,
                                    double conv_ratio,
                                    int num_iter,
                                    const CMatrixSquareSparse& mat,
                                    const CPreconditionerILU& ilu);

#endif 
