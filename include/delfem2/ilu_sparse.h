//
//  ilu_sparse.h
//  internal_cloth_sparse
//
//  Created by Nobuyuki Umetani on 11/18/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef __internal_cloth_sparse__ilu_sparse__
#define __internal_cloth_sparse__ilu_sparse__

#include <iostream>

#include "matrix_sparse.h"

class CPreconditionerILU
{
public:
  CPreconditionerILU();
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


void Solve_PCG
(double& conv_ratio,
 int& iteration,
 const CMatrixSquareSparse& mat,
 const CPreconditionerILU& ilu,
 std::vector<double>& r_vec,
 std::vector<double>& u_vec);

void Solve_PBiCGSTAB
(double& conv_ratio, int& num_iter,
 const CMatrixSquareSparse& mat,
 const CPreconditionerILU& ilu,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec);

#endif /* defined(__internal_cloth_sparse__ilu_sparse__) */
