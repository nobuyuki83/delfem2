#ifndef __internal_cloth_sparse__ilu_sparse__
#define __internal_cloth_sparse__ilu_sparse__

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


void Solve_PCG
(double* r_vec,
 double* u_vec,
 double conv_ratio,
 int iteration,
 const CMatrixSquareSparse& mat,
 const CPreconditionerILU& ilu);

void Solve_PBiCGSTAB
(double& conv_ratio, int& num_iter,
 const CMatrixSquareSparse& mat,
 const CPreconditionerILU& ilu,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec);

#endif /* defined(__internal_cloth_sparse__ilu_sparse__) */
