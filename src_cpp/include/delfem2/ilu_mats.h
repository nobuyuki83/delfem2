/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef ILU_SPARSE
#define ILU_SPARSE

#include <iostream>

#include "mats.h"

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
  void ForwardSubstitution(  std::vector<double>& vec ) const;
  void BackwardSubstitution( std::vector<double>& vec ) const;
public:
  CMatrixSparse<T> mat;
  std::vector<unsigned int> m_diaInd;
};

template <typename T>
std::vector<double> Solve_PCG(double* r_vec,
                              double* u_vec,
                              double conv_ratio,
                              int iteration,
                              const CMatrixSparse<T>& mat,
                              const CPreconditionerILU<T>& ilu);

template <typename T>
std::vector<double> Solve_PBiCGStab(double* r_vec,
                                    double* x_vec,
                                    double conv_ratio,
                                    int num_iter,
                                    const CMatrixSparse<T>& mat,
                                    const CPreconditionerILU<T>& ilu);

template <typename T>
void SolveLinSys_PCG
(const CMatrixSparse<T>& mat_A,
 std::vector<double>& vec_b,
 std::vector<double>& vec_x,
 CPreconditionerILU<T>& ilu_A,
 double& conv_ratio,
 int& iteration)
{
  // set ILU preconditioner
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  // solve linear system
  //Solve_CG(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_BiCGSTAB(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_PBiCGSTAB(conv_ratio, iteration, mat_A, ilu_A, vec_b, vec_x);
  vec_x.resize(vec_b.size());
  Solve_PCG(vec_b.data(), vec_x.data(), conv_ratio, iteration,
            mat_A, ilu_A);
  std::cout<<"  conv_ratio:"<<conv_ratio<<"  iteration:"<<iteration<<std::endl;
}

template <typename T>
bool SolveLinSys_BiCGStab
(CMatrixSparse<T>& mat_A,
 std::vector<double>& vec_b,
 std::vector<double>& vec_x,
 CPreconditionerILU<T>& ilu_A,
 double& conv_ratio,
 int& iteration)
{
  // set ILU preconditioner
  ilu_A.SetValueILU(mat_A);
  bool res_ilu = ilu_A.DoILUDecomp();
  if( !res_ilu ){ return false; }
  // solve linear system
  //  double conv_ratio = 1.0e-4;
  //  int iteration = 1000;
  //  Solve_CG(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_BiCGSTAB(conv_ratio, iteration, mat_A, vec_b, vec_x);
  vec_x.resize(vec_b.size());
  Solve_PBiCGStab(vec_b.data(), vec_x.data(), conv_ratio, iteration, mat_A, ilu_A);
  /// Solve_PCG(conv_ratio, iteration, mat_A, ilu_A, vec_b, vec_x);
  //  std::cout<<"  interative solver --- conv_ratio:"<<conv_ratio<<"  iteration:"<<iteration<<std::endl;
  return true;
}

#endif 
