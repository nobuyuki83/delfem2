/**
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include <vector>

class CMatrixSquareSparse
{
public:
	CMatrixSquareSparse();
	virtual ~CMatrixSquareSparse();
  
  void Initialize(int nblk, int len, bool is_dia);
  void operator = (const CMatrixSquareSparse& m);
  void SetPattern(const int* colind, int ncolind,
                  const int* rowptr, int nrowptr);

	bool SetZero();
	bool Mearge(int nblkel_col, const int* blkel_col,
              int nblkel_row, const int* blkel_row,
              int blksize, const double* emat,
              std::vector<int>& m_marge_tmp_buffer);
  // Calc Matrix Vector Product
  // {y} = alpha * [A]{x} + beta * {y}  
	void MatVec(double alpha,
              const std::vector<double>& x,
              double beta,
              std::vector<double>& y) const;
  void SetBoundaryCondition(const int* pBCFlag, int nP, int ndimVal);
  void SetMasterSlave(const int* aMSFlag);
  double CheckSymmetry() const;
public:
	int m_nblk_col;
	int m_len_col;
  int m_nblk_row;
  int m_len_row;
  
	int  m_ncrs;
	int* m_colInd;
	int* m_rowPtr;
  
	double* m_valCrs;
  
  bool is_dia;
	double* m_valDia;
};

double InnerProduct(const std::vector<double>& r_vec,
                    const std::vector<double>& u_vec);
double InnerProduct(const double* r_vec,
                    const double* u_vec,
                    int ndof);

void AXPY(double a,
          const std::vector<double>& x,
          std::vector<double>& y);
void AXPY(double a,
          const double* x,
          double* y,
          int n);

void Solve_CG
(double& conv_ratio,
 int& iteration,
 const CMatrixSquareSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& u_vec);

bool Solve_BiCGSTAB
(double& conv_ratio,
 int& num_iter,
 const CMatrixSquareSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec);


#endif // MATDIA_CRS_H
