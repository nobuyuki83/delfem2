/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include <vector>

class CMatrixSparse
{
public:
	CMatrixSparse();
	virtual ~CMatrixSparse();
  
  void Initialize(int nblk, int len, bool is_dia);
  void operator = (const CMatrixSparse& m);
  void SetPattern(const int* colind, unsigned int ncolind,
                  const int* rowptr, unsigned int nrowptr);

	bool SetZero();
	bool Mearge(unsigned int nblkel_col, const unsigned int* blkel_col,
              unsigned int nblkel_row, const unsigned int* blkel_row,
              unsigned int blksize, const double* emat,
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
  void ScaleLeftRight(const double* scale);
  void AddDia(double eps);
public:
	unsigned int nblk_col;
  unsigned int nblk_row;
  unsigned int len_col;
  unsigned int len_row;
  
	unsigned int  ncrs;
  std::vector<unsigned int> colInd;
  std::vector<unsigned int> rowPtr;
  
	double* valCrs;
  
  bool is_dia;
	double* valDia;
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

void Solve_CG(double& conv_ratio,
              int& iteration,
              const CMatrixSparse& mat,
              std::vector<double>& r_vec,
              std::vector<double>& u_vec);

bool Solve_BiCGSTAB(double& conv_ratio,
                    int& num_iter,
                    const CMatrixSparse& mat,
                    std::vector<double>& r_vec,
                    std::vector<double>& x_vec);


#endif // MATDIA_CRS_H
