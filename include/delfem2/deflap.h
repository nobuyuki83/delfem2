/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_DEFLAP_H
#define DFM2_DEFLAP_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/vec3.h"

// ---------------------------

namespace delfem2 {

class CDef_LaplacianLinearAsym
{
public:
  void Init(
      const std::vector<double>& aXYZ0,
      const std::vector<unsigned int>& aTri);
  void Deform(
      std::vector<double>& aXYZ1,
      const std::vector<double>& aXYZ0,
      const std::vector<int>& aBCFlag);
public:
  CMatrixSparse<double> mat_A;
  std::vector<double> aRhs0, aRhs1;
  std::vector<double> aHistConv;
};

// =====================================

/**
 * @brief deformation classs where the matrix is computed as L^TL
 */
class CDef_LaplacianLinearGram{
public:
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri,
            bool is_preconditioner);
  void Deform(std::vector<double>& aXYZ1,
              const std::vector<double>& aXYZ0) const;
  void SetBoundaryConditionToPreconditioner();
  // -----------
  // called from solver
  void MatVec(double* y,
              double alpha, const double* vec,  double beta) const;
  void SolvePrecond(double* v) const;
public:
  CMatrixSparse<double> Mat;
  bool is_preconditioner;
  double weight_bc = 100.0;
  std::vector<int> aBCFlag;
  std::vector<double> aRes0;
  // preconditioner
  std::vector<double> aDiaInv;
  // temprary vectors for solver
  mutable std::vector<double> vec_tmp0, vec_tmp1, vec_tmp2;
  mutable std::vector<double> aConvHist;
};

// =====================================

/**
 * @brief Laplacian deformation classs
 */
class CDef_LaplacianLinear{
public:
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri,
            bool is_preconditioner);
  void Deform(std::vector<double>& aXYZ1,
              const std::vector<double>& aXYZ0) const;
  void SetValueToPreconditioner();
  // -----------
  void MatVec(
      double* y,
      double alpha, const double* vec,  double beta) const;
public:
  //! penalty coefficient for the sliding boundary condition
  double weight_nrm = 100;
  std::vector< std::pair<unsigned int, CVec3d> > aIdpNrm;

  //! penalty coefficient for the fixed boundary condition
  double weight_bc = 100;
  std::vector<int> aBCFlag;

  //! system matrix of the elastic energy
  CMatrixSparse<double> Mat;

  //! resisual force of the residual force
  std::vector<double> aRes0;

  //! maximum number of iteration for the linear solver
  unsigned int max_itr = 300;

  //! convergence tolerance for the linear solver
  double conv_tol = 1.0e-5;

  //! swich to activate preconditioer
  bool is_preconditioner;

  //! preconditioiner
  CPreconditionerILU<double> Prec;

  // 
  mutable std::vector<double> vec_tmp0, vec_tmp1;

  //! convergence history of the linear solver
  mutable std::vector<double> aConvHist;
};

// =====================================

/**
 * @brief Laplacian deformation classs lightweight system matrix
 * @detail because the 3x3 block matrix in the sparse system matrix is diagonal,
 * we just have non-block sparse matrix.
 * functionality is the same as CDef_LaplacianLinear
 * TODO: probably I can make this faster by introducing ILU preconditioner
 */
class CDef_LaplacianLinearDegenerate{
public:
  void Init(
      const std::vector<double>& aXYZ0,
      const std::vector<unsigned int>& aTri,
      bool is_preconditinoner);
  void Deform(
      std::vector<double>& aXYZ1,
      const std::vector<double>& aXYZ0) const;
  void SetBoundaryConditionToPreconditioner();
  // -----------
  void MatVec(double* y,
      double alpha, const double* vec,  double beta) const;
  void SolvePrecond(double* y) const;
public:
  //! penalty coefficient for the sliding boundary condition
  double weight_nrm = 100;
  //! normal for the sliding boundary condition
  std::vector< std::pair<unsigned int, CVec3d> > aIdpNrm;

  //! penalty coefficient for the fixed boundary condition
  double weight_bc = 100;
  std::vector<int> aBCFlag;

  //! system matrix of the elastic energy
  CMatrixSparse<double> Mat;

  //! resisual force of the residual force
  std::vector<double> aRes0;

  //! maximum number of iteration for the linear solver
  unsigned int max_itr = 300;

  //! convergence tolerance for the linear solver
  double conv_tol = 1.0e-5;

  //! swich to activate preconditioer
  bool is_preconditioner;

  //! @brief preconditioiner
  std::vector<double> aDiaInv;

  //
  mutable std::vector<double> vec_tmp0, vec_tmp1;

  //! convergence history of the linear solver
  mutable std::vector<double> aConvHist;
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/deflap.cpp"
#endif

#endif /* def_h */
