/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_DEF_H
#define DFM2_DEF_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/ilu_mats.h"

// ---------------------------

namespace delfem2 {

class CDef_LaplacianLinearAsym
{
public:
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri);
  void Deform(std::vector<double>& aXYZ1,
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
  void SetBoundaryCondition(const std::vector<int>& aBCFlag);
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
 * @brief deformation classs where the matrix is computed as L^TL
 */
class CDef_LaplacianLinear{
public:
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri,
            bool is_preconditioner);
  void Deform(std::vector<double>& aXYZ1,
              const std::vector<double>& aXYZ0) const;
  void SetBoundaryCondition(const std::vector<int>& aBCFlag);
  // -----------
  void MatVec(double* y,
              double alpha, const double* vec,  double beta) const;
public:
  CMatrixSparse<double> Mat;
  double weight_bc = 100;
  std::vector<double> aRes0;
  std::vector<int> aBCFlag;
  //
  bool is_preconditioner;
  CPreconditionerILU<double> Prec;
  // 
  mutable std::vector<double> vec_tmp0, vec_tmp1;
  mutable std::vector<double> aConvHist;
};

// =====================================

class CDef_ArapEdgeLinearDisponly {
public:
  CDef_ArapEdgeLinearDisponly(const std::vector<double>& aXYZ0,
                          const std::vector<unsigned int>& aTri,
                          double weight_bc0,
                          const std::vector<int>& aBCFlag0);
  void Deform(std::vector<double>& aXYZ1,
              const std::vector<double>& aXYZ0);
  void MatVec(double* y,
              double alpha, const double* vec,  double beta) const;
private:
  void MakeLinearSystem(double* aRhs,
                        const double* aXYZ0,
                        const double* aXYZ1) const;
  void JacobiTVecTmp(double*y ,
                     double alpha, double beta) const;
public:
  std::vector<unsigned int> psup_ind;
  std::vector<unsigned int> psup;
  const double weight_bc;
  const std::vector<int> aBCFlag;
  // -------------
  std::vector<double> aMatEdge;
  mutable std::vector<double> vec_tmp;
};


/**
 * @brief Edge based As-Rigid-As-Possible shape deformation
 */
class CDef_ArapEdge {
public:
  CDef_ArapEdge(){}
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri,
            double weight_bc0,
            const std::vector<int>& aBCFlag,
            bool is_preconditioner);
  void Deform(std::vector<double>& aXYZ1,
              std::vector<double>& aQuat,
              const std::vector<double>& aXYZ0);
  void MatVec(double* y,
              double alpha, const double* vec,  double beta) const;
  void SolvePrecond(double* v) const;
private:
  void JacobiTVecTmp(double*y,
                     double alpha, double beta) const;
  void MakeLinearSystem(double* aRhs,
                        const double* aXYZ0,
                        const double* aXYZ1,
                        const double* aQuat);
  void MakePreconditionerJacobi();
public:
  std::vector<unsigned int> psup_ind, psup;
  double weight_bc;
  std::vector<int> aBCFlag;
  bool is_preconditioner;
  // -------------
  std::vector<double> aMatEdge;
  std::vector<double> aDiaInv; // for jacobi preconditining
  mutable std::vector<double> vec_tmp;
};

// =============================================

/**
 * @brief As-Rigid-As-Possible shape deformation
 */
class CDef_Arap {
public:
  CDef_Arap(){}
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri,
            bool is_preconditioner);
  void Deform(std::vector<double>& aXYZ1,
              std::vector<double>& aQuat,
              const std::vector<double>& aXYZ0,
              const std::vector<int>& aBCFlag);
public:
  mutable std::vector<double> aConvHist;
  std::vector<unsigned int> psup_ind, psup;
private:
  bool is_preconditioner;
  std::vector<double> Precomp;
  CMatrixSparse<double> Mat;
  std::vector<double> aRes1, aUpd1;
  CPreconditionerILU<double> Prec;
};

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/def.cpp"
#endif

#endif /* def_h */
