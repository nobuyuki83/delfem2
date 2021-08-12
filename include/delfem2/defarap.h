/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_DEFARAP_H
#define DFM2_DEFARAP_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/lsilu_mats.h"

// ---------------------------

namespace delfem2 {

/**
 * @brief Edge-based As-Rigid-As-Possible shape deformation wihtout rotation
 */
class CDef_ArapEdgeLinearDisponly {
public:
  CDef_ArapEdgeLinearDisponly(
      const std::vector<double>& aXYZ0,
      const std::vector<unsigned int>& aTri,
      double weight_bc0,
      std::vector<int>  aBCFlag0);
  void Deform(
      std::vector<double>& aXYZ1,
      const std::vector<double>& aXYZ0);
  void MatVec(
      double* y,
      double alpha,
      const double* vec,
      double beta) const;
private:
  void MakeLinearSystem(
      double* aRhs,
      const double* aXYZ0,
      const double* aXYZ1) const;
  void JacobiTVecTmp(
      double*y ,
      double alpha,
      double beta) const;
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
 * @brief Edge-based As-Rigid-As-Possible shape deformation
 */
class CDef_ArapEdge {
public:
  CDef_ArapEdge(){}
  void Init(
      const std::vector<double>& aXYZ0,
      const std::vector<unsigned int>& aTri,
      double weight_bc0,
      const std::vector<int>& aBCFlag,
      bool is_preconditioner);
  void Deform(
      std::vector<double>& aXYZ1,
      std::vector<double>& aQuat,
      const std::vector<double>& aXYZ0);
  void MatVec(
      double* y,
      double alpha,
      const double* vec,
      double beta) const;
  void SolvePrecond(double* v) const;
private:
  void JacobiTVecTmp(
      double*y,
      double alpha,
      double beta) const;
  void MakeLinearSystem(
      double* aRhs,
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
  void Init(
      const std::vector<double>& aXYZ0,
      const std::vector<unsigned int>& aTri,
      bool is_preconditioner);
  void Deform(
      std::vector<double>& aXYZ1,
      std::vector<double>& aQuat,
      const std::vector<double>& aXYZ0,
      const std::vector<int>& aBCFlag);
  void UpdateQuats_SVD(
      std::vector<double>& aXYZ1,
      std::vector<double>& aQuat1,
      const std::vector<double>& aXYZ0) const;
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

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/defarap.cpp"
#endif

#endif /* def_h */
