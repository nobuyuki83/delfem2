/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_DEFARAP_H
#define DFM2_DEFARAP_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/ls_ilu_block_sparse.h"

// ---------------------------

namespace delfem2 {


/**
 * @discussion interface doesn't have CQuat, CMat3, CVec3, but the implementation has it.
 */
DFM2_INLINE void UpdateRotationsByMatchingCluster_Linear(
    std::vector<double> &aQuat1,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

DFM2_INLINE void UpdateRotationsByMatchingCluster_SVD(
    std::vector<double> &aQuat1,
    unsigned int ip,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

/**
 * @brief Edge-based As-Rigid-As-Possible shape deformation wihtout rotation
 */
class CDef_ArapEdgeLinearDisponly {
 public:
  CDef_ArapEdgeLinearDisponly(
      const std::vector<double> &aXYZ0,
      const std::vector<unsigned int> &aTri,
      double weight_bc0,
      std::vector<int> aBCFlag0);
  void Deform(
      std::vector<double> &aXYZ1,
      const std::vector<double> &aXYZ0);
  void MatVec(
      double *y,
      double alpha,
      const double *vec,
      double beta) const;
 private:
  void MakeLinearSystem(
      double *aRhs,
      const double *aXYZ0,
      const double *aXYZ1) const;
  void JacobiTVecTmp(
      double *y,
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
  CDef_ArapEdge() {}
  void Init(
      const std::vector<double> &aXYZ0,
      const std::vector<unsigned int> &aTri,
      double weight_bc0,
      const std::vector<int> &aBCFlag,
      bool is_preconditioner);
  void Deform(
      std::vector<double> &aXYZ1,
      std::vector<double> &aQuat,
      const std::vector<double> &aXYZ0);
  void MatVec(
      double *y,
      double alpha,
      const double *vec,
      double beta) const;
  void SolvePrecond(double *v) const;
 private:
  void JacobiTVecTmp(
      double *y,
      double alpha,
      double beta) const;
  void MakeLinearSystem(
      double *aRhs,
      const double *aXYZ0,
      const double *aXYZ1,
      const double *aQuat);
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

void UpdateQuaternions_Svd(
    std::vector<double> &vtx_quaternion,
    const std::vector<double> &vtx_xyz_ini,
    const std::vector<double> &vtx_xyz_def,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

/**
 * @brief As-Rigid-As-Possible shape deformation
 */
class Deformer_Arap {
 public:
  Deformer_Arap() {}
  void Init(
      const std::vector<double> &vtx_xyz_ini,
      const std::vector<unsigned int> &tri_vtx,
      bool is_preconditioner);
  void Deform(
      std::vector<double> &vtx_xyz_def,
      std::vector<double> &vtx_quaternion,
      const std::vector<double> &vtx_xyz_ini,
      const std::vector<int> &aBCFlag);
 public:
  mutable std::vector<double> convergence_history;
  std::vector<unsigned int> psup_ind, psup;
 private:
  bool is_preconditioner_; // use preconditioner or not
  std::vector<double> precomp_; // size: np * 9, precomputed component of sparse
  CMatrixSparse<double> sparse_;
  std::vector<double> residual_, update_;
  std::vector<double> tmp_vec0_, tmp_vec1_;
  std::vector<unsigned int> tmp_buffer_for_merge_;
  CPreconditionerILU<double> precond_;
};

/**
 * @brief As-Rigid-As-Possible shape deformation
 */
class Deformer_Arap2 {
 public:
  Deformer_Arap2() {}
  void Init(
      const std::vector<double> &vtx_xyz_ini,
      const std::vector<unsigned int> &tri_vtx,
      const std::vector<double> &vtx_quaternion,
      const std::vector<int> &dof_bcflag);
  void Deform(
      std::vector<double> &vtx_xyz_def,
      std::vector<double> &vtx_quaternion,
      const std::vector<double> &vtx_xyz_ini,
      const std::vector<int> &aBCFlag);
 public:
  mutable std::vector<double> convergence_history;
  std::vector<unsigned int> psup_ind, psup;
 private:
  std::vector<double> precomp_;
  CMatrixSparse<double> sparse_;
  std::vector<double> residual_, update_;
  std::vector<double> tmp_vec0_, tmp_vec1_;
  std::vector<unsigned int> tmp_buffer_for_merge_;
  CPreconditionerILU<double> precond_;
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/defarap.cpp"
#endif

#endif /* DFM2_DEFARAP_H */
