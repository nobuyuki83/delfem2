/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/deflap.h"

#include <cstring> // memcpy

#include "delfem2/geo3_v23m34q.h" // update rotation by matching cluster
#include "delfem2/lsvecx.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"

namespace delfem2::defarap {

void SetLinSys_LaplaceGraph_MeshTri3(
    CMatrixSparse<double> &mat_A) {
  mat_A.setZero();
  for (unsigned int ip = 0; ip < mat_A.nrowblk_; ++ip) {
    const auto dn = static_cast<double>(mat_A.col_ind_[ip + 1] - mat_A.col_ind_[ip]);
    for (unsigned int icrs = mat_A.col_ind_[ip]; icrs < mat_A.col_ind_[ip + 1]; ++icrs) {
      mat_A.val_crs_[icrs * 9 + 0 * 3 + 0] = -1.0;
      mat_A.val_crs_[icrs * 9 + 1 * 3 + 1] = -1.0;
      mat_A.val_crs_[icrs * 9 + 2 * 3 + 2] = -1.0;
    }
    mat_A.val_dia_[ip * 9 + 0 * 3 + 0] = dn;
    mat_A.val_dia_[ip * 9 + 1 * 3 + 1] = dn;
    mat_A.val_dia_[ip * 9 + 2 * 3 + 2] = dn;
  }
}

}

// ==================================================

void delfem2::CDef_LaplacianLinearAsym::Init(
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri) {
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri.data(),
      aTri.size() / 3, 3,
      aXYZ0.size() / 3);
  JArray_Sort(psup_ind, psup);
  mat_A.Initialize(
      static_cast<unsigned int>(aXYZ0.size() / 3), 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(), psup.size());
  // ---
  defarap::SetLinSys_LaplaceGraph_MeshTri3(mat_A);
  aRhs0.resize(aXYZ0.size());
  mat_A.MatVec(aRhs0.data(),
               1.0, aXYZ0.data(), 0.0);
}

void delfem2::CDef_LaplacianLinearAsym::Deform(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0,
    const std::vector<int> &aBCFlag) {    // ----------
  aRhs1 = aRhs0;
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) {
    if (aBCFlag[i] == 0) { continue; }
    aRhs1[i] = aXYZ1[i];
  }
  mat_A.SetFixedBC_Dia(aBCFlag.data(), 1.0);
  mat_A.SetFixedBC_Row(aBCFlag.data());
  aXYZ1 = aXYZ0;
  aHistConv = Solve_BiCGStab(
      aRhs1, aXYZ1,
      1.0e-5, 100, mat_A);
}


// ===================================================
// below: implementation of CDef_LaplacianLinearGram

void delfem2::CDef_LaplacianLinearGram::Init(
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri,
    bool is_preconditioner0) {
  this->is_preconditioner = is_preconditioner0;
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri.data(), aTri.size() / 3, 3,
      aXYZ0.size() / 3);
  JArray_Sort(psup_ind, psup);
  Mat.Initialize(aXYZ0.size() / 3, 3, true);
  Mat.SetPattern(
      psup_ind.data(), psup_ind.size(),
      psup.data(), psup.size());
  defarap::SetLinSys_LaplaceGraph_MeshTri3(Mat);
  const unsigned int np = Mat.nrowblk_;
  aRes0.assign(np * 3, 0.0);
  Mat.MatVec(
      aRes0.data(),
      -1.0, aXYZ0.data(), 0.0);
}

void delfem2::CDef_LaplacianLinearGram::SetBoundaryConditionToPreconditioner() {
  if (!is_preconditioner) { return; }
  // ---------
  // make jacobi preconditioner
  const unsigned int np = Mat.nrowblk_;
  aDiaInv.assign(np * 9, 0.0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int icrs = Mat.col_ind_[ip]; icrs < Mat.col_ind_[ip + 1]; ++icrs) {
      unsigned int jp0 = Mat.row_ptr_[icrs];
      MatTMat3_ScaleAdd(
          aDiaInv.data() + jp0 * 9,
          Mat.val_crs_.data() + icrs * 9,
          Mat.val_crs_.data() + icrs * 9,
          1.0, 1.0); // del. prev. value and set new vaue
    }
    {
      MatTMat3_ScaleAdd(
          aDiaInv.data() + ip * 9,
          Mat.val_dia_.data() + ip * 9,
          Mat.val_dia_.data() + ip * 9,
          1.0, 1.0); // del. prev. value and set new vaue
    }
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (int i = 0; i < 3; ++i) {
      if (aBCFlag[ip * 3 + i] == 0) { continue; }
      aDiaInv[ip * 9 + i * 3 + i] += weight_bc;
    }
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    Inverse_Mat3(aDiaInv.data() + ip * 9);
  }
}

void delfem2::CDef_LaplacianLinearGram::Deform(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0) const {
  vec_tmp0.resize(aXYZ0.size());
  vec_tmp1.resize(aXYZ0.size());
  vec_tmp2.resize(aXYZ0.size());
  //
  std::vector<double> &aRhs = vec_tmp0;
  { // making RHS vector for elastic deformation
    vec_tmp2 = aRes0;
    Mat.MatVec(vec_tmp2.data(),
               +1.0, aXYZ1.data(), 1.0);
    Mat.MatTVec(aRhs.data(),
                -1.0, vec_tmp2.data(), 0.0);
  }
  std::vector<double> &aUpd = vec_tmp1;
  aUpd.assign(aXYZ0.size(), 0.0);
  if (is_preconditioner) {
    const std::size_t n = aRhs.size();
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_PCG(
        CVecXd(aRhs), CVecXd(aUpd), CVecXd(tmp0), CVecXd(tmp1),
        1.0e-7, 300, *this, *this);
  } else {
    const std::size_t n = aRhs.size();
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_CG(
        CVecXd(aRhs), CVecXd(aUpd), CVecXd(tmp0), CVecXd(tmp1),
        1.0e-7, 300, *this);
  }
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) { aXYZ1[i] += aUpd[i]; }
}

void delfem2::CDef_LaplacianLinearGram::MatVec
    (double *y,
     double alpha, const double *vec, double beta) const {
  Mat.MatVec(vec_tmp2.data(),
             1, vec, 0.0);
  Mat.MatTVec(y,
              alpha, vec_tmp2.data(), beta);
  // add diagonal for fixed boundary condition
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) {
    if (aBCFlag[i] == 0) { continue; }
    y[i] += weight_bc * vec[i];
  }
}

// for preconditioner
void delfem2::CDef_LaplacianLinearGram::SolvePrecond(double *v) const {
  const auto np = static_cast<unsigned int>(aBCFlag.size() / 3);
  for (unsigned int ip = 0; ip < np; ++ip) {
    double tmp[3];
    MatVec3(tmp, aDiaInv.data() + ip * 9, v + ip * 3);
    v[ip * 3 + 0] = tmp[0];
    v[ip * 3 + 1] = tmp[1];
    v[ip * 3 + 2] = tmp[2];
  }
}

// above: delfem2::CDef_LaplacianLinearGram
// =========================================================================
// below: delfem2::CDef_LaplacianLinear


namespace delfem2::defarap {

void DualLaplacianSymbolic_3x3(
    std::vector<double> &eM,
    const std::vector<unsigned int> &aIP) {
  const auto nIP = static_cast<unsigned int>(aIP.size());
  const unsigned int nNg = nIP - 1; // number of neighbor
  double dn = (double) nNg;
  eM.assign(nIP * nIP * 9, 0.0);
  const CMat3d L1 = CMat3d::Identity();
  L1.AddToScale(eM.data() + (nNg * nIP + nNg) * 9, +dn * dn);
  for (unsigned int jjp = 0; jjp < nNg; ++jjp) {
    L1.AddToScale(eM.data() + (nNg * nIP + jjp) * 9, -dn);
    L1.AddToScale(eM.data() + (jjp * nIP + nNg) * 9, -dn);
    for (unsigned int kkp = 0; kkp < nNg; ++kkp) {
      L1.AddToScale(eM.data() + (jjp * nIP + kkp) * 9, +1.0);
    }
  }
}

}

void delfem2::CDef_LaplacianLinear::Init(
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri,
    bool is_preconditioner_) {
  const auto np = static_cast<unsigned int>(aXYZ0.size() / 3);
  this->is_preconditioner = is_preconditioner_;
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri.data(), aTri.size() / 3, 3,
      np);
  JArray_Sort(psup_ind, psup);
  {
    std::vector<unsigned int> psup_ind1, psup1;
    JArray_Extend(
        psup_ind1, psup1,
        psup_ind.data(), psup_ind.size(), psup.data());
    JArray_Sort(psup_ind1, psup1);
    Mat.Initialize(np, 3, true);
    assert(psup_ind1.size() == np + 1);
    Mat.SetPattern(
        psup_ind1.data(), psup_ind1.size(),
        psup1.data(), psup1.size());
  }

  std::vector<unsigned int> tmp_buffer;
  Mat.setZero();
  for (unsigned int ip = 0; ip < np; ++ip) {
    std::vector<unsigned int> aIP;
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      aIP.push_back(psup[ipsup]);
    }
    aIP.push_back(ip);
    std::vector<double> eM;
    delfem2::defarap::DualLaplacianSymbolic_3x3(eM, aIP);
    Mearge(Mat,
           aIP.size(), aIP.data(),
           aIP.size(), aIP.data(),
           9, eM.data(),
           tmp_buffer);
  }

  aRes0.resize(aXYZ0.size());
  Mat.MatVec(
      aRes0.data(),
      -1.0, aXYZ0.data(), 0.0);

  aBCFlag.assign(aXYZ0.size(), 0);

  this->Prec.Clear();
  if (is_preconditioner) {
    this->Prec.Initialize_ILUk(Mat, 20);
  }
}

void delfem2::CDef_LaplacianLinear::SetValueToPreconditioner() {
  if (!is_preconditioner) { return; }
  //
  const unsigned int np = Mat.nrowblk_;
  assert(aBCFlag.size() == np * 3);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (int idim = 0; idim < 3; ++idim) {
      if (aBCFlag[ip * 3 + idim] == 0) { continue; }
      Mat.val_dia_[ip * 9 + idim * 3 + idim] += weight_bc;
    }
  }
  for (auto &iip: aIdpNrm) {
    const unsigned int ip0 = iip.first;
    const double *n0 = iip.second.p;
    Mat.val_dia_[ip0 * 9 + 0 * 3 + 0] += weight_nrm * n0[0] * n0[0];
    Mat.val_dia_[ip0 * 9 + 0 * 3 + 1] += weight_nrm * n0[0] * n0[1];
    Mat.val_dia_[ip0 * 9 + 0 * 3 + 2] += weight_nrm * n0[0] * n0[2];
    Mat.val_dia_[ip0 * 9 + 1 * 3 + 0] += weight_nrm * n0[1] * n0[0];
    Mat.val_dia_[ip0 * 9 + 1 * 3 + 1] += weight_nrm * n0[1] * n0[1];
    Mat.val_dia_[ip0 * 9 + 1 * 3 + 2] += weight_nrm * n0[1] * n0[2];
    Mat.val_dia_[ip0 * 9 + 2 * 3 + 0] += weight_nrm * n0[2] * n0[0];
    Mat.val_dia_[ip0 * 9 + 2 * 3 + 1] += weight_nrm * n0[2] * n0[1];
    Mat.val_dia_[ip0 * 9 + 2 * 3 + 2] += weight_nrm * n0[2] * n0[2];
  }
  // --------
  this->Prec.CopyValue(Mat);
  this->Prec.Decompose();
  /*
  for(int ip=0;ip<np;++ip){
    for(unsigned int icrs=Prec.mat.colInd[ip];icrs<Prec.mat.colInd[ip+1];++icrs){
      unsigned int jp = Prec.mat.rowPtr[icrs];
      const double* p = Prec.mat.valCrs.data()+icrs*9;
      std::cout << ip << " " << jp << "  --> " << p[0] << " " << p[4] << " " << p[8] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
    }
  }
   */
  // -------
  for (auto &iip: aIdpNrm) {
    const unsigned int ip0 = iip.first;
    const double *n0 = iip.second.p;
    Mat.val_dia_[ip0 * 9 + 0 * 3 + 0] -= weight_nrm * n0[0] * n0[0];
    Mat.val_dia_[ip0 * 9 + 0 * 3 + 1] -= weight_nrm * n0[0] * n0[1];
    Mat.val_dia_[ip0 * 9 + 0 * 3 + 2] -= weight_nrm * n0[0] * n0[2];
    Mat.val_dia_[ip0 * 9 + 1 * 3 + 0] -= weight_nrm * n0[1] * n0[0];
    Mat.val_dia_[ip0 * 9 + 1 * 3 + 1] -= weight_nrm * n0[1] * n0[1];
    Mat.val_dia_[ip0 * 9 + 1 * 3 + 2] -= weight_nrm * n0[1] * n0[2];
    Mat.val_dia_[ip0 * 9 + 2 * 3 + 0] -= weight_nrm * n0[2] * n0[0];
    Mat.val_dia_[ip0 * 9 + 2 * 3 + 1] -= weight_nrm * n0[2] * n0[1];
    Mat.val_dia_[ip0 * 9 + 2 * 3 + 2] -= weight_nrm * n0[2] * n0[2];
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (int idim = 0; idim < 3; ++idim) {
      if (aBCFlag[ip * 3 + idim] == 0) { continue; }
      Mat.val_dia_[ip * 9 + idim * 3 + idim] -= weight_bc;
    }
  }
}

void delfem2::CDef_LaplacianLinear::Deform(
    std::vector<double> &vtx_xyz_def,
    const std::vector<double> &vtx_xyz_ini) const {
  vec_tmp0.resize(vtx_xyz_ini.size());
  vec_tmp1.resize(vtx_xyz_ini.size());
  // ------------------------------
  std::vector<double> &aRhs = vec_tmp0;
  std::memcpy(aRhs.data(), aRes0.data(), aRes0.size() * sizeof(double));
  Mat.MatVec(aRhs.data(),
             -1.0, vtx_xyz_def.data(), -1.0);
  std::vector<double> &aUpd = vec_tmp1;
  aUpd.assign(vtx_xyz_ini.size(), 0.0);
  if (is_preconditioner) {
    std::size_t n = aRhs.size();
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_PCG(
        CVecXd(aRhs),
        CVecXd(aUpd),
        CVecXd(tmp0),
        CVecXd(tmp1),
        this->conv_tol, this->max_itr, *this, Prec);
  } else {
    std::size_t n = aRhs.size();
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_CG(
        CVecXd(aRhs),
        CVecXd(aUpd),
        CVecXd(tmp0),
        CVecXd(tmp1),
        this->conv_tol, this->max_itr, *this);
  }
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) { vtx_xyz_def[i] += aUpd[i]; }
}

void delfem2::CDef_LaplacianLinear::MatVec(
    double *y,
    double alpha, const double *vec, double beta) const {
  Mat.MatVec(y,
             alpha, vec, beta);
  //
  for (const auto &iip: aIdpNrm) {
    const unsigned int ip0 = iip.first;
    const double *n0 = iip.second.p;
    const double d = Dot3(n0, vec + ip0 * 3);
    y[ip0 * 3 + 0] += weight_nrm * alpha * d * n0[0];
    y[ip0 * 3 + 1] += weight_nrm * alpha * d * n0[1];
    y[ip0 * 3 + 2] += weight_nrm * alpha * d * n0[2];
  }
  // add diagonal for fixed boundary condition
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) {
    if (aBCFlag[i] == 0) { continue; }
    y[i] += weight_bc * vec[i];
  }
}

// above: delfem2::CDef_LaplacianLinear
// ============================================
// below: delfem2::CDef_LaplacianLinearLightweight


namespace delfem2::defarap {

void DualLaplacianSymbolic_1x1(
    std::vector<double> &eM,
    const std::vector<unsigned int> &aIP) {
  const auto nIP = static_cast<unsigned int>(aIP.size());
  assert(nIP > 0);
  const unsigned int nNg = nIP - 1; // number of neighbor
  const auto dn = static_cast<double>(nNg);
  eM.assign(nIP * nIP, 0.0);
  eM[nNg * nIP + nNg] = dn * dn;
  for (unsigned int jjp = 0; jjp < nNg; ++jjp) {
    eM[nNg * nIP + jjp] = -dn;
    eM[jjp * nIP + nNg] = -dn;
    for (unsigned int kkp = 0; kkp < nNg; ++kkp) {
      eM[jjp * nIP + kkp] = 1.0;
    }
  }
}

}

void delfem2::CDef_LaplacianLinearDegenerate::Init(
    const std::vector<double> &vtx_xyz_ini,
    const std::vector<unsigned int> &tri_vtx,
    bool is_preconditioner_) {
  this->is_preconditioner = is_preconditioner_;
  const auto np = static_cast<unsigned int>(vtx_xyz_ini.size() / 3);
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(
      psup_ind, psup,
      tri_vtx.data(), tri_vtx.size() / 3, 3,
      np);
  JArray_Sort(psup_ind, psup);
  {
    std::vector<unsigned int> psup_ind1, psup1;
    JArray_Extend(psup_ind1, psup1,
                  psup_ind.data(), psup_ind.size(), psup.data());
    JArray_Sort(psup_ind1, psup1);
    Mat.Initialize(np, 1, true);
    assert(psup_ind1.size() == np + 1);
    Mat.SetPattern(psup_ind1.data(), psup_ind1.size(),
                   psup1.data(), psup1.size());
  }

  std::vector<unsigned int> tmp_buffer;
  Mat.setZero();
  for (unsigned int ip = 0; ip < np; ++ip) {
    std::vector<unsigned int> aIP;
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      aIP.push_back(psup[ipsup]);
    }
    aIP.push_back(ip);
    std::vector<double> eM;
    delfem2::defarap::DualLaplacianSymbolic_1x1(eM, aIP);
    Mearge(
        Mat,
        aIP.size(), aIP.data(),
        aIP.size(), aIP.data(),
        1, eM.data(),
        tmp_buffer);
  }

  aRes0.resize(vtx_xyz_ini.size());
  Mat.MatVecDegenerate(
      aRes0.data(),
      3, -1.0, vtx_xyz_ini.data(), 0.0);
  aBCFlag.assign(vtx_xyz_ini.size(), 0);

  this->aDiaInv.clear();
  if (is_preconditioner) {
    aDiaInv.resize(np * 9, 0);
    for (unsigned int ip = 0; ip < np; ++ip) {
      aDiaInv[ip * 9 + 0] = 1.0;
      aDiaInv[ip * 9 + 4] = 1.0;
      aDiaInv[ip * 9 + 8] = 1.0;
    }
  }
}

void delfem2::CDef_LaplacianLinearDegenerate::SetBoundaryConditionToPreconditioner() {
  if (!is_preconditioner) { return; }
  // ---------
  // make jacobi preconditioner
  const unsigned int np = Mat.nrowblk_;
  aDiaInv.assign(np * 9, 0.0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    double v0 = Mat.val_dia_[ip];
    aDiaInv[ip * 9 + 0] = v0;
    aDiaInv[ip * 9 + 4] = v0;
    aDiaInv[ip * 9 + 8] = v0;
  }
  for (auto &iip: aIdpNrm) {
    const unsigned int ip0 = iip.first;
    const double *n0 = iip.second.p;
    aDiaInv[ip0 * 9 + 0 * 3 + 0] += weight_nrm * n0[0] * n0[0];
    aDiaInv[ip0 * 9 + 0 * 3 + 1] += weight_nrm * n0[0] * n0[1];
    aDiaInv[ip0 * 9 + 0 * 3 + 2] += weight_nrm * n0[0] * n0[2];
    aDiaInv[ip0 * 9 + 1 * 3 + 0] += weight_nrm * n0[1] * n0[0];
    aDiaInv[ip0 * 9 + 1 * 3 + 1] += weight_nrm * n0[1] * n0[1];
    aDiaInv[ip0 * 9 + 1 * 3 + 2] += weight_nrm * n0[1] * n0[2];
    aDiaInv[ip0 * 9 + 2 * 3 + 0] += weight_nrm * n0[2] * n0[0];
    aDiaInv[ip0 * 9 + 2 * 3 + 1] += weight_nrm * n0[2] * n0[1];
    aDiaInv[ip0 * 9 + 2 * 3 + 2] += weight_nrm * n0[2] * n0[2];
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (int i = 0; i < 3; ++i) {
      if (aBCFlag[ip * 3 + i] == 0) { continue; }
      aDiaInv[ip * 9 + i * 3 + i] += weight_bc;
    }
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    Inverse_Mat3(aDiaInv.data() + ip * 9);
  }
}

void delfem2::CDef_LaplacianLinearDegenerate::Deform(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0) const {
  vec_tmp0.resize(aXYZ0.size());
  vec_tmp1.resize(aXYZ0.size());
  // ------------------------------
  std::vector<double> &aRhs = vec_tmp0;
  std::memcpy(aRhs.data(), aRes0.data(), aRes0.size() * sizeof(double));
  Mat.MatVecDegenerate(
      aRhs.data(),
      3, -1.0, aXYZ1.data(), -1.0);
  std::vector<double> &aUpd = vec_tmp1;
  aUpd.assign(aXYZ0.size(), 0.0);
  if (is_preconditioner) {
    const std::size_t n = aRhs.size();
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_PCG(
        CVecXd(aRhs),
        CVecXd(aUpd),
        CVecXd(tmp0),
        CVecXd(tmp1),
        this->conv_tol, this->max_itr, *this, *this);
  } else {
    const std::size_t n = aRhs.size();
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_CG(
        CVecXd(aRhs),
        CVecXd(aUpd),
        CVecXd(tmp0),
        CVecXd(tmp1),
        this->conv_tol, this->max_itr, *this);
  }
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) { aXYZ1[i] += aUpd[i]; }
}

void delfem2::CDef_LaplacianLinearDegenerate::MatVec(
    double *y,
    double alpha, const double *vec, double beta) const {
  Mat.MatVecDegenerate(y,
                       3, alpha, vec, beta);
  // sliding boundary condition
  for (const auto &iip: aIdpNrm) {
    const unsigned int ip0 = iip.first;
    const double *n0 = iip.second.p;
    const double d = Dot3(n0, vec + ip0 * 3);
    y[ip0 * 3 + 0] += weight_nrm * alpha * d * n0[0];
    y[ip0 * 3 + 1] += weight_nrm * alpha * d * n0[1];
    y[ip0 * 3 + 2] += weight_nrm * alpha * d * n0[2];
  }
  // add diagonal for fixed boundary condition
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) {
    if (aBCFlag[i] == 0) { continue; }
    y[i] += weight_bc * vec[i];
  }
}

void delfem2::CDef_LaplacianLinearDegenerate::SolvePrecond(double *v) const {
  const auto np = static_cast<unsigned int>(aBCFlag.size() / 3);
  for (unsigned int ip = 0; ip < np; ++ip) {
    double tmp[3];
    MatVec3(tmp, aDiaInv.data() + ip * 9, v + ip * 3);
    v[ip * 3 + 0] = tmp[0];
    v[ip * 3 + 1] = tmp[1];
    v[ip * 3 + 2] = tmp[2];
  }
}

