/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/lsmats.h"

#include <cassert>
#include <vector>
#include <climits>
#include <complex>

namespace delfem2 {
namespace mats {

DFM2_INLINE double MatNorm_Assym(
    const double *V0,
    unsigned int n0,
    unsigned int m0,
    const double *V1) {
  double s = 0.0;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < m0; ++j) {
      double v0 = V0[i * m0 + j];
      double v1 = V1[j * n0 + i];
      s += (v0 - v1) * (v0 - v1);
    }
  }
  return s;
}

DFM2_INLINE double MatNorm(
    const double *V,
    unsigned int n,
    unsigned int m) {
  double s = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < m; ++j) {
      double v = V[i * m + j];
      s += v * v;
    }
  }
  return s;
}

DFM2_INLINE double MatNorm_Assym(
    const double *V,
    unsigned int n) {
  double s = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      double v0 = V[i * n + j];
      double v1 = V[j * n + i];
      s += (v0 - v1) * (v0 - v1);
    }
  }
  return s;
}

template<typename T>
void MatVec_MatSparseCRS_Blk11(
    T *y,
    T alpha,
    unsigned nrowblk,
    const T *vcrs,
    const T *vdia,
    const unsigned int *colind,
    const unsigned int *rowptr,
    const T *x) {
  for (unsigned int iblk = 0; iblk < nrowblk; iblk++) {
    const unsigned int colind0 = colind[iblk];
    const unsigned int colind1 = colind[iblk + 1];
    for (unsigned int icrs = colind0; icrs < colind1; icrs++) {
      const unsigned int jblk0 = rowptr[icrs];
      y[iblk] += alpha * vcrs[icrs] * x[jblk0];
    }
    y[iblk] += alpha * vdia[iblk] * x[iblk];
  }
}

template<typename T>
void MatVec_MatSparseCRS_Blk22(
    T *y,
    T alpha,
    unsigned nrowblk,
    const T *vcrs,
    const T *vdia,
    const unsigned int *colind,
    const unsigned int *rowptr,
    const T *x) {
  for (unsigned int iblk = 0; iblk < nrowblk; iblk++) {
    const unsigned int icrs0 = colind[iblk];
    const unsigned int icrs1 = colind[iblk + 1];
    for (unsigned int icrs = icrs0; icrs < icrs1; icrs++) {
      const unsigned int jblk0 = rowptr[icrs];
      y[iblk * 2 + 0] += alpha * (vcrs[icrs * 4] * x[jblk0 * 2 + 0] + vcrs[icrs * 4 + 1] * x[jblk0 * 2 + 1]);
      y[iblk * 2 + 1] += alpha * (vcrs[icrs * 4 + 2] * x[jblk0 * 2 + 0] + vcrs[icrs * 4 + 3] * x[jblk0 * 2 + 1]);
    }
    y[iblk * 2 + 0] += alpha * (vdia[iblk * 4 + 0] * x[iblk * 2 + 0] + vdia[iblk * 4 + 1] * x[iblk * 2 + 1]);
    y[iblk * 2 + 1] += alpha * (vdia[iblk * 4 + 2] * x[iblk * 2 + 0] + vdia[iblk * 4 + 3] * x[iblk * 2 + 1]);
  }
}

template<typename T>
void MatVec_MatSparseCRS_Blk33(
    T *y,
    T alpha,
    unsigned nrowblk,
    const T *vcrs,
    const T *vdia,
    const unsigned int *colind,
    const unsigned int *rowptr,
    const T *x) {
  for (unsigned int iblk = 0; iblk < nrowblk; iblk++) {
    const unsigned int icrs0 = colind[iblk];
    const unsigned int icrs1 = colind[iblk + 1];
    for (unsigned int icrs = icrs0; icrs < icrs1; icrs++) {
      const unsigned int jblk0 = rowptr[icrs];
      const unsigned int i0 = iblk * 3;
      const unsigned int j0 = jblk0 * 3;
      const unsigned int k0 = icrs * 9;
      y[i0 + 0] += alpha * (vcrs[k0 + 0] * x[j0 + 0] + vcrs[k0 + 1] * x[j0 + 1] + vcrs[k0 + 2] * x[j0 + 2]);
      y[i0 + 1] += alpha * (vcrs[k0 + 3] * x[j0 + 0] + vcrs[k0 + 4] * x[j0 + 1] + vcrs[k0 + 5] * x[j0 + 2]);
      y[i0 + 2] += alpha * (vcrs[k0 + 6] * x[j0 + 0] + vcrs[k0 + 7] * x[j0 + 1] + vcrs[k0 + 8] * x[j0 + 2]);
    }
    {
      const unsigned int i0 = iblk * 3;
      const unsigned int k0 = iblk * 9;
      y[i0 + 0] += alpha * (vdia[k0 + 0] * x[i0 + 0] + vdia[k0 + 1] * x[i0 + 1] + vdia[k0 + 2] * x[i0 + 2]);
      y[i0 + 1] += alpha * (vdia[k0 + 3] * x[i0 + 0] + vdia[k0 + 4] * x[i0 + 1] + vdia[k0 + 5] * x[i0 + 2]);
      y[i0 + 2] += alpha * (vdia[k0 + 6] * x[i0 + 0] + vdia[k0 + 7] * x[i0 + 1] + vdia[k0 + 8] * x[i0 + 2]);
    }
  }
}

template<typename T>
void MatVec_MatSparseCRS_Blk44(
    T *y,
    T alpha,
    unsigned nrowblk,
    const T *vcrs,
    const T *vdia,
    const unsigned int *colind,
    const unsigned int *rowptr,
    const T *x) {
  for (unsigned int iblk = 0; iblk < nrowblk; iblk++) {
    const unsigned int icrs0 = colind[iblk];
    const unsigned int icrs1 = colind[iblk + 1];
    for (unsigned int icrs = icrs0; icrs < icrs1; icrs++) {
      const unsigned int jblk0 = rowptr[icrs];
      const unsigned int i0 = iblk * 4;
      const unsigned int j0 = jblk0 * 4;
      const unsigned int k0 = icrs * 16;
      y[i0 + 0] +=
          alpha * (
              vcrs[k0 + 0] * x[j0 + 0] +
                  vcrs[k0 + 1] * x[j0 + 1] +
                  vcrs[k0 + 2] * x[j0 + 2] +
                  vcrs[k0 + 3] * x[j0 + 3]);
      y[i0 + 1] +=
          alpha * (
              vcrs[k0 + 4] * x[j0 + 0] +
                  vcrs[k0 + 5] * x[j0 + 1] +
                  vcrs[k0 + 6] * x[j0 + 2] +
                  vcrs[k0 + 7] * x[j0 + 3]);
      y[i0 + 2] +=
          alpha * (
              vcrs[k0 + 8] * x[j0 + 0] +
                  vcrs[k0 + 9] * x[j0 + 1] +
                  vcrs[k0 + 10] * x[j0 + 2] +
                  vcrs[k0 + 11] * x[j0 + 3]);
      y[i0 + 3] +=
          alpha * (
              vcrs[k0 + 12] * x[j0 + 0] +
                  vcrs[k0 + 13] * x[j0 + 1] +
                  vcrs[k0 + 14] * x[j0 + 2] +
                  vcrs[k0 + 15] * x[j0 + 3]);
    }
    {
      const unsigned int i0 = iblk * 4;
      const unsigned int k0 = iblk * 16;
      y[i0 + 0] +=
          alpha * (
              vdia[k0 + 0] * x[i0 + 0] +
                  vdia[k0 + 1] * x[i0 + 1] +
                  vdia[k0 + 2] * x[i0 + 2] +
                  vdia[k0 + 3] * x[i0 + 3]);
      y[i0 + 1] +=
          alpha * (
              vdia[k0 + 4] * x[i0 + 0] +
                  vdia[k0 + 5] * x[i0 + 1] +
                  vdia[k0 + 6] * x[i0 + 2] +
                  vdia[k0 + 7] * x[i0 + 3]);
      y[i0 + 2] +=
          alpha * (
              vdia[k0 + 8] * x[i0 + 0] +
                  vdia[k0 + 9] * x[i0 + 1] +
                  vdia[k0 + 10] * x[i0 + 2] +
                  vdia[k0 + 11] * x[i0 + 3]);
      y[i0 + 3] +=
          alpha * (
              vdia[k0 + 12] * x[i0 + 0] +
                  vdia[k0 + 13] * x[i0 + 1] +
                  vdia[k0 + 14] * x[i0 + 2] +
                  vdia[k0 + 15] * x[i0 + 3]);
    }
  }
}

}
}

// -------------------------------------------------------

// Calc Matrix Vector Product
// {y} = alpha*[A]{x} + beta*{y}
template<typename T>
void delfem2::CMatrixSparse<T>::MatVec(
    T *y,
    T alpha,
    const T *x,
    T beta) const {
  const unsigned int ndofcol = nrowdim_ * nrowblk_;
  for (unsigned int i = 0; i < ndofcol; ++i) { y[i] *= beta; }
  // --------
  if (nrowdim_ == 1 && ncoldim_ == 1) {
    mats::MatVec_MatSparseCRS_Blk11(
        y,
        alpha, nrowblk_, val_crs_.data(), val_dia_.data(),
        col_ind_.data(), row_ptr_.data(), x);
  } else if (nrowdim_ == 2 && ncoldim_ == 2) {
    mats::MatVec_MatSparseCRS_Blk22(
        y,
        alpha, nrowblk_, val_crs_.data(), val_dia_.data(),
        col_ind_.data(), row_ptr_.data(), x);
  } else if (nrowdim_ == 3 && ncoldim_ == 3) {
    mats::MatVec_MatSparseCRS_Blk33(
        y,
        alpha, nrowblk_, val_crs_.data(), val_dia_.data(),
        col_ind_.data(), row_ptr_.data(), x);
  } else if (nrowdim_ == 4 && ncoldim_ == 4) {
    mats::MatVec_MatSparseCRS_Blk44(
        y,
        alpha, nrowblk_, val_crs_.data(), val_dia_.data(),
        col_ind_.data(), row_ptr_.data(), x);
  } else {
    const unsigned int blksize = nrowdim_ * ncoldim_;
    const T *vcrs = val_crs_.data();
    const T *vdia = val_dia_.data();
    const unsigned int *colind = col_ind_.data();
    const unsigned int *rowptr = row_ptr_.data();
    //
    for (unsigned int iblk = 0; iblk < nrowblk_; iblk++) {
      const unsigned int colind0 = colind[iblk];
      const unsigned int colind1 = colind[iblk + 1];
      for (unsigned int icrs = colind0; icrs < colind1; icrs++) {
        assert(icrs < row_ptr_.size());
        const unsigned int jblk0 = rowptr[icrs];
        assert(jblk0 < ncolblk_);
        for (unsigned int idof = 0; idof < nrowdim_; idof++) {
          for (unsigned int jdof = 0; jdof < ncoldim_; jdof++) {
            y[iblk * nrowdim_ + idof] +=
                alpha * vcrs[icrs * blksize + idof * ncoldim_ + jdof] * x[jblk0 * ncoldim_ + jdof];
          }
        }
      }
      for (unsigned int idof = 0; idof < nrowdim_; idof++) {
        for (unsigned int jdof = 0; jdof < ncoldim_; jdof++) {
          y[iblk * nrowdim_ + idof] +=
              alpha * vdia[iblk * blksize + idof * ncoldim_ + jdof] * x[iblk * ncoldim_ + jdof];
        }
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMatrixSparse<float>::MatVec(
    float *y, float alpha, const float *x, float beta) const;
template void delfem2::CMatrixSparse<double>::MatVec(
    double *y, double alpha, const double *x, double beta) const;
template void delfem2::CMatrixSparse<std::complex<double>>::MatVec(
    std::complex<double> *y, std::complex<double> alpha,
    const std::complex<double> *x, std::complex<double> beta) const;
#endif


// -------------------------------------------------------

/**
 * @brief Calc Matrix Vector Product {y} = alpha*[A]{x} + beta*{y}
 * the 1x1 sparse matrix is expanded as the len x len sparse matrix
 */

template<typename T>
void delfem2::CMatrixSparse<T>::MatVecDegenerate(
    T *y,
    unsigned int len,
    T alpha,
    const T *x,
    T beta) const {
  assert(nrowdim_ == 1 && ncoldim_ == 1);
  const unsigned int ndofcol = len * nrowblk_;
  for (unsigned int i = 0; i < ndofcol; ++i) { y[i] *= beta; }

  const T *vcrs = val_crs_.data();
  const T *vdia = val_dia_.data();
  const unsigned int *colind = col_ind_.data();
  const unsigned int *rowptr = row_ptr_.data();

  for (unsigned int iblk = 0; iblk < nrowblk_; iblk++) {
    const unsigned int colind0 = colind[iblk];
    const unsigned int colind1 = colind[iblk + 1];
    for (unsigned int icrs = colind0; icrs < colind1; icrs++) {
      assert(icrs < row_ptr_.size());
      const unsigned int jblk0 = rowptr[icrs];
      assert(jblk0 < ncolblk_);
      const T mval0 = alpha * vcrs[icrs];
      for (unsigned int ilen = 0; ilen < len; ilen++) {
        y[iblk * len + ilen] += mval0 * x[jblk0 * len + ilen];
      }
    }
    { // compute diagonal
      const T mval0 = alpha * vdia[iblk];
      for (unsigned int ilen = 0; ilen < len; ilen++) {
        y[iblk * len + ilen] += mval0 * x[iblk * len + ilen];
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMatrixSparse<float>::MatVecDegenerate(float *y,
                                                              unsigned int len,
                                                              float alpha,
                                                              const float *x,
                                                              float beta) const;
template void delfem2::CMatrixSparse<double>::MatVecDegenerate(double *y,
                                                               unsigned int len,
                                                               double alpha,
                                                               const double *x,
                                                               double beta) const;
#endif

// -------------------------------------------------------

// Calc Matrix Vector Product
// {y} = alpha*[A]^T{x} + beta*{y}
template<typename T>
void delfem2::CMatrixSparse<T>::MatTVec(
    T *y,
    T alpha,
    const T *x,
    T beta) const {
  const unsigned int ndofrow = ncoldim_ * ncolblk_;
  for (unsigned int i = 0; i < ndofrow; ++i) { y[i] *= beta; }
  const unsigned int blksize = nrowdim_ * ncoldim_;
  {
    const T *vcrs = val_crs_.data();
    const T *vdia = val_dia_.data();
    const unsigned int *colind = col_ind_.data();
    const unsigned int *rowptr = row_ptr_.data();
    //
    for (unsigned int iblk = 0; iblk < nrowblk_; iblk++) {
      const unsigned int colind0 = colind[iblk];
      const unsigned int colind1 = colind[iblk + 1];
      for (unsigned int icrs = colind0; icrs < colind1; icrs++) {
        assert(icrs < row_ptr_.size());
        const unsigned int jblk0 = rowptr[icrs];
        assert(jblk0 < ncolblk_);
        for (unsigned int idof = 0; idof < nrowdim_; idof++) {
          for (unsigned int jdof = 0; jdof < ncoldim_; jdof++) {
            y[jblk0 * ncoldim_ + jdof] +=
                alpha * vcrs[icrs * blksize + idof * ncoldim_ + jdof] * x[iblk * nrowdim_ + idof];
          }
        }
      }
      for (unsigned int jdof = 0; jdof < ncoldim_; jdof++) {
        for (unsigned int idof = 0; idof < nrowdim_; idof++) {
          y[iblk * ncoldim_ + jdof] +=
              alpha * vdia[iblk * blksize + idof * ncoldim_ + jdof] * x[iblk * nrowdim_ + idof];
        }
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMatrixSparse<float>::MatTVec(
    float *y, float alpha, const float *x, float beta) const;
template void delfem2::CMatrixSparse<double>::MatTVec(
    double *y, double alpha, const double *x, double beta) const;
template void delfem2::CMatrixSparse<std::complex<double>>::MatTVec(
    std::complex<double> *y, std::complex<double> alpha,
    const std::complex<double> *x, std::complex<double> beta) const;
#endif

// -----------------------------------------------------------------

template<typename T>
void delfem2::CMatrixSparse<T>::SetFixedBC_Dia(
    const int *bc_flag,
    T val_dia) {
  assert(!this->val_dia_.empty());
  assert(this->ncolblk_ == this->nrowblk_);
  assert(this->ncoldim_ == this->nrowdim_);
  const int blksize = nrowdim_ * ncoldim_;
  for (unsigned int iblk = 0; iblk < nrowblk_; iblk++) {  // set diagonal
    for (unsigned int ilen = 0; ilen < nrowdim_; ilen++) {
      if (bc_flag[iblk * nrowdim_ + ilen] == 0) continue;
      for (unsigned int jlen = 0; jlen < ncoldim_; jlen++) {
        val_dia_[iblk * blksize + ilen * nrowdim_ + jlen] = 0.0;
        val_dia_[iblk * blksize + jlen * nrowdim_ + ilen] = 0.0;
      }
      val_dia_[iblk * blksize + ilen * nrowdim_ + ilen] = val_dia;
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMatrixSparse<float>::SetFixedBC_Dia(
    const int *bc_flag, float val_dia);
template void delfem2::CMatrixSparse<double>::SetFixedBC_Dia(
    const int *bc_flag, double val_dia);
template void delfem2::CMatrixSparse<std::complex<double>>::SetFixedBC_Dia(
    const int *bc_flag, std::complex<double> val_dia);
#endif

// --------------------------

template<typename T>
void delfem2::CMatrixSparse<T>::SetFixedBC_Row(
    const int *bc_flag) {
  assert(!this->val_dia_.empty());
  assert(this->ncolblk_ == this->nrowblk_);
  assert(this->ncoldim_ == this->nrowdim_);
  const int blksize = nrowdim_ * ncoldim_;
  for (unsigned int iblk = 0; iblk < nrowblk_; iblk++) {  // set row
    for (unsigned int icrs = col_ind_[iblk]; icrs < col_ind_[iblk + 1]; icrs++) {
      for (unsigned int ilen = 0; ilen < nrowdim_; ilen++) {
        if (bc_flag[iblk * nrowdim_ + ilen] == 0) continue;
        for (unsigned int jlen = 0; jlen < ncoldim_; jlen++) {
          val_crs_[icrs * blksize + ilen * nrowdim_ + jlen] = 0.0;
        }
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMatrixSparse<float>::SetFixedBC_Row(
    const int *bc_flag);
template void delfem2::CMatrixSparse<double>::SetFixedBC_Row(
    const int *bc_flag);
template void delfem2::CMatrixSparse<std::complex<double>>::SetFixedBC_Row(
    const int *bc_flag);
#endif

// ---------------------------

template<typename T>
void delfem2::CMatrixSparse<T>::SetFixedBC_Col(
    const int *bc_flag) {
  assert(!this->val_dia_.empty());
  assert(this->ncolblk_ == this->nrowblk_);
  assert(this->ncoldim_ == this->nrowdim_);
  const int blksize = nrowdim_ * ncoldim_;
  for (unsigned int icrs = 0; icrs < row_ptr_.size(); icrs++) {  // set column
    const int jblk1 = row_ptr_[icrs];
    for (unsigned int jlen = 0; jlen < ncoldim_; jlen++) {
      if (bc_flag[jblk1 * ncoldim_ + jlen] == 0) continue;
      for (unsigned int ilen = 0; ilen < nrowdim_; ilen++) {
        val_crs_[icrs * blksize + ilen * nrowdim_ + jlen] = 0.0;
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMatrixSparse<float>::SetFixedBC_Col(const int *bc_flag);
template void delfem2::CMatrixSparse<double>::SetFixedBC_Col(const int *bc_flag);
template void delfem2::CMatrixSparse<std::complex<double>>::SetFixedBC_Col(const int *bc_flag);
#endif

// -----------------------------------------------------------------

DFM2_INLINE void delfem2::MatSparse_ScaleBlk_LeftRight(
    delfem2::CMatrixSparse<double> &mat,
    const double *scale) {
  assert(mat.ncolblk_ == mat.nrowblk_);
  assert(mat.ncoldim_ == mat.nrowdim_);
  const unsigned int nblk = mat.nrowblk_;
  const unsigned int len = mat.nrowdim_;
  const unsigned int blksize = len * len;
  for (unsigned int ino = 0; ino < nblk; ++ino) {
    for (unsigned int icrs0 = mat.col_ind_[ino]; icrs0 < mat.col_ind_[ino + 1]; ++icrs0) {
      const unsigned int jno = mat.row_ptr_[icrs0];
      const double s0 = scale[ino] * scale[jno];
      for (unsigned int i = 0; i < blksize; ++i) { mat.val_crs_[icrs0 * blksize + i] *= s0; }
    }
  }
  if (!mat.val_dia_.empty()) {
    for (unsigned int ino = 0; ino < nblk; ++ino) {
      double s0 = scale[ino] * scale[ino];
      for (unsigned int i = 0; i < blksize; ++i) { mat.val_dia_[ino * blksize + i] *= s0; }
    }
  }
}

DFM2_INLINE void delfem2::MatSparse_ScaleBlkLen_LeftRight(
    delfem2::CMatrixSparse<double> &mat,
    const double *scale) {
  assert(mat.ncolblk_ == mat.nrowblk_);
  assert(mat.ncoldim_ == mat.nrowdim_);
  const unsigned int nblk = mat.nrowblk_;
  const unsigned int len = mat.nrowdim_;
  const unsigned int blksize = len * len;
  for (unsigned int ino = 0; ino < nblk; ++ino) {
    for (unsigned int icrs0 = mat.col_ind_[ino]; icrs0 < mat.col_ind_[ino + 1]; ++icrs0) {
      const unsigned int jno = mat.row_ptr_[icrs0];
      for (unsigned int ilen = 0; ilen < len; ++ilen) {
        for (unsigned int jlen = 0; jlen < len; ++jlen) {
          mat.val_crs_[icrs0 * blksize + ilen * len + jlen] *=
              scale[ino * len + ilen] * scale[jno * len + jlen];
        }
      }
    }
  }
  if (!mat.val_dia_.empty()) {
    for (unsigned int ino = 0; ino < nblk; ++ino) {
      for (unsigned int ilen = 0; ilen < len; ++ilen) {
        for (unsigned int jlen = 0; jlen < len; ++jlen) {
          mat.val_dia_[ino * blksize + ilen * len + jlen] *=
              scale[ino * len + ilen] * scale[ino * len + jlen];
        }
      }
    }
  }
}

DFM2_INLINE double delfem2::CheckSymmetry(
    const delfem2::CMatrixSparse<double> &mat) {
  assert(mat.ncolblk_ == mat.nrowblk_);
  assert(mat.ncoldim_ == mat.nrowdim_);
  const unsigned int blksize = mat.nrowdim_ * mat.ncoldim_;
  const unsigned int nlen = mat.nrowdim_;
  //
  double sum = 0;
  for (unsigned int ino = 0; ino < mat.nrowblk_; ++ino) {
    for (unsigned int icrs0 = mat.col_ind_[ino]; icrs0 < mat.col_ind_[ino + 1]; ++icrs0) {
      unsigned int jno = mat.row_ptr_[icrs0];
      unsigned int icrs1 = mat.col_ind_[jno];
      for (; icrs1 < mat.col_ind_[jno + 1]; ++icrs1) {
        if (mat.row_ptr_[icrs1] == ino) { break; }
      }
      if (icrs1 == mat.col_ind_[jno + 1]) {  // no counterpart
        sum += mats::MatNorm(
            mat.val_crs_.data() + blksize * icrs0, mat.nrowdim_, mat.ncoldim_);
      } else {
        sum += mats::MatNorm_Assym(
            mat.val_crs_.data() + blksize * icrs0, mat.nrowdim_, mat.ncoldim_,
            mat.val_crs_.data() + blksize * icrs1);
      }
    }
    sum += mats::MatNorm_Assym(mat.val_dia_.data() + blksize * ino, nlen);
  }
  return sum;
}
