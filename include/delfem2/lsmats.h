/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_LSMATS_H
#define DFM2_LSMATS_H

#include <vector>
#include <cassert>
#include <complex>
#include <climits>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 * @class sparse matrix class
 * @tparam T float, double and std::complex<double>
 */
template<typename T>
class CMatrixSparse {
 public:
  CMatrixSparse() noexcept
      : nrowblk_(0), ncolblk_(0), nrowdim_(0), ncoldim_(0) {}

  virtual ~CMatrixSparse() {
    this->Clear();
  }

  void Clear() {
    col_ind_.clear();
    row_ptr_.clear();
    val_crs_.clear();
    val_dia_.clear();
    this->nrowblk_ = 0;
    this->nrowdim_ = 0;
    this->ncolblk_ = 0;
    this->ncoldim_ = 0;
  }

  void Initialize(size_t nblk, unsigned int len, bool is_dia) {
    this->nrowblk_ = static_cast<unsigned int>(nblk);
    this->nrowdim_ = len;
    this->ncolblk_ = static_cast<unsigned int>(nblk);
    this->ncoldim_ = len;
    col_ind_.assign(nblk + 1, 0);
    row_ptr_.clear();
    val_crs_.clear();
    if (is_dia) { val_dia_.assign(nblk * len * len, 0.0); }
    else { val_dia_.clear(); }
  }

  /**
   * deep copy (duplicate values) the non-zero pattern and their values
   * @param m
   */
  void operator=(const CMatrixSparse &m) {
    this->nrowblk = m.nrowblk;
    this->nrowdim = m.nrowdim;
    this->ncolblk = m.ncolblk;
    this->ncoldim = m.ncoldim;
    col_ind_ = m.colInd;
    row_ptr_ = m.rowPtr;
    val_crs_ = m.valCrs;
    val_dia_ = m.valDia; // copy value
  }

  void SetPattern(
      const unsigned int *colind,
      [[maybe_unused]] size_t ncolind,
      const unsigned int *rowptr,
      [[maybe_unused]] size_t nrowptr) {
    assert(row_ptr_.empty());
    assert(ncolind == nrowblk_ + 1);
    for (unsigned int iblk = 0; iblk < nrowblk_ + 1; iblk++) { col_ind_[iblk] = colind[iblk]; }
    const unsigned int ncrs = colind[nrowblk_];
    assert(ncrs == nrowptr);
    row_ptr_.resize(ncrs);
    for (unsigned int icrs = 0; icrs < ncrs; icrs++) { row_ptr_[icrs] = rowptr[icrs]; }
    val_crs_.resize(ncrs * nrowdim_ * ncoldim_);
  }

  /**
   * TODO: return void instead of bool
   * @detail the name is the same as the Eigen library
   * @return
   */
  bool setZero() {
    if (val_dia_.size() != 0) {
      assert(nrowdim_ == ncoldim_ && nrowblk_ == ncolblk_);
      const size_t n = val_dia_.size();
      assert(n == nrowdim_ * nrowdim_ * nrowblk_);
      for (unsigned int i = 0; i < n; ++i) { val_dia_[i] = 0; }
    }
    {
      const size_t n = val_crs_.size();
      assert(n == nrowdim_ * ncoldim_ * row_ptr_.size());
      for (unsigned int i = 0; i < n; i++) { val_crs_[i] = 0.0; }
    }
    return true;
  }

  /**
   * @func Matrix vector product as: {y} = alpha * [A]{x} + beta * {y}
   */
  void MatVec(
      T *y,
      T alpha, const T *x,
      T beta) const;

  /**
   * @func Matrix vector product as: {y} = alpha * [A]{x} + beta * {y}.
   *  the sparse matrix is regared as block sparse matrix where each blcok is diagonal
   */
  void MatVecDegenerate(
      T *y,
      unsigned nlen,
      T alpha, const T *x,
      T beta) const;

  /**
   * @func Matrix vector product as: {y} = alpha * [A]^T{x} + beta * {y}
   */
  void MatTVec(
      T *y,
      T alpha, const T *x,
      T beta) const;

  /**
   * @func set fixed bc for diagonal block matrix where( pBCFlag[i] != 0).
   */
  void SetFixedBC_Dia(const int *pBCFlag, T val_dia);

  void SetFixedBC_Col(const int *pBCFlag);

  void SetFixedBC_Row(const int *pBCFlag);

  /**
   * @func     if pBCFlag is *not* 0 for a dof, set all the off-diagonal componenet to zero and set diagonal to one.
   * @details pBCFlag need to have memory at least larger than nlen*nblk
   * This matrix need to be a squared matrix
   */
  void SetFixedBC(const int *pBCFlag) {
    this->SetFixedBC_Dia(pBCFlag, 1.0);
    this->SetFixedBC_Row(pBCFlag);
    this->SetFixedBC_Col(pBCFlag);
  }

  void AddDia(T eps) {
    assert(this->ncolblk_ == this->nrowblk_);
    assert(this->ncoldim_ == this->nrowdim_);
    const int blksize = nrowdim_ * ncoldim_;
    const int nlen = this->nrowdim_;
    if (val_dia_.empty()) { return; }
    for (unsigned int ino = 0; ino < nrowblk_; ++ino) {
      for (int ilen = 0; ilen < nlen; ++ilen) {
        val_dia_[ino * blksize + ilen * nlen + ilen] += eps;
      }
    }
  }

  /**
   * @func add vector to diagonal component
   * @param[in] lm a lumped mass vector with size of nblk
   * @param[in] scale scaling factor for the lumped mass (typically 1/dt^2).
   * @details the matrix need to be square matrix
   */
  void AddDia_LumpedMass(const T *lm, double scale) {
    assert(this->nblk_row == this->nblk_col);
    assert(this->len_row == this->nrowdim);
    const int blksize = nrowdim_ * ncoldim_;
    const int nlen = this->nrowdim;
    if (val_dia_.empty()) { return; }
    for (unsigned int iblk = 0; iblk < nrowblk_; ++iblk) {
      for (int ilen = 0; ilen < nlen; ++ilen) {
        val_dia_[iblk * blksize + ilen * nlen + ilen] += lm[iblk];
      }
    }
  }

 public:
  unsigned int nrowblk_;
  unsigned int ncolblk_;
  unsigned int nrowdim_;
  unsigned int ncoldim_;

  /**
   * @param colInd indeces where the row starts in CRS data structure
   */
  std::vector<unsigned int> col_ind_;

  /**
   * @param rowPtr index of CRS data structure
   */
  std::vector<unsigned int> row_ptr_;
  std::vector<T> val_crs_;
  std::vector<T> val_dia_;
};

// ----------------------------------------------

template<typename T>
bool Mearge(
    CMatrixSparse<T> &A,
    size_t nrow, const unsigned int *aIpRow,
    size_t ncol, const unsigned int *aIpCol,
    unsigned int blksize, const T *emat,
    std::vector<unsigned int> &merge_buffer) {
  assert(!A.val_crs_.empty());
  assert(!A.val_dia_.empty());
  assert(blksize == A.nrowdim_ * A.ncoldim_);
  merge_buffer.resize(A.ncolblk_);
  const unsigned int *colind = A.col_ind_.data();
  const unsigned int *rowptr = A.row_ptr_.data();
  T *vcrs = A.val_crs_.data();
  T *vdia = A.val_dia_.data();
  for (unsigned int irow = 0; irow < nrow; irow++) {
    const unsigned int iblk1 = aIpRow[irow];
    assert(iblk1 < A.nrowblk_);
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.row_ptr_.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jcol = 0; jcol < ncol; jcol++) {
      const unsigned int jblk1 = aIpCol[jcol];
      assert(jblk1 < A.ncolblk_);
      if (iblk1 == jblk1) {  // Marge Diagonal
        const T *pval_in = &emat[(irow * ncol + irow) * blksize];
        T *pval_out = &vdia[iblk1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      } else {  // Marge Non-Diagonal
        if (merge_buffer[jblk1] == UINT_MAX) {
          assert(0);
          return false;
        }
        assert(merge_buffer[jblk1] < A.row_ptr_.size());
        const unsigned int jpsup1 = merge_buffer[jblk1];
        assert(A.row_ptr_[jpsup1] == jblk1);
        const T *pval_in = &emat[(irow * ncol + jcol) * blksize];
        T *pval_out = &vcrs[jpsup1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
    }
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.row_ptr_.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = UINT_MAX;
    }
  }
  return true;
}

template<int nrow, int ncol, int ndimrow, int ndimcol, typename T>
bool Merge(
    CMatrixSparse<T> &A,
    const unsigned int aIpRow[nrow],
    const unsigned int aIpCol[ncol],
    const T emat[nrow][ncol][ndimrow][ndimcol],
    std::vector<unsigned int> &merge_buffer) {
  assert(!A.val_crs_.empty());
  assert(!A.val_dia_.empty());
  const unsigned int blksize = ndimrow * ndimcol;
  merge_buffer.resize(A.ncolblk_, UINT_MAX);
  const unsigned int *colind = A.col_ind_.data();
  const unsigned int *rowptr = A.row_ptr_.data();
  T *vcrs = A.val_crs_.data();
  T *vdia = A.val_dia_.data();
  for (unsigned int irow = 0; irow < nrow; irow++) {
    const unsigned int iblk1 = aIpRow[irow];
    assert(iblk1 < A.nrowblk_);
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.row_ptr_.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jcol = 0; jcol < ncol; jcol++) {
      const unsigned int jblk1 = aIpCol[jcol];
      assert(jblk1 < A.ncolblk_);
      if (iblk1 == jblk1) {  // Marge Diagonal
        const T *pval_in = &emat[irow][irow][0][0];
        T *pval_out = &vdia[iblk1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      } else {  // Marge Non-Diagonal
        if (merge_buffer[jblk1] == UINT_MAX) {
          assert(0);
          return false;
        }
        assert(merge_buffer[jblk1] < A.row_ptr_.size());
        const int jpsup1 = merge_buffer[jblk1];
        assert(A.row_ptr_[jpsup1] == jblk1);
        const T *pval_in = &emat[irow][jcol][0][0];
        T *pval_out = &vcrs[jpsup1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
    }
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.row_ptr_.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = UINT_MAX;
    }
  }
  return true;
}

template<int nrow, int ncol, typename T>
bool Merge(
    CMatrixSparse<T> &A,
    const unsigned int *aIpRow,
    const unsigned int *aIpCol,
    const T emat[nrow][ncol],
    std::vector<unsigned int> &merge_buffer) {
  assert(!A.val_crs_.empty());
  assert(!A.val_dia_.empty());
  merge_buffer.resize(A.ncolblk_);
  const unsigned int *colind = A.col_ind_.data();
  const unsigned int *rowptr = A.row_ptr_.data();
  T *vcrs = A.val_crs_.data();
  T *vdia = A.val_dia_.data();
  for (unsigned int irow = 0; irow < nrow; irow++) {
    const unsigned int iblk1 = aIpRow[irow];
    assert(iblk1 < A.nrowblk_);
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.row_ptr_.size());
      const unsigned int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jcol = 0; jcol < ncol; jcol++) {
      const unsigned int jblk1 = aIpCol[jcol];
      assert(jblk1 < A.ncolblk_);
      if (iblk1 == jblk1) {  // Marge Diagonal
        vdia[iblk1] += emat[irow][irow];
      } else {  // Marge Non-Diagonal
        if (merge_buffer[jblk1] == UINT_MAX) {
          assert(0);
          return false;
        }
        assert(merge_buffer[jblk1] < A.row_ptr_.size());
        const unsigned int jpsup1 = merge_buffer[jblk1];
        assert(A.row_ptr_[jpsup1] == jblk1);
        vcrs[jpsup1] += emat[irow][jcol];
      }
    }
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.row_ptr_.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = UINT_MAX;
    }
  }
  return true;
}

DFM2_INLINE double CheckSymmetry(
    const delfem2::CMatrixSparse<double> &mat);

DFM2_INLINE void MatSparse_ScaleBlk_LeftRight(
    delfem2::CMatrixSparse<double> &mat,
    const double *scale);

DFM2_INLINE void MatSparse_ScaleBlkLen_LeftRight(
    delfem2::CMatrixSparse<double> &mat,
    const double *scale);

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/lsmats.cpp"
#endif

#endif // MATDIA_CRS_H
