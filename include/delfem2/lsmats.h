/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_LSMATS_H
#define DFM2_LSMATS_H

#include "delfem2/dfm2_inline.h"
#include <vector>
#include <cassert>
#include <complex>
#include <climits>

namespace delfem2 {

/**
 * @class sparse matrix class
 * @tparam T float, double and std::complex<double>
 */
template<typename T>
class CMatrixSparse {
public:
  CMatrixSparse() noexcept
  : nrowblk(0), ncolblk(0), nrowdim(0), ncoldim(0)
  {}

  virtual ~CMatrixSparse() {
    this->Clear();
  }
  
  void Clear(){
    colInd.clear();
    rowPtr.clear();
    valCrs.clear();
    valDia.clear();
    this->nrowblk = 0;
    this->nrowdim = 0;
    this->ncolblk = 0;
    this->ncoldim = 0;
  }

  void Initialize(unsigned int nblk, unsigned int len, bool is_dia) {
    this->nrowblk = nblk;
    this->nrowdim = len;
    this->ncolblk = nblk;
    this->ncoldim = len;
    colInd.assign(nblk + 1, 0);
    rowPtr.clear();
    valCrs.clear();
    if (is_dia) { valDia.assign(nblk * len * len, 0.0); }
    else { valDia.clear(); }
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
    colInd = m.colInd;
    rowPtr = m.rowPtr;
    valCrs = m.valCrs;
    valDia = m.valDia; // copy value
  }

  void SetPattern(
      const unsigned int *colind,
      size_t ncolind,
      const unsigned int *rowptr,
      size_t nrowptr)
  {
    assert(rowPtr.empty());
    assert(ncolind == nrowblk + 1);
    for (unsigned int iblk = 0; iblk < nrowblk + 1; iblk++) { colInd[iblk] = colind[iblk]; }
    const unsigned int ncrs = colind[nrowblk];
    assert(ncrs == nrowptr);
    rowPtr.resize(ncrs);
    for (unsigned int icrs = 0; icrs < ncrs; icrs++) { rowPtr[icrs] = rowptr[icrs]; }
    valCrs.resize(ncrs * nrowdim * ncoldim);
  }

  /**
   * @detail the name is the same as the Eigen library
   * @return
   */
  bool setZero() {
    if (valDia.size() != 0) {
      assert(nrowdim == ncoldim && nrowblk == ncolblk);
      const unsigned int n = valDia.size();
      assert(n == nrowdim * nrowdim * nrowblk);
      for (unsigned int i = 0; i < n; ++i) { valDia[i] = 0; }
    }
    {
      const unsigned int n = valCrs.size();
      assert(n == nrowdim * ncoldim * rowPtr.size());
      for (unsigned int i = 0; i < n; i++) { valCrs[i] = 0.0; }
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
  void SetFixedBC(const int *pBCFlag){
    this->SetFixedBC_Dia(pBCFlag,1.0);
    this->SetFixedBC_Row(pBCFlag);
    this->SetFixedBC_Col(pBCFlag);
  }

  void AddDia(T eps) {
    assert(this->ncolblk == this->nrowblk);
    assert(this->ncoldim == this->nrowdim);
    const int blksize = nrowdim * ncoldim;
    const int nlen = this->nrowdim;
    if (valDia.empty()) { return; }
    for (unsigned int ino = 0; ino < nrowblk; ++ino) {
      for (int ilen = 0; ilen < nlen; ++ilen) {
        valDia[ino * blksize + ilen * nlen + ilen] += eps;
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
    const int blksize = nrowdim * ncoldim;
    const int nlen = this->nrowdim;
    if (valDia.empty()) { return; }
    for (unsigned int iblk = 0; iblk < nrowblk; ++iblk) {
      for (int ilen = 0; ilen < nlen; ++ilen) {
        valDia[iblk * blksize + ilen * nlen + ilen] += lm[iblk];
      }
    }
  }

public:
  unsigned int nrowblk;
  unsigned int ncolblk;
  unsigned int nrowdim;
  unsigned int ncoldim;
  
  /**
   * @param colInd indeces where the row starts in CRS data structure
   */
  std::vector<unsigned int> colInd;
  
  /**
   * @param rowPtr index of CRS data structure
   */
  std::vector<unsigned int> rowPtr;
  std::vector<T> valCrs;
  std::vector<T> valDia;
};

// ----------------------------------------------

template <typename T>
bool Mearge(
    CMatrixSparse<T>& A,
    unsigned int nrow, const unsigned int *aIpRow,
    unsigned int ncol, const unsigned int *aIpCol,
    unsigned int blksize, const T *emat,
    std::vector<unsigned int> &merge_buffer)
{
  assert(!A.valCrs.empty());
  assert(!A.valDia.empty());
  assert(blksize == A.nrowdim * A.ncoldim);
  merge_buffer.resize(A.ncolblk);
  const unsigned int *colind = A.colInd.data();
  const unsigned int *rowptr = A.rowPtr.data();
  T *vcrs = A.valCrs.data();
  T *vdia = A.valDia.data();
  for (unsigned int irow = 0; irow < nrow; irow++) {
    const unsigned int iblk1 = aIpRow[irow];
    assert(iblk1 < A.nrowblk);
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jcol = 0; jcol < ncol; jcol++) {
      const unsigned int jblk1 = aIpCol[jcol];
      assert(jblk1 < A.ncolblk);
      if (iblk1 == jblk1) {  // Marge Diagonal
        const T *pval_in = &emat[(irow * ncol + irow) * blksize];
        T *pval_out = &vdia[iblk1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
      else {  // Marge Non-Diagonal
        if (merge_buffer[jblk1] == UINT_MAX) {
          assert(0);
          return false;
        }
        assert( merge_buffer[jblk1] < A.rowPtr.size());
        const unsigned int jpsup1 = merge_buffer[jblk1];
        assert(A.rowPtr[jpsup1] == jblk1);
        const T *pval_in = &emat[(irow * ncol + jcol) * blksize];
        T *pval_out = &vcrs[jpsup1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
    }
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = UINT_MAX;
    }
  }
  return true;
}

template <int nrow, int ncol, int ndimrow, int ndimcol, typename T>
bool Merge(
    CMatrixSparse<T>& A,
    const unsigned int* aIpRow,
    const unsigned int* aIpCol,
    const T emat[nrow][ncol][ndimrow][ndimcol],
    std::vector<unsigned int>& merge_buffer)
{
  assert(!A.valCrs.empty());
  assert(!A.valDia.empty());
  const unsigned int blksize = ndimrow * ndimcol;
  merge_buffer.resize(A.ncolblk,UINT_MAX);
  const unsigned int *colind = A.colInd.data();
  const unsigned int *rowptr = A.rowPtr.data();
  T *vcrs = A.valCrs.data();
  T *vdia = A.valDia.data();
  for (unsigned int irow = 0; irow < nrow; irow++) {
    const unsigned int iblk1 = aIpRow[irow];
    assert(iblk1 < A.nrowblk);
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jcol = 0; jcol < ncol; jcol++) {
      const unsigned int jblk1 = aIpCol[jcol];
      assert(jblk1 < A.ncolblk);
      if (iblk1 == jblk1) {  // Marge Diagonal
        const T *pval_in = &emat[irow][irow][0][0];
        T *pval_out = &vdia[iblk1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
      else {  // Marge Non-Diagonal
        if (merge_buffer[jblk1] == UINT_MAX) {
          assert(0);
          return false;
        }
        assert( merge_buffer[jblk1] < A.rowPtr.size() );
        const int jpsup1 = merge_buffer[jblk1];
        assert(A.rowPtr[jpsup1] == jblk1);
        const T *pval_in = &emat[irow][jcol][0][0];
        T *pval_out = &vcrs[jpsup1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
    }
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = UINT_MAX;
    }
  }
  return true;
}

template <int nrow, int ncol, typename T>
bool Merge(
    CMatrixSparse<T>& A,
    const unsigned int* aIpRow,
    const unsigned int* aIpCol,
    const T emat[nrow][ncol],
    std::vector<unsigned int>& merge_buffer)
{
  assert(!A.valCrs.empty());
  assert(!A.valDia.empty());
  merge_buffer.resize(A.ncolblk);
  const unsigned int *colind = A.colInd.data();
  const unsigned int *rowptr = A.rowPtr.data();
  T *vcrs = A.valCrs.data();
  T *vdia = A.valDia.data();
  for (unsigned int irow = 0; irow < nrow; irow++) {
    const unsigned int iblk1 = aIpRow[irow];
    assert(iblk1 < A.nrowblk);
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const unsigned int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jcol = 0; jcol < ncol; jcol++) {
      const unsigned int jblk1 = aIpCol[jcol];
      assert(jblk1 < A.ncolblk);
      if (iblk1 == jblk1) {  // Marge Diagonal
        vdia[iblk1] += emat[irow][irow];
      }
      else {  // Marge Non-Diagonal
        if (merge_buffer[jblk1] == UINT_MAX) {
          assert(0);
          return false;
        }
        assert( merge_buffer[jblk1] < A.rowPtr.size() );
        const unsigned int jpsup1 = merge_buffer[jblk1];
        assert(A.rowPtr[jpsup1] == jblk1);
        vcrs[jpsup1] += emat[irow][jcol];
      }
    }
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const int jblk1 = rowptr[jpsup];
      merge_buffer[jblk1] = UINT_MAX;
    }
  }
  return true;
}

DFM2_INLINE double CheckSymmetry(
    const delfem2::CMatrixSparse<double> &mat);

DFM2_INLINE void SetMasterSlave(
    delfem2::CMatrixSparse<double> &mat,
    const unsigned int *aMSFlag);

DFM2_INLINE void MatSparse_ScaleBlk_LeftRight(
    delfem2::CMatrixSparse<double> &mat,
    const double *scale);

DFM2_INLINE void MatSparse_ScaleBlkLen_LeftRight(
    delfem2::CMatrixSparse<double> &mat,
    const double *scale);

} // delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/lsmats.cpp"
#endif
  
#endif // MATDIA_CRS_H
