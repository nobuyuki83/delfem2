/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_MATS_H
#define DFM2_MATS_H

#include <vector>
#include <cassert>
#include <complex>

namespace delfem2 {

/**
 * @class sparse matrix class
 * @tparam T float, double and std::complex<double>
 */
template<typename T>
class CMatrixSparse {
public:
  CMatrixSparse() : nblk_col(0), nblk_row(0), len_col(0), len_row(0) {}

  virtual ~CMatrixSparse() {
    colInd.clear();
    rowPtr.clear();
    valCrs.clear();
    valDia.clear();
  }

  void Initialize(unsigned int nblk, unsigned int len, bool is_dia) {
    this->nblk_col = nblk;
    this->len_col = len;
    this->nblk_row = nblk;
    this->len_row = len;
    colInd.assign(nblk + 1, 0);
    rowPtr.clear();
    valCrs.clear();
    if (is_dia) { valDia.assign(nblk * len * len, 0.0); }
    else { valDia.clear(); }
  }

  void operator=(const CMatrixSparse &m) {
    this->nblk_col = m.nblk_col;
    this->len_col = m.len_col;
    this->nblk_row = m.nblk_row;
    this->len_row = m.len_row;
    colInd = m.colInd;
    rowPtr = m.rowPtr;
    valCrs = m.valCrs;
    valDia = m.valDia; // copy value
  }

  void SetPattern(const unsigned int *colind, unsigned int ncolind,
                  const unsigned int *rowptr, unsigned int nrowptr) {
    assert(rowPtr.empty());
    assert(ncolind == nblk_col + 1);
    for (unsigned int iblk = 0; iblk < nblk_col + 1; iblk++) { colInd[iblk] = colind[iblk]; }
    const unsigned int ncrs = colind[nblk_col];
    assert(ncrs == nrowptr);
    rowPtr.resize(ncrs);
    for (unsigned int icrs = 0; icrs < ncrs; icrs++) { rowPtr[icrs] = rowptr[icrs]; }
    valCrs.resize(ncrs * len_col * len_row);
  }

  bool SetZero() {
    if (valDia.size() != 0) {
      assert(len_col == len_row);
      assert(nblk_col == nblk_row);
      const unsigned int n = valDia.size();
      assert(n == len_col * len_col * nblk_col);
      for (unsigned int i = 0; i < n; ++i) { valDia[i] = 0; }
    }
    {
      const unsigned int n = valCrs.size();
      assert(n == len_col * len_row * rowPtr.size());
      for (unsigned int i = 0; i < n; i++) { valCrs[i] = 0.0; }
    }
    return true;
  }

  bool Mearge(unsigned int nblkel_col, const unsigned int *blkel_col,
              unsigned int nblkel_row, const unsigned int *blkel_row,
              unsigned int blksize, const T *emat,
              std::vector<int> &m_marge_tmp_buffer);

  /**
   * @func Matrix vector product as: {y} = alpha * [A]{x} + beta * {y}
   */
  void MatVec(T *y,
              T alpha, const T *x,
              T beta) const;
  /**
   * @func Matrix vector product as: {y} = alpha * [A]^T{x} + beta * {y}
   */
  void MatTVec(T *y,
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
    assert(this->nblk_row == this->nblk_col);
    assert(this->len_row == this->len_col);
    const int blksize = len_col * len_row;
    const int nlen = this->len_col;
    if (valDia.empty()) { return; }
    for (unsigned int ino = 0; ino < nblk_col; ++ino) {
      for (int ilen = 0; ilen < nlen; ++ilen) {
        valDia[ino * blksize + ilen * nlen + ilen] += eps;
      }
    }
  }

  /**
   * @func add vector to diagonal component
   * @param lm        (in) a lumped mass vector with size of nblk
   * @param scale (in) scaling factor for the lumped mass (typically 1/dt^2).
   * @details the matrix need to be square matrix
   */
  void AddDia_LumpedMass(const T *lm, double scale) {
    assert(this->nblk_row == this->nblk_col);
    assert(this->len_row == this->len_col);
    const int blksize = len_col * len_row;
    const int nlen = this->len_col;
    if (valDia.empty()) { return; }
    for (unsigned int iblk = 0; iblk < nblk_col; ++iblk) {
      for (int ilen = 0; ilen < nlen; ++ilen) {
        valDia[iblk * blksize + ilen * nlen + ilen] += lm[iblk];
      }
    }
  }

public:
  unsigned int nblk_col;
  unsigned int nblk_row;
  unsigned int len_col;
  unsigned int len_row;
  /**
   * @param colInd indeces where the row starts in CRS data structure
   */
  std::vector<unsigned int> colInd;
  /**
   * @param row index of CRS data structure
   */
  std::vector<unsigned int> rowPtr;
  std::vector<T> valCrs;
  std::vector<T> valDia;
};

double CheckSymmetry(const delfem2::CMatrixSparse<double> &mat);
  
void SetMasterSlave(delfem2::CMatrixSparse<double> &mat, const int *aMSFlag);

void MatSparse_ScaleBlk_LeftRight(delfem2::CMatrixSparse<double> &mat,
                                  const double *scale);

void MatSparse_ScaleBlkLen_LeftRight(delfem2::CMatrixSparse<double> &mat,
                                     const double *scale);

} // delfem2
  
#endif // MATDIA_CRS_H
