
#ifndef DFM2_EIGEN_LS_SPARSE_H
#define DFM2_EIGEN_LS_SPARSE_H

#include <Eigen/Core>
#include <vector>

namespace delfem2 {

template<class MAT, class ALLOCATOR, int RowsActive = MAT::RowsAtCompileTime>
class CMatrixSparseBlock {
public:
  CMatrixSparseBlock() noexcept
      : nrowblk(0), ncolblk(0) {
    static_assert( RowsActive <= MAT::RowsAtCompileTime, "active dim should be no greater than the actual dim");
  }

  void Initialize(unsigned int nblk) {
    this->nrowblk = nblk;
    this->ncolblk = nblk;
    colInd.assign(nblk + 1, 0);
    rowPtr.clear();
    valCrs.clear();
    valDia.resize(nblk);
  }

  void SetPattern(
      const unsigned int *colind,
      size_t ncolind,
      const unsigned int *rowptr,
      size_t nrowptr) {
    assert(rowPtr.empty() && ncolind == nrowblk + 1);
    for (unsigned int iblk = 0; iblk < nrowblk + 1; iblk++) { colInd[iblk] = colind[iblk]; }
    const unsigned int ncrs = colind[nrowblk];
    assert(ncrs == nrowptr);
    rowPtr.resize(ncrs);
    for (unsigned int icrs = 0; icrs < ncrs; icrs++) { rowPtr[icrs] = rowptr[icrs]; }
    valCrs.resize(ncrs);
  }

  /**
   * @detail the name is the same as the Eigen library
   * @return
   */
  void setZero() {
    assert(nrowblk == ncolblk && valCrs.size() == rowPtr.size());
    for (auto &m : valDia) { m.setZero(); }
    for (auto &m : valCrs) { m.setZero(); }
    // --------
    if( RowsActive == MAT::RowsAtCompileTime ){ return; }
    for (auto &m : valDia) {
      for (int i = RowsActive; i < MAT::RowsAtCompileTime; ++i) {
        m(i, i) = 1;
      }
    }
  }

public:
  unsigned int nrowblk;
  unsigned int ncolblk;
  std::vector<unsigned int> colInd;
  std::vector<unsigned int> rowPtr;
  std::vector<MAT, ALLOCATOR> valCrs;
  std::vector<MAT, ALLOCATOR> valDia;
};

template<int nrow, int ncol,
    int ndimrow, int ndimcol,
    typename REAL, class MAT, class ALLOCATOR, int RowsActive=MAT::RowsAtCompileTime>
bool Merge(
    CMatrixSparseBlock<MAT, ALLOCATOR, RowsActive> &A,
    const unsigned int aIpRow[nrow],
    const unsigned int aIpCol[ncol],
    const REAL emat[nrow][ncol][ndimrow][ndimcol],
    std::vector<unsigned int> &merge_buffer) {
  static_assert(RowsActive<=MAT::RowsAtCompileTime,
      "the active dim must be not greater than the matrix size");
  merge_buffer.resize(A.ncolblk, UINT_MAX);
  for (unsigned int irow = 0; irow < nrow; irow++) {
    const unsigned int iblk1 = aIpRow[irow];
    assert(iblk1 < A.nrowblk);
    for (unsigned int jpsup = A.colInd[iblk1]; jpsup < A.colInd[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const int jblk1 = A.rowPtr[jpsup];
      merge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jcol = 0; jcol < ncol; jcol++) {
      const unsigned int jblk1 = aIpCol[jcol];
      assert(jblk1 < A.ncolblk);
      if (iblk1 == jblk1) {  // Marge Diagonal
        for (int i = 0; i < ndimrow; ++i) {
          for (int j = 0; j < ndimcol; ++j) {
            A.valDia[iblk1](i, j) += emat[irow][irow][i][j];
          }
        }
      } else {  // Marge Non-Diagonal
        assert(merge_buffer[jblk1] < A.rowPtr.size());
        const unsigned int jpsup1 = merge_buffer[jblk1];
        assert(A.rowPtr[jpsup1] == jblk1);
        for (int i = 0; i < ndimrow; ++i) {
          for (int j = 0; j < ndimcol; ++j) {
            A.valCrs[jpsup1](i, j) += emat[irow][jcol][i][j];
          }
        }
      }
    }
    for (unsigned int jpsup = A.colInd[iblk1]; jpsup < A.colInd[iblk1 + 1]; jpsup++) {
      assert(jpsup < A.rowPtr.size());
      const int jblk1 = A.rowPtr[jpsup];
      merge_buffer[jblk1] = UINT_MAX;
    }
  }
  return true;
}

template<typename REAL,
    class MAT, class ALLOCATOR, int RowsActive = MAT::RowsAtCompileTime >
void SetFixedBC_Dia(
    CMatrixSparseBlock<MAT, ALLOCATOR, RowsActive> &A,
    const int *bc_flag,
    REAL val_dia) {
  constexpr int nrowdim = MAT::RowsAtCompileTime;
  constexpr int ncoldim = MAT::ColsAtCompileTime;
  static_assert(nrowdim == ncoldim,
      "matrix should be diagonal");
  assert(A.ncolblk == A.nrowblk && ncoldim == nrowdim);
  for (unsigned int iblk = 0; iblk < A.nrowblk; iblk++) { // set diagonal
    for (unsigned int ilen = 0; ilen < RowsActive; ++ilen) {
      if (bc_flag[iblk * RowsActive + ilen] == 0){ continue; }
      for (unsigned int jlen = 0; jlen < RowsActive; ++jlen) {
        A.valDia[iblk](ilen, jlen) = 0;
        A.valDia[iblk](jlen, ilen) = 0;
      }
      A.valDia[iblk](ilen, ilen) = val_dia;
    }
  }
}

template<class MAT, class ALLOCATOR, int RowsActive = MAT::RowsAtCompileTime >
void SetFixedBC_Row(
    CMatrixSparseBlock<MAT, ALLOCATOR, RowsActive> &A,
    const int *bc_flag) {
  static_assert(RowsActive<=MAT::RowsAtCompileTime,
      "active dimension should be smaller than the matrix row");
  for (unsigned int iblk = 0; iblk < A.nrowblk; iblk++) { // set row
    for (unsigned int icrs = A.colInd[iblk]; icrs < A.colInd[iblk + 1]; icrs++) {
      for (unsigned int ilen = 0; ilen < RowsActive; ilen++) {
        if (bc_flag[iblk * RowsActive + ilen] == 0) continue;
        for (unsigned int jlen = 0; jlen < RowsActive; jlen++) {
          A.valCrs[icrs](ilen, jlen) = 0;
        }
      }
    }
  }
}

template<class MAT, class ALLOCATOR, int RowsActive = MAT::RowsAtCompileTime>
void SetFixedBC_Col(
    CMatrixSparseBlock<MAT, ALLOCATOR, RowsActive> &A,
    const int *bc_flag) {
  static_assert(RowsActive<=MAT::RowsAtCompileTime,
      "number of active rows should not be greater than the actual row");
  for (unsigned int icrs = 0; icrs < A.rowPtr.size(); icrs++) { // set column
    const unsigned int jblk1 = A.rowPtr[icrs];
    for (unsigned int jlen = 0; jlen < RowsActive; jlen++) {
      if (bc_flag[jblk1 * RowsActive + jlen] == 0){ continue; }
      for (unsigned int ilen = 0; ilen < RowsActive; ilen++) {
        A.valCrs[icrs](ilen, jlen) = 0;
      }
    }
  }
}

// --------------------------------
// Eigen dependency from here

template<typename REAL,
    class MAT, class ALLOCATOR, int RowsActive = MAT::RowsAtCompileTime >
void AddMatVec(
    Eigen::Matrix<REAL,-1,1> &lhs,
    REAL beta,
    REAL alpha,
    const CMatrixSparseBlock<MAT, ALLOCATOR, RowsActive> &A,
    const Eigen::Matrix<REAL,-1,1> &rhs)
{
  assert(lhs.rows() == rhs.rows());
  constexpr int nrowdim = MAT::RowsAtCompileTime;
  constexpr int ncoldim = MAT::ColsAtCompileTime;
  assert(lhs.rows()%nrowdim==0);
  lhs *= beta;
  for (unsigned int iblk = 0; iblk < A.nrowblk; iblk++) {
    for (unsigned int icrs = A.colInd[iblk]; icrs < A.colInd[iblk + 1]; icrs++) {
      assert(icrs < A.rowPtr.size());
      const unsigned int jblk0 = A.rowPtr[icrs];
      assert(jblk0 < A.ncolblk);
      lhs.template segment<nrowdim>(iblk * nrowdim) += alpha * A.valCrs[icrs] * rhs.template segment<ncoldim>(jblk0 * ncoldim); // SIMD?
    }
    lhs.template segment<nrowdim>(iblk * nrowdim) += alpha * A.valDia[iblk] * rhs.template segment<ncoldim>(iblk * ncoldim); // SIMD?
  }
}

template<typename REAL,
    class MAT, class ALLOCATOR, int RowsActive = MAT::RowsAtCompileTime >
void AddMatVec(
    Eigen::Matrix<REAL,-1,MAT::RowsAtCompileTime,Eigen::RowMajor,-1,MAT::RowsAtCompileTime> &lhs,
    REAL beta,
    REAL alpha,
    const CMatrixSparseBlock<MAT, ALLOCATOR, RowsActive> &A,
    const Eigen::Matrix<REAL,-1,MAT::ColsAtCompileTime,Eigen::RowMajor,-1,MAT::ColsAtCompileTime> &rhs) {
  assert(lhs.rows() == rhs.rows());
  assert(lhs.rows()%MAT::RowsAtCompileTime==0);
  lhs *= beta;
  for (unsigned int iblk = 0; iblk < A.nrowblk; iblk++) {
    for (unsigned int icrs = A.colInd[iblk]; icrs < A.colInd[iblk + 1]; icrs++) {
      assert(icrs < A.rowPtr.size());
      const unsigned int jblk0 = A.rowPtr[icrs];
      assert(jblk0 < A.ncolblk);
      lhs.row(iblk) += alpha * A.valCrs[icrs] * rhs.row(jblk0).transpose(); // SIMD?
    }
    lhs.row(iblk) += alpha * A.valDia[iblk] * rhs.row(iblk).transpose(); // SIMD?
  }
}

}

#endif