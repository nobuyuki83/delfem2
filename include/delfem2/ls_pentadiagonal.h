#ifndef DFM2_LS_PENTADIAGONAL_H
#define DFM2_LS_PENTADIAGONAL_H

#include "delfem2/matn.h"

namespace delfem2 {

/**
 *
 * @tparam ndim size of block
 */
template<unsigned int ndim>
class BlockPentaDiagonalMatrix {
 public:
  // initialize with block size n
  BlockPentaDiagonalMatrix() { nblk = 0; }
  ~BlockPentaDiagonalMatrix() {}

  void Initialize(size_t n) {
    this->nblk = n;
    assert(nblk >= 4);
    const unsigned int nblks = 3 + 4 + 5 * (nblk - 4) + 4 + 3;
    v.resize(nblks * ndim * ndim);
  }

  void Clear() {
    nblk = 0;
    v.clear();
  }

  /**
   * 0  1  2  #  #  #
   * 3  4  5  6  #  #
   * 7  8  9 10 11  #
   * # 12 13 14 15 16
   * #  # 17 18 19 20
   * #  #  # 21 22 23
   */
  double *GetValuePointer(int iblk, int jblk) {
    constexpr unsigned int blksize = ndim * ndim;
    if (iblk < 0 || iblk >= nblk) { return nullptr; }
    if (jblk < 0 || jblk >= nblk) { return nullptr; }
    if (iblk - jblk < -2 || iblk - jblk > +2) { return nullptr; }
    if (iblk < 2) { return v.data() + (iblk * 3 + jblk) * blksize; }
    if (iblk == nblk - 1) { return v.data() + (iblk * 4 - 2 + jblk) * blksize; }
    return v.data() + (iblk * 4 - 1 + jblk) * blksize;
  }

  /**
   * named after Eigen library
   */
  void setZero() { v.assign(v.size(), 0.0); }

  /**
   * marge element stiffness matrix to the position (idiv,idiv)
   */
  void Merge(unsigned int iblk, double eM[][3][ndim][ndim]) {
    assert(iblk < nblk - 2);
    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        double *p = GetValuePointer(iblk + i, iblk + j);
        if (p == nullptr) { continue; }
        const double *q = (&eM[i][j][0][0]);
        for (unsigned int k = 0; k < ndim * ndim; ++k) { p[k] += q[k]; }
      }
    }
  }

  //! define fixed boudnary condition
  void FixBC(unsigned int iblk, unsigned int idim);

  // execute LU factorization
  void Decompose();

  // solve matrix
  void Solve(std::vector<double> &res);

 private:
  int nblk;
  std::vector<double> v;
};

//! define fixed boudnary condition
template<unsigned int ndim>
void BlockPentaDiagonalMatrix<ndim>::FixBC(
    unsigned int iblk,
    unsigned int idim) {
  assert(idim < ndim && int(iblk) < nblk);
  const int iblk0 = static_cast<int>(iblk);
  for (int jblk = iblk0 - 2; jblk <= iblk0 + 2; ++jblk) {
    double *p = this->GetValuePointer(iblk, jblk);
    if (p == nullptr) { continue; }
    for (unsigned int jdim = 0; jdim < ndim; ++jdim) {
      p[idim * ndim + jdim] = 0;
    }
  }
  for (int jblk = iblk0 - 2; jblk <= iblk0 + 2; ++jblk) {
    double *p = this->GetValuePointer(jblk, iblk);
    if (p == nullptr) { continue; }
    for (unsigned int jdim = 0; jdim < ndim; ++jdim) {
      p[jdim * ndim + idim] = 0;
    }
  }
  GetValuePointer(iblk, iblk)[idim * ndim + idim] = 1;
}

template<unsigned int ndim>
void BlockPentaDiagonalMatrix<ndim>::Decompose() {
  constexpr unsigned int blksize = ndim * ndim;
  double tmpBlk[blksize];
  for (int iblk = 0; iblk < nblk; iblk++) {
    if (iblk != 0) {
      const double *pVal_ik = GetValuePointer(iblk, iblk - 1);
      const double *pVal_kj = GetValuePointer(iblk - 1, iblk);
      double *pVal_ij = GetValuePointer(iblk, iblk);
      Sub_MatMat<double, ndim>(pVal_ij, pVal_ik, pVal_kj);
    }
    {   // calc inverse of diagonal
      double *pVal_ii = GetValuePointer(iblk, iblk);
      Inverse_Matrix<double, ndim>(pVal_ii);
    }
    // [U] = [1/D][U]
    if (iblk != nblk - 1) {
      double *pVal_ij = GetValuePointer(iblk, iblk + 1);
      const double *pVal_ii = GetValuePointer(iblk, iblk);
      for (unsigned int i = 0; i < blksize; i++) { tmpBlk[i] = pVal_ij[i]; }
      MatMat<double, ndim>(pVal_ij, pVal_ii, tmpBlk);
    }
  } // end iblk
}

template<unsigned int ndim>
void BlockPentaDiagonalMatrix<ndim>::Solve(std::vector<double> &res) {
  double pTmpVec[ndim];
  for (int iblk = 0; iblk < nblk; iblk++) {
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      pTmpVec[idim] = res[iblk * ndim + idim];
    }
    if (iblk != 0) {
      const double *pVal_ij = GetValuePointer(iblk, iblk - 1);
      const double *valj = res.data() + (iblk - 1) * ndim;
      Sub_MatVec<double, ndim, ndim>(pTmpVec, pVal_ij, valj);
    }
    const double *pVal_ii = GetValuePointer(iblk, iblk);
    MatVec<double, ndim, ndim>(res.data() + iblk * ndim, pVal_ii, pTmpVec);
  }
  for (int iblk = nblk - 1; iblk >= 0; iblk--) {
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      pTmpVec[idim] = res[iblk * ndim + idim];
    }
    if (iblk != nblk - 1) {
      const double *pVal_ij = GetValuePointer(iblk, iblk + 1);
      const double *valj = res.data() + (iblk + 1) * ndim;
      Sub_MatVec<double, ndim, ndim>(res.data() + iblk * ndim, pVal_ij, valj);
    }
  }
}

}

#endif 
