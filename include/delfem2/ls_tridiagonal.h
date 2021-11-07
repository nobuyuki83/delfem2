#ifndef DFM2_LS_TRIDIAGONAL_H
#define DFM2_LS_TRIDIAGONAL_H

#include "delfem2/matn.hpp"

namespace delfem2 {

/**
 * 0 1 #  #  #
 * 2 3 4  #  #
 * # 5 6  7  #
 * # # 8  9 10
 * # # # 11 12
 */
template<unsigned int ndim>
class BlockTriDiagonalMatrix {
 public:
  // initialize with block size n
  BlockTriDiagonalMatrix() { nblk = 0; }
  ~BlockTriDiagonalMatrix() {}

  void Initialize(size_t n) {
    this->nblk = n;
    v.resize((nblk * 3 - 2) * ndim * ndim );
  }

  void Clear() {
    nblk = 0;
    v.clear();
  }

  double *GetValuePointer(int iblk, int jblk) {
    constexpr unsigned int blksize = ndim * ndim;
    if (iblk < 0 || iblk >= nblk) { return nullptr; }
    if (jblk < 0 || jblk >= nblk) { return nullptr; }
    if (iblk - jblk < -1 || iblk - jblk > +1) { return nullptr; }
    return v.data() + (iblk * 3 + jblk - iblk) * blksize;
  }

  //! clear value
  void SetZero() { v.assign((nblk * 3 - 2) * ndim * ndim, 0.0); }

  //! marge element stiffness matrix to the position (idiv,idiv+1)
  void Merge(unsigned int idiv, double eM[][2][ndim][ndim]) {
    for (unsigned int i = 0; i < 4*ndim*ndim; i++) {
      v[idiv * ndim*ndim*3 + i] += (&eM[0][0][0][0])[i];
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
void BlockTriDiagonalMatrix<ndim>::FixBC(
    unsigned int iblk,
    unsigned int idim) {
  assert(idim < ndim && int(iblk) < nblk);
  const int iblk0 = static_cast<int>(iblk);
  for (int jblk = iblk0 - 1; jblk <= iblk0 + 1; ++jblk) {
    double *p = this->GetValuePointer(iblk, jblk);
    if (p == nullptr) { continue; }
    for (unsigned int jdim = 0; jdim < ndim; ++jdim) {
      p[idim * ndim + jdim] = 0;
    }
  }
  for (int jblk = iblk0 - 1; jblk <= iblk0 + 1; ++jblk) {
    double *p = this->GetValuePointer(jblk, iblk);
    if (p == nullptr) { continue; }
    for (unsigned int jdim = 0; jdim < ndim; ++jdim) {
      p[jdim * ndim + idim] = 0;
    }
  }
  GetValuePointer(iblk, iblk)[idim * ndim + idim] = 1;
}

template<unsigned int ndim>
void BlockTriDiagonalMatrix<ndim>::Decompose() {
  constexpr unsigned int blksize = ndim * ndim;
  double tmpBlk[blksize];
  for (int iblk = 0; iblk < nblk; iblk++) {
    if (iblk != 0) {
      const double *pVal_ik = GetValuePointer(iblk, iblk - 1);
      const double *pVal_kj = GetValuePointer(iblk - 1, iblk);
      double *pVal_ij = GetValuePointer(iblk, iblk);
      Sub_MatMat<double,ndim>(pVal_ij, pVal_ik, pVal_kj);
    }
    {   // calc inverse of diagonal
      double *pVal_ii = GetValuePointer(iblk, iblk);
      Inverse_Matrix<double,ndim>(pVal_ii);
    }
    // [U] = [1/D][U]
    if (iblk != nblk - 1) {
      double *pVal_ij = GetValuePointer(iblk, iblk + 1);
      const double *pVal_ii = GetValuePointer(iblk, iblk);
      for (unsigned int i = 0; i < blksize; i++) { tmpBlk[i] = pVal_ij[i]; }
      MatMat<double,ndim>(pVal_ij, pVal_ii, tmpBlk);
    }
  } // end iblk
}

template<unsigned int ndim>
void BlockTriDiagonalMatrix<ndim>::Solve(std::vector<double> &res) {
  double pTmpVec[ndim];
  for (int iblk = 0; iblk < nblk; iblk++) {
    for(unsigned int idim=0;idim<ndim;++idim){
      pTmpVec[idim] = res[iblk * ndim + idim];
    }
    if (iblk != 0) {
      const double *pVal_ij = GetValuePointer(iblk,iblk-1);
      const double *valj = res.data() + (iblk - 1) * ndim;
      Sub_MatVec<double,ndim,ndim>(pTmpVec, pVal_ij, valj);
    }
    const double *pVal_ii = GetValuePointer(iblk,iblk);
    MatVec<double,ndim,ndim>(res.data() + iblk * ndim, pVal_ii, pTmpVec);
  }
  for (int iblk = nblk - 1; iblk >= 0; iblk--) {
    for(unsigned int idim=0;idim<ndim;++idim){
      pTmpVec[idim] = res[iblk * ndim + idim];
    }
    if (iblk != nblk - 1) {
      const double *pVal_ij = GetValuePointer(iblk,iblk+1);
      const double *valj = res.data() + (iblk + 1) * ndim;
      Sub_MatVec<double,ndim,ndim>(res.data() + iblk * ndim, pVal_ij, valj);
    }
  }
}

}

#endif 
