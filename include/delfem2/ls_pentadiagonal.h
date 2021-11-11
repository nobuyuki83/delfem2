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
  BlockPentaDiagonalMatrix() = default;
  ~BlockPentaDiagonalMatrix() {}

  unsigned int nblk() const { return nblk_; }

  void Initialize(size_t n) {
    this->nblk_ = n;
    assert(nblk_ >= 4);
    const unsigned int nblks = 3 + 4 + 5 * (nblk_ - 4) + 4 + 3;
    v.resize(nblks * ndim * ndim);
  }

  void Clear() {
    nblk_ = 0;
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
    if (iblk < 0 || iblk >= nblk_) { return nullptr; }
    if (jblk < 0 || jblk >= nblk_) { return nullptr; }
    if (iblk - jblk < -2 || iblk - jblk > +2) { return nullptr; }
    if (iblk < 2) { return v.data() + (iblk * 3 + jblk) * blksize; }
    if (iblk == nblk_ - 1) { return v.data() + (iblk * 4 - 2 + jblk) * blksize; }
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
    assert(iblk < nblk_ - 2);
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
  void DecompIJK(
      int iblk,
      int jblk,
      int kblk) {
    if (iblk < 0 || jblk < 0 || kblk < 0) { return; }
    assert(kblk < iblk && kblk < jblk);
    const double *pVal_ik = GetValuePointer(iblk, kblk);
    if (pVal_ik == nullptr) { return; }
    const double *pVal_kj = GetValuePointer(kblk, jblk);
    if (pVal_kj == nullptr) { return; }
    double *pVal_ij = GetValuePointer(iblk, jblk);
    if (pVal_ij == nullptr) { return; }
    Sub_MatMat<double, ndim>(pVal_ij, pVal_ik, pVal_kj);
  }
  void UpdateU(int iblk, int jblk) {
    assert(jblk > iblk);
    double *pVal_ij = GetValuePointer(iblk, jblk);
    if (pVal_ij == nullptr) { return; }
    constexpr unsigned int blksize = ndim * ndim;
    const double *pVal_ii = GetValuePointer(iblk, iblk);
    assert(pVal_ii != nullptr);
    double tmpBlk[blksize];
    for (unsigned int i = 0; i < blksize; i++) { tmpBlk[i] = pVal_ij[i]; }
    MatMat<double, ndim>(pVal_ij, pVal_ii, tmpBlk);
  }
  void Substitution(double *pTmpVec, const double *res, int iblk, int jblk) {
    const double *pVal_ij = GetValuePointer(iblk, jblk);
    if (pVal_ij == nullptr) { return; }
    Sub_MatVec<double, ndim, ndim>(pTmpVec, pVal_ij, res);
  }
 private:
  int nblk_ = 0;
  std::vector<double> v;
};

template<unsigned int ndim_>
class LinearSystemSolver_BlockPentaDiagonal {
 public:
  LinearSystemSolver_BlockPentaDiagonal() = default;

  void Initialize(
      unsigned int nblk) {
    this->dia.Initialize(nblk);
    dof_bcflag.assign(ndof(), 0);
  }

  [[nodiscard]] size_t nblk() const { return dia.nblk(); }
  [[nodiscard]] size_t ndim() const { return ndim_; }
  [[nodiscard]] size_t ndof() const { return nblk() * ndim(); }

  void BeginMerge() {
    dia.setZero();
    vec_r.assign(ndof(), 0.);
  }

  template<int nrow, int ncol, int ndimrow, int ndimcol>
  void Merge(
      const unsigned int *aIpRow,
      const unsigned int *aIpCol,
      const double emat[nrow][ncol][ndimrow][ndimcol]) {
    assert(ndimrow <= ndim_ && ndimcol <= ndim_);
    for (unsigned int irow = 0; irow < nrow; ++irow) {
      for (unsigned int icol = 0; icol < ncol; ++icol) {
        assert(aIpRow[irow] < dia.nblk() && aIpCol[icol] < dia.nblk());
        double *p = dia.GetValuePointer(aIpRow[irow], aIpCol[icol]);
        assert(p != nullptr);
        for (int i = 0; i < ndimrow; ++i) {
          for (int j = 0; j < ndimrow; ++j) {
            p[i * ndim_ + j] += emat[irow][icol][i][j];
          }
        }
      }
    }
  }

  void AddValueToDiagonal(
      unsigned int ip,
      unsigned int idim,
      double val) {
    double *p = dia.GetValuePointer(ip, ip);
    p[idim * ndim_ + idim] += val;
  }

  void Solve() {
    for (unsigned int iblk = 0; iblk < nblk(); ++iblk) {
      for (unsigned int idim = 0; idim < ndim_; ++idim) {
        if (dof_bcflag[iblk * ndim_ + idim] == 0) { continue; }
        dia.FixBC(iblk, idim);
        vec_r[iblk * ndim_ + idim] = 0.;
      }
    }
    // BlockPentaDiagonalMatrix<ndim_> dia0 = dia;
    dia.Decompose();
    vec_x = vec_r;
    dia.Solve(vec_x);
    /*
    for(int iblk=0;iblk<dia.nblk();++iblk){
      double tmp[ndim_];
      std::fill_n(tmp,ndim_,0.);
      for(int jblk=iblk-2;jblk<=iblk+2;++jblk) {
        double *p = dia0.GetValuePointer(iblk, jblk);
        if( p == nullptr ){ continue; }
        Sub_MatVec<double, ndim_, ndim_>(tmp, p, vec_x.data() + jblk*ndim_);
      }
      for(unsigned int idim=0;idim<ndim_;++idim) {
        std::cout << iblk << " " << idim << " " << tmp[idim] << " " << tmp[idim] + vec_r[iblk * ndim_ + idim] << std::endl;
      }
    }
     */
  }

 public:
  std::vector<double> vec_r;
  std::vector<double> vec_x;
  std::vector<int> dof_bcflag;
  BlockPentaDiagonalMatrix<ndim_> dia;
};

}  // delfem2


// =======================================================

//! define fixed boudnary condition
template<unsigned int ndim>
void delfem2::BlockPentaDiagonalMatrix<ndim>::FixBC(
    unsigned int iblk,
    unsigned int idim) {
  assert(idim < ndim && int(iblk) < nblk_);
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
void delfem2::BlockPentaDiagonalMatrix<ndim>::Decompose() {
  // constexpr unsigned int blksize = ndim * ndim;
  for (int iblk = 0; iblk < nblk_; iblk++) {
    DecompIJK(iblk,iblk-1,iblk-2);
    DecompIJK(iblk,iblk,iblk-2);
    DecompIJK(iblk,iblk,iblk-1);
    DecompIJK(iblk,iblk+1,iblk-1);
    {   // calc inverse of diagonal
      double *pVal_ii = GetValuePointer(iblk, iblk);
      Inverse_Matrix<double, ndim>(pVal_ii);
    }
    // [U] = [1/D][U]
    UpdateU(iblk,iblk+1);
    UpdateU(iblk,iblk+2);
  } // end iblk
}

template<unsigned int ndim>
void delfem2::BlockPentaDiagonalMatrix<ndim>::Solve(std::vector<double> &res) {
  double pTmpVec[ndim];
  for (int iblk = 0; iblk < nblk_; iblk++) {
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      pTmpVec[idim] = res[iblk * ndim + idim];
    }
    Substitution(pTmpVec, res.data() + (iblk - 2) * ndim, iblk, iblk-2);
    Substitution(pTmpVec, res.data() + (iblk - 1) * ndim, iblk, iblk-1);
    {
      const double *pVal_ii = GetValuePointer(iblk, iblk);
      MatVec<double, ndim, ndim>(res.data() + iblk * ndim, pVal_ii, pTmpVec);
    }
  }
  for (int iblk = nblk_ - 2; iblk >= 0; --iblk) {
    Substitution(res.data() + iblk * ndim,res.data() + (iblk + 1) * ndim, iblk, iblk+1);
    Substitution(res.data() + iblk * ndim,res.data() + (iblk + 2) * ndim, iblk, iblk+2);
  }
}

#endif 
