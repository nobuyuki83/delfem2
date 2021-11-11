//
// Created by Nobuyuki Umetani on 2021-11-10.
//

#ifndef LS_SOLVER_BLOCK_SPARSE_H_
#define LS_SOLVER_BLOCK_SPARSE_H_

#include "delfem2/lsmats.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/vecxitrsol.h"

namespace delfem2 {

class LinearSystemSolver_BlockSparse {
 public:
  void Initialize(
      unsigned int nblk, unsigned int ndim,
      std::vector<unsigned int> &psup_ind,
      std::vector<unsigned int> &psup) {
    matrix.Initialize(nblk, ndim, true);
    matrix.SetPattern(
        psup_ind.data(), psup_ind.size(),
        psup.data(), psup.size());
    dof_bcflag.assign(ndof(), 0);
  }

  [[nodiscard]] size_t nblk() const { return matrix.nrowblk_; }
  [[nodiscard]] size_t ndim() const { return matrix.nrowdim_; }
  [[nodiscard]] size_t ndof() const { return nblk() * ndim(); }

  void BeginMerge() {
    matrix.setZero();
    vec_r.assign(ndof(), 0.0);
  }

  template<int nrow, int ncol, int ndimrow, int ndimcol>
  void Merge(
      const unsigned int *aIpRow,
      const unsigned int *aIpCol,
      const double emat[nrow][ncol][ndimrow][ndimcol]) {
    ::delfem2::Merge<nrow, ncol, ndimrow, ndimcol, double>(
        matrix, aIpRow, aIpCol, emat, merge_buffer);
  }

  void AddValueToDiagonal(unsigned int iblk, unsigned int idim, double val) {
    assert(iblk<nblk());
    const unsigned int n = ndim();
    assert(idim<n);
    matrix.val_dia_[iblk * n * n + idim * n + idim] += val;
  }

  void Solve() {
    assert(dof_bcflag.size() == ndof());
    matrix.SetFixedBC(dof_bcflag.data());
    setRHS_Zero(vec_r, dof_bcflag, 0);
    vec_x.assign(ndof(), 0.0);
    {
      tmp0.resize(ndof());
      tmp1.resize(ndof());
      conv_hist = Solve_CG(
          ViewAsVectorXd(vec_r),
          ViewAsVectorXd(vec_x),
          ViewAsVectorXd(tmp0),
          ViewAsVectorXd(tmp1),
          1.0e-4, 300, matrix);
    }
  }

 public:
  std::vector<double> vec_r;
  std::vector<double> vec_x;
  std::vector<int> dof_bcflag;
  std::vector<double> conv_hist;
  std::vector<unsigned int> merge_buffer;
  std::vector<double> tmp0, tmp1;
  CMatrixSparse<double> matrix;
};

}


#endif // LS_SOLVER_BLOCK_SPARSE_H_
