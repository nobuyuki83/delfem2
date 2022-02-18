//
// Created by Nobuyuki Umetani on 2022/02/18.
//

#ifndef DFM2_LS_SOLVER_BLOCK_SPARSE_ILU_H_
#define DFM2_LS_SOLVER_BLOCK_SPARSE_ILU_H_

#include <vector>

#include "delfem2/ls_ilu_block_sparse.h"
#include "delfem2/ls_block_sparse.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/vecxitrsol.h"

namespace delfem2 {

class LinearSystemSolver_BlockSparseILU {
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
    // initialize sparse solver
    ilu_sparse.Initialize_ILUk(matrix, 0);
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
    assert(iblk < nblk());
    const unsigned int n = ndim();
    assert(idim < n);
    matrix.val_dia_[iblk * n * n + idim * n + idim] += val;
  }

  void Solve_Cg() {
    assert(dof_bcflag.size() == ndof());
    matrix.SetFixedBC(dof_bcflag.data());
    setRHS_Zero(vec_r, dof_bcflag, 0);
    vec_x.assign(ndof(), 0.0);
    //
    tmp0.resize(ndof());
    tmp1.resize(ndof());
    conv_hist = ::delfem2::Solve_CG(
        ViewAsVectorXd(vec_r),
        ViewAsVectorXd(vec_x),
        ViewAsVectorXd(tmp0),
        ViewAsVectorXd(tmp1),
        1.0e-4, 300, matrix);
  }

  void Solve_PcgIlu() {
    namespace dfm2 = delfem2;
    matrix.SetFixedBC(dof_bcflag.data());
    dfm2::setRHS_Zero(vec_r, dof_bcflag, 0);
    //
    ilu_sparse.CopyValue(matrix);
    ilu_sparse.Decompose();
    vec_x.assign(vec_r.size(), 0.0);
    //
    tmp0.resize(ndof());
    tmp1.resize(ndof());
    std::vector<double> conv = dfm2::Solve_PCG(
        dfm2::ViewAsVectorXd(vec_r),
        dfm2::ViewAsVectorXd(vec_x),
        dfm2::ViewAsVectorXd(tmp0),
        dfm2::ViewAsVectorXd(tmp1),
        1.0e-5, 1000, matrix, ilu_sparse);
    std::cout << "convergence   nitr:" << conv.size() << "    res:" << conv[conv.size() - 1] << std::endl;
  }

 public:
  std::vector<double> vec_r;
  std::vector<double> vec_x;
  std::vector<double> tmp0, tmp1;
  std::vector<int> dof_bcflag;
  std::vector<double> conv_hist;
  std::vector<unsigned int> merge_buffer;
  CMatrixSparse<double> matrix;
  CPreconditionerILU<double> ilu_sparse;
};

}

#endif //DFM2_LS_SOLVER_BLOCK_SPARSE_ILU_H_
