/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/ls_masterslave.h"

#include <cassert>
#include <vector>
#include <climits>

DFM2_INLINE void delfem2::setRHS_MasterSlave(
    double* vec_b,
    unsigned int nDoF,
    const unsigned int* aMSFlag)
{
  for(unsigned int idof=0;idof<nDoF;++idof){
    unsigned int jdof = aMSFlag[idof];
    if( jdof == UINT_MAX ) continue;
    vec_b[jdof] += vec_b[idof];
    vec_b[idof] = 0;
  }
}

DFM2_INLINE void delfem2::JArray_AddMasterSlavePattern(
    std::vector<unsigned int> &index,
    std::vector<unsigned int> &array,
    const unsigned int* aMSFlag,
    unsigned int ndim,
    const unsigned int *psup_ind0,
    size_t npsup_ind0,
    const unsigned int *psup0)
{
  assert(npsup_ind0>0);
  const size_t nno = npsup_ind0-1;
  std::vector< std::vector<int> > mapM2S(nno);
  for(unsigned int ino1=0;ino1<nno;++ino1){
    for(int idim1=0;idim1<ndim;++idim1){
      unsigned int idof0 = aMSFlag[ino1*ndim+idim1];
      if( idof0 == UINT_MAX ){ continue; }
      unsigned int ino0 = idof0/ndim;
      assert( ino0 < nno && idof0 < ndim + ino0*ndim );
      assert( idim1 + ino0*ndim == idof0  );
      mapM2S[ino0].push_back(ino1);
    }
  }
  //
  index.assign(nno+1,0);
  array.clear();
  std::vector<int> aflg(nno,-1);
  //
  for(unsigned int ino0=0;ino0<nno;++ino0){
    aflg[ino0] = ino0;
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      if( aflg[jno] == ino0 ){ continue; }
      aflg[jno] = ino0;
      index[ino0+1]++;
    }
    for(int iino1=0;iino1<(int)mapM2S[ino0].size();++iino1){
      const int ino1 = mapM2S[ino0][iino1];
      if( aflg[ino1] != ino0 ){
        aflg[ino1] = ino0;
        index[ino0+1]++;
      }
      for(unsigned int jcrs=psup_ind0[ino1];jcrs<psup_ind0[ino1+1];++jcrs){
        const unsigned int jno1 = psup0[jcrs];
        if( aflg[jno1] == ino0 ){ continue; }
        aflg[jno1] = ino0;
        index[ino0+1]++;
      }
    }
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      for(int jdim=0;jdim<ndim;++jdim){
        unsigned int kdof = aMSFlag[jno*ndim+jdim];
        if( kdof == UINT_MAX ) continue;
        unsigned int kno = kdof/ndim;
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        index[ino0+1]++;
      }
    }
  }
  //
  for(unsigned int ino=0;ino<nno;ino++){ index[ino+1] += index[ino]; }
  const unsigned int narray = index[nno];
  array.resize(narray);
  for(unsigned int ino=0;ino<nno;ino++){ aflg[ino] = -1; }
  //
  for(unsigned int ino0=0;ino0<nno;++ino0){
    aflg[ino0] = ino0;
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      if( aflg[jno] == ino0 ){ continue; }
      aflg[jno] = ino0;
      const unsigned int ind = index[ino0];
      array[ind] = jno;
      index[ino0]++;
    }
    for(std::size_t jjno=0;jjno<mapM2S[ino0].size();++jjno){
      const int jno = mapM2S[ino0][jjno];
      if( aflg[jno] != ino0 ){
        aflg[jno] = ino0;
        const unsigned int ind = index[ino0];
        array[ind] = jno;
        index[ino0]++;
      }
      for(unsigned int jcrs=psup_ind0[jno];jcrs<psup_ind0[jno+1];++jcrs){
        const unsigned int kno = psup0[jcrs];
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        const unsigned int ind = index[ino0];
        array[ind] = kno;
        index[ino0]++;
      }
    }
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      for(int jdim=0;jdim<ndim;++jdim){
        unsigned int kdof = aMSFlag[jno*ndim+jdim];
        if( kdof == UINT_MAX ) continue;
        unsigned int kno = kdof/ndim;
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        const unsigned int ind = index[ino0];
        array[ind] = kno;
        index[ino0]++;
      }
    }
  }
  // ---------
  for(int ino=static_cast<int>(nno);ino>0;--ino){ index[ino] = index[ino-1]; }
  index[0] = 0;
}


DFM2_INLINE void delfem2::SetMasterSlave(
    delfem2::CMatrixSparse<double> &mat,
    const unsigned int *aMSFlag) {
  assert(!mat.val_dia_.empty());
  assert(mat.ncolblk_ == mat.nrowblk_);
  assert(mat.ncoldim_ == mat.nrowdim_);
  const unsigned int len = mat.nrowdim_;
  const unsigned int nblk = mat.nrowblk_;
  const unsigned int blksize = len * len;
  const unsigned int ndof = nblk * len;
  /////
  std::vector<unsigned int> row2crs(nblk, UINT_MAX);
  for (unsigned int idof1 = 0; idof1 < ndof; ++idof1) {  // add row
    unsigned int idof0 = aMSFlag[idof1];
    if (idof0 == UINT_MAX) continue;
    unsigned int ino0 = idof0 / len;
    unsigned int ilen0 = idof0 - ino0 * len;
    assert(ilen0 < len);
    assert(ino0 < nblk && ilen0 < len);
    unsigned int ino1 = idof1 / len;
    unsigned int ilen1 = idof1 - ino1 * len;
    assert(ino1 < nblk && ilen1 < len);
    assert(ilen0 == ilen1);
    for (unsigned int icrs0 = mat.col_ind_[ino0]; icrs0 < mat.col_ind_[ino0 + 1]; ++icrs0) {
      unsigned int jno0 = mat.row_ptr_[icrs0];
      assert(jno0 < nblk);
      row2crs[jno0] = icrs0;
    }
    for (unsigned int icrs1 = mat.col_ind_[ino1]; icrs1 < mat.col_ind_[ino1 + 1]; ++icrs1) {
      unsigned int jno1 = mat.row_ptr_[icrs1];
      assert(jno1 < nblk);
      assert(jno1 != ino1);
      if (jno1 != ino0) { // add non-diagonal 1 to non-diagonal 0
        const unsigned int icrs0 = row2crs[jno1];
        assert(icrs0 >= 0 && icrs0 < mat.row_ptr_.size());
        for (unsigned int jdim = 0; jdim < len; ++jdim) {
          mat.val_crs_[icrs0 * blksize + ilen0 * len + jdim] +=
              mat.val_crs_[icrs1 * blksize + ilen1 * len + jdim];
        }
      } else { // add non-diagonal 1 to diagonal 0
        for (unsigned int jdim = 0; jdim < len; ++jdim) {
          mat.val_dia_[ino0 * blksize + ilen0 * len + jdim] +=
              mat.val_crs_[icrs1 * blksize + ilen1 * len + jdim];
        }
      }
    }
    { // add diagonal 1 to non-diagonal 0
      const unsigned int icrs0 = row2crs[ino1];
      assert(icrs0 >= 0 && icrs0 < mat.row_ptr_.size());
      for (unsigned int jdim = 0; jdim < len; ++jdim) {
        mat.val_crs_[icrs0 * blksize + ilen0 * len + jdim]
            += mat.val_dia_[ino1 * blksize + ilen1 * len + jdim];
      }
    }
    for (unsigned int icrs0 = mat.col_ind_[ino0]; icrs0 < mat.col_ind_[ino0 + 1]; ++icrs0) {
      unsigned int jno0 = mat.row_ptr_[icrs0];
      assert(jno0 < nblk);
      row2crs[jno0] = UINT_MAX;
    }
  }
  // ---------------------------------------------
  row2crs.assign(nblk, UINT_MAX);
  for (unsigned int ino = 0; ino < nblk; ino++) {
    for (unsigned int icrs = mat.col_ind_[ino]; icrs < mat.col_ind_[ino + 1]; ++icrs) {
      unsigned int jno0 = mat.row_ptr_[icrs];
      assert(jno0 < nblk);
      row2crs[jno0] = icrs;
    }
    for (unsigned int jlen1 = 0; jlen1 < len; jlen1++) {
      unsigned int jdof0 = aMSFlag[ino * len + jlen1];
      if (jdof0 == UINT_MAX) continue;
      int jno0 = (int) (jdof0 / len);
      assert(jdof0 - jno0 * len == jlen1);
      const unsigned int icrs0 = row2crs[jno0];
      assert(icrs0 >= 0 && icrs0 < mat.row_ptr_.size());
      for (unsigned int ilen = 0; ilen < len; ilen++) {
        mat.val_crs_[icrs0 * blksize + ilen * len + jlen1] +=
            mat.val_dia_[ino * blksize + ilen * len + jlen1];
      }
    }
    for (unsigned int icrs1 = mat.col_ind_[ino]; icrs1 < mat.col_ind_[ino + 1]; icrs1++) {
      const unsigned int jno1 = mat.row_ptr_[icrs1];
      assert(jno1 < nblk);
      for (unsigned int jlen1 = 0; jlen1 < len; jlen1++) {
        if (aMSFlag[jno1 * len + jlen1] == UINT_MAX) continue;
        auto jdof0 = (unsigned int) aMSFlag[jno1 * len + jlen1];
        unsigned int jno0 = jdof0 / len;
        assert(jno0 < nblk);
        assert(jdof0 - jno0 * len == jlen1);
        if (ino == jno0) {
          for (unsigned int ilen = 0; ilen < len; ilen++) {
            mat.val_dia_[jno0 * blksize + ilen * len + jlen1] +=
                mat.val_crs_[icrs1 * blksize + ilen * len + jlen1];
          }
        } else {
          const unsigned int icrs0 = row2crs[jno0];
          assert(icrs0 >= 0 && icrs0 < mat.row_ptr_.size());
          for (unsigned int ilen = 0; ilen < len; ilen++) {
            mat.val_crs_[icrs0 * blksize + ilen * len + jlen1] +=
                mat.val_crs_[icrs1 * blksize + ilen * len + jlen1];
          }
        }
      }
    }
    for (unsigned int icrs = mat.col_ind_[ino]; icrs < mat.col_ind_[ino + 1]; ++icrs) {
      unsigned int jno0 = mat.row_ptr_[icrs];
      assert(jno0 < nblk);
      row2crs[jno0] = UINT_MAX;
    }
  }
  // --------------------------------------
  for (unsigned int iblk = 0; iblk < nblk; iblk++) {
    for (unsigned int ilen = 0; ilen < len; ilen++) {
      if (aMSFlag[iblk * len + ilen] == UINT_MAX) continue;
      for (unsigned int jlen = 0; jlen < len; jlen++) {
        mat.val_dia_[iblk * blksize + ilen * len + jlen] = 0.0;
        mat.val_dia_[iblk * blksize + jlen * len + ilen] = 0.0;
      }
      mat.val_dia_[iblk * blksize + ilen * len + ilen] = 1.0;
    }
  }
  // ---------------------------------------------
  for (unsigned int iblk = 0; iblk < nblk; iblk++) {
    for (unsigned int icrs = mat.col_ind_[iblk]; icrs < mat.col_ind_[iblk + 1]; icrs++) {
      for (unsigned int idim = 0; idim < len; idim++) {
        if (aMSFlag[iblk * len + idim] == UINT_MAX) continue;
        auto idof0 = (unsigned int) aMSFlag[iblk * len + idim];
        unsigned int jblk = mat.row_ptr_[icrs];
        for (unsigned int jdim = 0; jdim < len; jdim++) {
          unsigned int idof1 = jblk * len + jdim;
          if (idof0 != idof1) { mat.val_crs_[icrs * blksize + idim * len + jdim] = +0.0; }
          else { mat.val_crs_[icrs * blksize + idim * len + jdim] = -1.0; }
          mat.val_crs_[icrs * blksize + idim * len + jdim] = +0.0;
        }
      }
    }
  }
  // ---------------------------------------------
  for (unsigned int iblk = 0; iblk < nblk; iblk++) {
    for (unsigned int icrs = mat.col_ind_[iblk]; icrs < mat.col_ind_[iblk + 1]; icrs++) {
      const unsigned int jblk1 = mat.row_ptr_[icrs];
      for (unsigned int jdim = 0; jdim < len; jdim++) {
        if (aMSFlag[jblk1 * len + jdim] == UINT_MAX) continue;
        auto idof0 = (unsigned int) aMSFlag[jblk1 * len + jdim];
        for (unsigned int idim = 0; idim < len; idim++) {
          unsigned int idof1 = iblk * len + idim;
          if (idof0 != idof1) { mat.val_crs_[icrs * blksize + idim * len + jdim] = +0.0; }
          else { mat.val_crs_[icrs * blksize + idim * len + jdim] = -1.0; }
          mat.val_crs_[icrs * blksize + idim * len + jdim] = +0.0;
        }
      }
    }
  }
}
