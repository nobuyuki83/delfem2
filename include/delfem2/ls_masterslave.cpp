/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <vector>
#include <climits>
#include "delfem2/ls_masterslave.h"

DFM2_INLINE void
delfem2::setRHS_MasterSlave(
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

DFM2_INLINE void
delfem2::JArray_AddMasterSlavePattern(
    std::vector<unsigned int> &index,
    std::vector<unsigned int> &array,
    const unsigned int* aMSFlag,
    int ndim,
    const unsigned int *psup_ind0,
    int npsup_ind0,
    const unsigned int *psup0)
{
  assert(npsup_ind0>0);
  const int nno = npsup_ind0-1;
  //assert( aMSFlag.size() == nno*ndim );
  std::vector< std::vector<int> > mapM2S(nno);
  for(int ino1=0;ino1<nno;++ino1){
    for(int idim1=0;idim1<ndim;++idim1){
      int idof0 = aMSFlag[ino1*ndim+idim1];
      if( idof0 == -1 ){ continue; }
      int ino0 = idof0/ndim;
//      int idim0 = idof0 - ino0*ndim;
      assert( ino0 < nno && idof0 - ino0*ndim < ndim );
//      std::cout << idim1 << " " << idim0 << " " << ndim << std::endl;
      assert( idim1 == idof0 - ino0*ndim );
      mapM2S[ino0].push_back(ino1);
    }
  }
  //
  index.assign(nno+1,0);
  array.clear();
  std::vector<int> aflg(nno,-1);
  ///
  for(int ino0=0;ino0<nno;++ino0){
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
        int kdof = aMSFlag[jno*ndim+jdim];
        if( kdof == -1 ) continue;
        int kno = kdof/ndim;
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        index[ino0+1]++;
      }
    }
  }
  //
  for(int ino=0;ino<nno;ino++){ index[ino+1] += index[ino]; }
  const int narray = index[nno];
  array.resize(narray);
  for(int ino=0;ino<nno;ino++){ aflg[ino] = -1; }
  //
  for(int ino0=0;ino0<nno;++ino0){
    aflg[ino0] = ino0;
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      if( aflg[jno] == ino0 ){ continue; }
      aflg[jno] = ino0;
      const int ind = index[ino0];
      array[ind] = jno;
      index[ino0]++;
    }
    for(std::size_t jjno=0;jjno<mapM2S[ino0].size();++jjno){
      const int jno = mapM2S[ino0][jjno];
      if( aflg[jno] != ino0 ){
        aflg[jno] = ino0;
        const int ind = index[ino0];
        array[ind] = jno;
        index[ino0]++;
      }
      for(unsigned int jcrs=psup_ind0[jno];jcrs<psup_ind0[jno+1];++jcrs){
        const unsigned int kno = psup0[jcrs];
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        const int ind = index[ino0];
        array[ind] = kno;
        index[ino0]++;
      }
    }
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      for(int jdim=0;jdim<ndim;++jdim){
        int kdof = aMSFlag[jno*ndim+jdim];
        if( kdof == -1 ) continue;
        int kno = kdof/ndim;
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        const int ind = index[ino0];
        array[ind] = kno;
        index[ino0]++;
      }
    }
  }
  // ---------
  for(int ino=nno;ino>0;ino--){ index[ino] = index[ino-1]; }
  index[0] = 0;
}
