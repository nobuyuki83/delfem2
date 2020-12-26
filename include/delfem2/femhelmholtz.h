/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMHELMHOLTZ_H
#define DFM2_FEMHELMHOLTZ_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/femutil.h"
#include <vector>
#include <complex>
#include <cassert>

namespace delfem2 {

void EMat_Helmholtz_Tri2D(
    std::complex<double> eres[3],
    std::complex<double> emat[][3],
    const double wave_length,
    const double coords[3][2],
    const std::complex<double> value[3]);

void EMat_SommerfeltRadiationBC_Line2D(
    std::complex<double> eres[2],
    std::complex<double> emat[2][2],
    double wave_length,
    const double P[2][2],
    const std::complex<double> val[2]);

template <class MAT>
void MergeLinSys_Helmholtz_MeshTri2D(
    MAT& mat_A,
    std::complex<double>* vec_b,
    const double wave_length,
    const double* aXY1,
    unsigned int np,
    const unsigned int* aTri1,
    unsigned int nTri,
    const std::complex<double>* aVal)
{
  using COMPLEX = std::complex<double>;
  const unsigned int nDoF = np;
  std::vector<int> tmp_buffer(nDoF, -1);
  for (unsigned int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    const COMPLEX value[3] = { aVal[i0], aVal[i1], aVal[i2] };
    //
    COMPLEX eres[3];
    COMPLEX emat[3][3];
    EMat_Helmholtz_Tri2D(eres,emat,
                         wave_length,
                         coords, value);
    for(int ino=0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

template <class MAT>
void MergeLinSys_SommerfeltRadiationBC_Polyline2D(
    MAT& mat_A,
    std::complex<double>* vec_b,
    const double wave_length,
    const double* aXY1,
    unsigned int np,
    const unsigned int* aIPPolyline,
    unsigned int nIPPolyline,
    const std::complex<double>* aVal)
{
  using COMPLEX = std::complex<double>;
  const unsigned int nDoF = np;
  std::vector<int> tmp_buffer(nDoF, -1);
  assert( nIPPolyline >= 2 );
  for(unsigned int iel=0; iel < nIPPolyline - 1; ++iel){
    const unsigned int i0 = aIPPolyline[iel + 0];
    const unsigned int i1 = aIPPolyline[iel + 1];
    const unsigned int aip[2] = {i0,i1};
    double P[2][2]; FetchData(&P[0][0],2,2,aip, aXY1);
    const COMPLEX val[2] = { aVal[i0], aVal[i1] };
    //
    COMPLEX eres[2], emat[2][2];
    EMat_SommerfeltRadiationBC_Line2D(eres,emat,
                                      wave_length,P,val);
    for(int ino=0;ino<2;ino++){
      const unsigned int ip = aip[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(2, aip, 2, aip, 1, &emat[0][0], tmp_buffer);
  }
}


} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femhelmholtz.cpp"
#endif
  
#endif /* fem_ematrix_h */
