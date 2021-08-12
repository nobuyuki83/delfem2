/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femhelmholtz.h"

#include <cmath>
#include <cassert>
#include <complex>

DFM2_INLINE void delfem2::EMat_Helmholtz_Tri2D(
    std::complex<double> eres[3],
    std::complex<double> emat[3][3],
    const double wave_length,
    const double coords[3][2],
    const std::complex<double> value[3])
{
  const int nno = 3;
  const int ndim = 2;
  const double area = ::delfem2::femutil::TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx,const_term,
          coords[0],coords[1],coords[2]);
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat[ino][jno] = area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
    }
  }
  {
    double k = 2*3.1416/wave_length;
    double tmp_val = k*k*area/12.0;
    for(unsigned int ino=0;ino<nno;ino++){
      emat[ino][ino] -= tmp_val;
      for(unsigned int jno=0;jno<nno;jno++){
        emat[ino][jno] -= tmp_val;
      }
    }
  }
  for(unsigned int ino=0;ino<nno;ino++){
    eres[ino] = 0.0;
    for(unsigned int jno=0;jno<nno;jno++){
      eres[ino] -= emat[ino][jno]*value[jno];
    }
  }
}

DFM2_INLINE void delfem2::EMat_SommerfeltRadiationBC_Line2D(
    std::complex<double> eres[2],
    std::complex<double> emat[2][2],
    double wave_length,
    const double P[2][2],
    const std::complex<double> val[2])
{
  const double elen = sqrt( (P[0][0]-P[1][0])*(P[0][0]-P[1][0]) + (P[0][1]-P[1][1])*(P[0][1]-P[1][1]) );
  {
    const double k = 2*3.1416/wave_length;
    std::complex<double> tmp_val1 = (k/6.0*elen)*std::complex<double>(0,1);
    std::complex<double> tmp_val2 = -1/(2.0*elen*k)*std::complex<double>(0,1);
    //      Com::Complex tmp_val2 = 0.0;
    emat[0][0] = tmp_val1*2.0+tmp_val2;
    emat[0][1] = tmp_val1    -tmp_val2;
    emat[1][0] = tmp_val1    -tmp_val2;
    emat[1][1] = tmp_val1*2.0+tmp_val2;
  }
  for(unsigned int ino=0;ino<2;ino++){
    eres[ino] = 0.0;
    for(unsigned int jno=0;jno<2;jno++){
      eres[ino] -= emat[ino][jno]*val[jno];
    }
  }
}
