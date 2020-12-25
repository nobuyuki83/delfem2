/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femutil.h"
#include <complex>

// area of a triangle
DFM2_INLINE double delfem2::femutil::TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

// area coordinate inside a triangle
DFM2_INLINE void delfem2::femutil::TriAreaCoord(double vc_p[],
  const double p0[], const double p1[], const double p2[], const double pb[])
{
  vc_p[0] = TriArea2D(pb, p1, p2);
  vc_p[1] = TriArea2D(p0, pb, p2);
  vc_p[2] = TriArea2D(p0, p1, pb);

  const double area = TriArea2D(p0, p1, p2);
  const double inv_area = 1.0/area;

  vc_p[0] = vc_p[0]*inv_area;
  vc_p[1] = vc_p[1]*inv_area;
  vc_p[2] = vc_p[2]*inv_area;

  assert(fabs(vc_p[0]+vc_p[1]+vc_p[2]-1.0) < 1.0e-15);
}

// =======================================================================

// derivative of a shape function of a triangle and constant compornent 
DFM2_INLINE void delfem2::TriDlDx(
    double dldx[][2],
    double const_term[],
    const double p0[],
    const double p1[],
    const double p2[])
{
  const double area = ::delfem2::femutil::TriArea2D(p0, p1, p2);
  const double tmp1 = 0.5/area;

  const_term[0] = tmp1*(p1[0]*p2[1]-p2[0]*p1[1]);
  const_term[1] = tmp1*(p2[0]*p0[1]-p0[0]*p2[1]);
  const_term[2] = tmp1*(p0[0]*p1[1]-p1[0]*p0[1]);

  dldx[0][0] = tmp1*(p1[1]-p2[1]);
  dldx[1][0] = tmp1*(p2[1]-p0[1]);
  dldx[2][0] = tmp1*(p0[1]-p1[1]);

  dldx[0][1] = tmp1*(p2[0]-p1[0]);
  dldx[1][1] = tmp1*(p0[0]-p2[0]);
  dldx[2][1] = tmp1*(p1[0]-p0[0]);
}

DFM2_INLINE void delfem2::FetchData(
    double* val_to,
    int nno, int ndim,
    const unsigned int* aIP,
    const double* val_from,
    int nstride)
{
  if( nstride == -1 ){ nstride = ndim; }
  assert( nstride >= ndim );
  for(int ino=0;ino<nno;++ino){
    unsigned int ip = aIP[ino];
    for(int idim=0;idim<ndim;++idim){
      val_to[ino*ndim+idim] = val_from[ip*nstride+idim];
    }
  }
}
