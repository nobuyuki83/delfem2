/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femutil.h"
#include <complex>

// area of a triangle
DFM2_INLINE double delfem2::femutil::TriArea2D(
    const double p0[],
    const double p1[],
    const double p2[])
{
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

// area coordinate inside a triangle
DFM2_INLINE void delfem2::femutil::TriAreaCoord(
    double vc_p[],
    const double p0[],
    const double p1[],
    const double p2[],
    const double pb[])
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

DFM2_INLINE double delfem2::femutil::Dot3D(const double a[], const double b[]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

DFM2_INLINE void delfem2::femutil::Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
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


DFM2_INLINE double delfem2::femutil::TetVolume3D(
    const double v1[3],
    const double v2[3],
    const double v3[3],
    const double v4[3])
{
  return
      ((v2[0]-v1[0])*((v3[1]-v1[1])*(v4[2]-v1[2])-(v4[1]-v1[1])*(v3[2]-v1[2]))
       -(v2[1]-v1[1])*((v3[0]-v1[0])*(v4[2]-v1[2])-(v4[0]-v1[0])*(v3[2]-v1[2]))
       +(v2[2]-v1[2])*((v3[0]-v1[0])*(v4[1]-v1[1])-(v4[0]-v1[0])*(v3[1]-v1[1]))
      ) * 0.16666666666666666666666666666667;
}

// caluculate Derivative of Area Coord
DFM2_INLINE void delfem2::TetDlDx(
    double dldx[][3],
    double a[],
    const double p0[],
    const double p1[],
    const double p2[],
    const double p3[])
{
  const double vol = femutil::TetVolume3D(p0, p1, p2, p3);
  const double dtmp1 = 1.0/(vol * 6.0);

  a[0] = +dtmp1*(p1[0]*(p2[1]*p3[2]-p3[1]*p2[2])-p1[1]*(p2[0]*p3[2]-p3[0]*p2[2])+p1[2]*(p2[0]*p3[1]-p3[0]*p2[1]));
  a[1] = -dtmp1*(p2[0]*(p3[1]*p0[2]-p0[1]*p3[2])-p2[1]*(p3[0]*p0[2]-p0[0]*p3[2])+p2[2]*(p3[0]*p0[1]-p0[0]*p3[1]));
  a[2] = +dtmp1*(p3[0]*(p0[1]*p1[2]-p1[1]*p0[2])-p3[1]*(p0[0]*p1[2]-p1[0]*p0[2])+p3[2]*(p0[0]*p1[1]-p1[0]*p0[1]));
  a[3] = -dtmp1*(p0[0]*(p1[1]*p2[2]-p2[1]*p1[2])-p0[1]*(p1[0]*p2[2]-p2[0]*p1[2])+p0[2]*(p1[0]*p2[1]-p2[0]*p1[1]));

  dldx[0][0] = -dtmp1*((p2[1]-p1[1])*(p3[2]-p1[2])-(p3[1]-p1[1])*(p2[2]-p1[2]));
  dldx[0][1] = +dtmp1*((p2[0]-p1[0])*(p3[2]-p1[2])-(p3[0]-p1[0])*(p2[2]-p1[2]));
  dldx[0][2] = -dtmp1*((p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]));

  dldx[1][0] = +dtmp1*((p3[1]-p2[1])*(p0[2]-p2[2])-(p0[1]-p2[1])*(p3[2]-p2[2]));
  dldx[1][1] = -dtmp1*((p3[0]-p2[0])*(p0[2]-p2[2])-(p0[0]-p2[0])*(p3[2]-p2[2]));
  dldx[1][2] = +dtmp1*((p3[0]-p2[0])*(p0[1]-p2[1])-(p0[0]-p2[0])*(p3[1]-p2[1]));

  dldx[2][0] = -dtmp1*((p0[1]-p3[1])*(p1[2]-p3[2])-(p1[1]-p3[1])*(p0[2]-p3[2]));
  dldx[2][1] = +dtmp1*((p0[0]-p3[0])*(p1[2]-p3[2])-(p1[0]-p3[0])*(p0[2]-p3[2]));
  dldx[2][2] = -dtmp1*((p0[0]-p3[0])*(p1[1]-p3[1])-(p1[0]-p3[0])*(p0[1]-p3[1]));

  dldx[3][0] = +dtmp1*((p1[1]-p0[1])*(p2[2]-p0[2])-(p2[1]-p0[1])*(p1[2]-p0[2]));
  dldx[3][1] = -dtmp1*((p1[0]-p0[0])*(p2[2]-p0[2])-(p2[0]-p0[0])*(p1[2]-p0[2]));
  dldx[3][2] = +dtmp1*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));

  //  std::cout << dldx[0][0]+dldx[1][0]+dldx[2][0]+dldx[3][0] << std::endl;
  //  std::cout << dldx[0][1]+dldx[1][1]+dldx[2][1]+dldx[3][1] << std::endl;
  //  std::cout << dldx[0][2]+dldx[1][2]+dldx[2][2]+dldx[3][2] << std::endl;

  //  std::cout << a[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1]+dldx[0][2]*p0[2] << std::endl;
  //  std::cout << a[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1]+dldx[1][2]*p1[2] << std::endl;
  //  std::cout << a[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1]+dldx[2][2]*p2[2] << std::endl;
  //  std::cout << a[3]+dldx[3][0]*p3[0]+dldx[3][1]*p3[1]+dldx[3][2]*p3[2] << std::endl;
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
