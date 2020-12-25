/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>
#include <cassert>
#include <complex>

#include "delfem2/femdkt.h"

DFM2_INLINE void delfem2::MakeCurvetureDKT(
    double B1[][3], double B2[][3][2],
    const double coord0[], const double coord1[], const double coord2[],
    const double l1, const double l2 )
{
  const double l0 = 1-l1-l2;
  const double vec0[2] = { coord2[0]-coord1[0], coord2[1]-coord1[1] };
  const double vec1[2] = { coord0[0]-coord2[0], coord0[1]-coord2[1] };
  const double vec2[2] = { coord1[0]-coord0[0], coord1[1]-coord0[1] };
  const double invsqlen0 = 1.0/(vec0[0]*vec0[0]+vec0[1]*vec0[1]);
  const double invsqlen1 = 1.0/(vec1[0]*vec1[0]+vec1[1]*vec1[1]);
  const double invsqlen2 = 1.0/(vec2[0]*vec2[0]+vec2[1]*vec2[1]);
  double p0=-6*vec0[0]*invsqlen0, q0=3*vec0[0]*vec0[1]*invsqlen0, r0=3*vec0[1]*vec0[1]*invsqlen0, t0=-6*vec0[1]*invsqlen0;
  double p1=-6*vec1[0]*invsqlen1, q1=3*vec1[0]*vec1[1]*invsqlen1, r1=3*vec1[1]*vec1[1]*invsqlen1, t1=-6*vec1[1]*invsqlen1;
  double p2=-6*vec2[0]*invsqlen2, q2=3*vec2[0]*vec2[1]*invsqlen2, r2=3*vec2[1]*vec2[1]*invsqlen2, t2=-6*vec2[1]*invsqlen2;

  double H1[4][3];
  H1[0][0]=-(l1-l0)*t2+l2*t1;  H1[0][1]=+(l1-l0)*t2+l2*t0;  H1[0][2]=-l2*(t0+t1);
  H1[1][0]= (l2-l0)*t1-l1*t2;  H1[1][1]=+l1*(t0+t2);    H1[1][2]=-l1*t0-(l2-l0)*t1;
  H1[2][0]= (l1-l0)*p2-l2*p1;  H1[2][1]=-(l1-l0)*p2-l2*p0;  H1[2][2]=+l2*(p0+p1);
  H1[3][0]=-(l2-l0)*p1+l1*p2;  H1[3][1]=-l1*(p0+p2);    H1[3][2]=+l1*p0+(l2-l0)*p1;

  double H2[4][3][2];
  H2[0][0][0]=-1+(l1-l0)*r2+l2*r1;    H2[0][0][1]=-(l1-l0)*q2-l2*q1;
  H2[0][1][0]= 1+(l1-l0)*r2-l2*r0;    H2[0][1][1]=-(l1-l0)*q2+l2*q0;
  H2[0][2][0]=-l2*(r0-r1);        H2[0][2][1]= l2*(q0-q1);

  H2[1][0][0]=-1+l1*r2+(l2-l0)*r1;    H2[1][0][1]=-l1*q2-(l2-l0)*q1;
  H2[1][1][0]=-l1*(r0-r2);        H2[1][1][1]= l1*(q0-q2);
  H2[1][2][0]= 1-l1*r0+(l2-l0)*r1;    H2[1][2][1]= l1*q0-(l2-l0)*q1;

  H2[2][0][0]=-(l1-l0)*q2-l2*q1;      H2[2][0][1]= 2-6*l0-(l1-l0)*r2-l2*r1;
  H2[2][1][0]=-(l1-l0)*q2+l2*q0;      H2[2][1][1]=-2+6*l1-(l1-l0)*r2+l2*r0;
  H2[2][2][0]= l2*(q0-q1);        H2[2][2][1]= l2*(r0-r1);

  H2[3][0][0]=-l1*q2-(l2-l0)*q1;      H2[3][0][1]= 2-6*l0-l1*r2-(l2-l0)*r1;
  H2[3][1][0]= l1*(q0-q2);        H2[3][1][1]= l1*(r0-r2);
  H2[3][2][0]= l1*q0-(l2-l0)*q1;      H2[3][2][1]=-2+6*l2+l1*r0-(l2-l0)*r1;

  double dldx[3][2];
  double const_term[3];
  TriDlDx(dldx,const_term,coord0,coord1,coord2);

  for(unsigned int i=0;i<3;i++){
    B1[0][i] =  dldx[1][0]*H1[2][i]+dldx[2][0]*H1[3][i];
    B1[1][i] = -dldx[1][1]*H1[0][i]-dldx[2][1]*H1[1][i];
    B1[2][i] =  dldx[1][1]*H1[2][i]+dldx[2][1]*H1[3][i] - dldx[1][0]*H1[0][i]-dldx[2][0]*H1[1][i];
  }
  for(unsigned int i=0;i<3;i++){
    B2[0][i][0] =  dldx[1][0]*H2[2][i][0]+dldx[2][0]*H2[3][i][0];
    B2[0][i][1] =  dldx[1][0]*H2[2][i][1]+dldx[2][0]*H2[3][i][1];
    B2[1][i][0] = -dldx[1][1]*H2[0][i][0]-dldx[2][1]*H2[1][i][0];
    B2[1][i][1] = -dldx[1][1]*H2[0][i][1]-dldx[2][1]*H2[1][i][1];
    B2[2][i][0] =  dldx[1][1]*H2[2][i][0]+dldx[2][1]*H2[3][i][0] - dldx[1][0]*H2[0][i][0]-dldx[2][0]*H2[1][i][0];
    B2[2][i][1] =  dldx[1][1]*H2[2][i][1]+dldx[2][1]*H2[3][i][1] - dldx[1][0]*H2[0][i][1]-dldx[2][0]*H2[1][i][1];
  }
}

DFM2_INLINE void delfem2::MakeMat_PlateBendingDKT
(double emat_ww[3][3],
 double emat_wr[3][3][2],
 double emat_rw[3][3][2],
 double emat_rr[3][3][2][2],
 double eres_w[3],
 double eres_r[3][2],
 const double young, const double poisson, const double thickness,
 const double coord[][2], const double w[], const double rot[][2])
{
  const unsigned int ndim = 2;
  const unsigned int nno = 3;
  for(unsigned int i=0;i<nno*nno;    i++){ *(&emat_ww[0][0]      +i) = 0.0; }
  for(unsigned int i=0;i<nno*nno*2;  i++){ *(&emat_wr[0][0][0]   +i) = 0.0; }
  for(unsigned int i=0;i<nno*nno*2;  i++){ *(&emat_rw[0][0][0]   +i) = 0.0; }
  for(unsigned int i=0;i<nno*nno*2*2;i++){ *(&emat_rr[0][0][0][0]+i) = 0.0; }
  double dmat[3][3];
  {
    const double dtmp1 = young*thickness*thickness*thickness/(12.0*(1.0-poisson*poisson));
    dmat[0][0] = dtmp1;      dmat[0][1] = dtmp1*poisson;  dmat[0][2] = 0;
    dmat[1][0] = dtmp1*poisson;  dmat[1][1] = dtmp1;      dmat[1][2] = 0;
    dmat[2][0] = 0.0;      dmat[2][1] = 0.0;      dmat[2][2] = dtmp1*(1-poisson)*0.5;
  }
  double B1[3][nno];
  double B2[3][nno][ndim];
  const double area = femutil::TriArea2D(coord[0],coord[1],coord[2]);
  const double dtmp1 = area/3.0;
  for(unsigned int iw=0;iw<3;iw++){
    if(      iw == 0 ){ MakeCurvetureDKT(B1,B2,coord[0],coord[1],coord[2],0.5,0.5); }
    else if( iw == 1 ){ MakeCurvetureDKT(B1,B2,coord[0],coord[1],coord[2],0.0,0.5); }
    else if( iw == 2 ){ MakeCurvetureDKT(B1,B2,coord[0],coord[1],coord[2],0.5,0.0); }
    for(unsigned int ino=0;ino<nno;ino++){
      for(unsigned int jno=0;jno<nno;jno++){
        for(unsigned int k=0;k<3;k++){
          for(unsigned int l=0;l<3;l++){
            emat_ww[ino][jno] += dtmp1*B1[k][ino]*dmat[k][l]*B1[l][jno];
            for(unsigned int idim=0;idim<ndim;idim++){
              for(unsigned int jdim=0;jdim<ndim;jdim++){
                emat_rr[ino][jno][idim][jdim] += dtmp1*B2[k][ino][idim]*dmat[k][l]*B2[l][jno][jdim];
              }
            }
            for(unsigned int idim=0;idim<ndim;idim++){
              emat_rw[ino][jno][idim] += dtmp1*B2[k][ino][idim]*dmat[k][l]*B1[l][jno];
              emat_wr[ino][jno][idim] += dtmp1*B1[k][ino]*dmat[k][l]*B2[l][jno][idim];
            }
          }
        }
      }
    }
  }
  
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int idim=0;idim<2;idim++){
      eres_r[ino][idim] = 0.0;
      for(unsigned int jno=0;jno<nno;jno++){
        eres_r[ino][idim] -= emat_rr[ino][jno][idim][0]*rot[jno][0]
        + emat_rr[ino][jno][idim][1]*rot[jno][1]
        + emat_rw[ino][jno][idim]*w[jno];
      }
    }
  }
  for(unsigned int ino=0;ino<nno;ino++){
    eres_w[ino] = 0.0;
    for(unsigned int jno=0;jno<nno;jno++){
      eres_w[ino] -= emat_ww[ino][jno]*w[jno]
      + emat_wr[ino][jno][0]*rot[jno][0]
      + emat_wr[ino][jno][1]*rot[jno][1];
    }
  }
}