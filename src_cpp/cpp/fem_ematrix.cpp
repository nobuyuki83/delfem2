/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <math.h>
#include <assert.h>
#include <iostream>

#include "delfem2/fem_ematrix.h"

/////////////////////////////////////////////////////////////////////////////////

const static unsigned int NIntLineGauss[4] = {
  1, 2, 3, 4
};
const static double LineGauss[4][4][2] =
{
  {
    { 0.0, 2.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
  },
  {
    { -0.577350269189626, 1.0 },
    {  0.577350269189626, 1.0 },
    {  0.0,               0.0 },
    {  0.0,               0.0 },
  },
  {
    { -0.774596669241483, 0.555555555555556 },
    {  0.0,               0.888888888888889 },
    {  0.774596669241483, 0.555555555555556 },
    {  0.0,               0.0               },
  },
  {
    { -0.861136311594053, 0.347854845137454 },
    { -0.339981043584856, 0.652145154862546 },
    {  0.339981043584856, 0.652145154862546 },
    {  0.861136311594053, 0.347854845137454 },
  }
};



void ShapeFunc_Hex8
(const double& r0, const double& r1,	const double& r2,
 const double coords[][3],
 double& detjac,
 double dndx[][3],
 double an[] )
{
  an[0] = 0.125*(1.0-r0)*(1.0-r1)*(1.0-r2);
  an[1] = 0.125*(1.0+r0)*(1.0-r1)*(1.0-r2);
  an[2] = 0.125*(1.0-r0)*(1.0+r1)*(1.0-r2);
  an[3] = 0.125*(1.0+r0)*(1.0+r1)*(1.0-r2);
  an[4] = 0.125*(1.0-r0)*(1.0-r1)*(1.0+r2);
  an[5] = 0.125*(1.0+r0)*(1.0-r1)*(1.0+r2);
  an[6] = 0.125*(1.0-r0)*(1.0+r1)*(1.0+r2);
  an[7] = 0.125*(1.0+r0)*(1.0+r1)*(1.0+r2);
  
  double dndr[8][3];
  dndr[0][0] = -0.125*(1.0-r1)*(1.0-r2);
  dndr[1][0] = -dndr[0][0];
  dndr[2][0] = -0.125*(1.0+r1)*(1.0-r2);
  dndr[3][0] = -dndr[2][0];
  dndr[4][0] = -0.125*(1.0-r1)*(1.0+r2);
  dndr[5][0] = -dndr[4][0];
  dndr[6][0] = -0.125*(1.0+r1)*(1.0+r2);
  dndr[7][0] = -dndr[6][0];
  
  dndr[0][1] = -0.125*(1.0-r0)*(1.0-r2);
  dndr[1][1] = -0.125*(1.0+r0)*(1.0-r2);
  dndr[2][1] = -dndr[0][1];
  dndr[3][1] = -dndr[1][1];
  dndr[4][1] = -0.125*(1.0-r0)*(1.0+r2);
  dndr[5][1] = -0.125*(1.0+r0)*(1.0+r2);
  dndr[6][1] = -dndr[4][1];
  dndr[7][1] = -dndr[5][1];
  
  dndr[0][2] = -0.125*(1.0-r0)*(1.0-r1);
  dndr[1][2] = -0.125*(1.0+r0)*(1.0-r1);
  dndr[2][2] = -0.125*(1.0-r0)*(1.0+r1);
  dndr[3][2] = -0.125*(1.0+r0)*(1.0+r1);
  dndr[4][2] = -dndr[0][2];
  dndr[5][2] = -dndr[1][2];
  dndr[6][2] = -dndr[2][2];
  dndr[7][2] = -dndr[3][2];
  
  double dxdr[3][3]  = {
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
  };
  
  for(int inode=0;inode<8;inode++){
    dxdr[0][0] += coords[inode][0]*dndr[inode][0];
    dxdr[0][1] += coords[inode][0]*dndr[inode][1];
    dxdr[0][2] += coords[inode][0]*dndr[inode][2];
    dxdr[1][0] += coords[inode][1]*dndr[inode][0];
    dxdr[1][1] += coords[inode][1]*dndr[inode][1];
    dxdr[1][2] += coords[inode][1]*dndr[inode][2];
    dxdr[2][0] += coords[inode][2]*dndr[inode][0];
    dxdr[2][1] += coords[inode][2]*dndr[inode][1];
    dxdr[2][2] += coords[inode][2]*dndr[inode][2];
  }
  
  detjac = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]
  + dxdr[1][0]*dxdr[2][1]*dxdr[0][2]
  + dxdr[2][0]*dxdr[0][1]*dxdr[1][2]
  - dxdr[0][0]*dxdr[2][1]*dxdr[1][2]
  - dxdr[1][0]*dxdr[0][1]*dxdr[2][2]
  - dxdr[2][0]*dxdr[1][1]*dxdr[0][2];
  
  const double inv_jac = 1.0 / detjac;
  
  double drdx[3][3];
  drdx[0][0] = inv_jac*( dxdr[1][1]*dxdr[2][2]-dxdr[1][2]*dxdr[2][1] );
  drdx[0][1] = inv_jac*( dxdr[0][2]*dxdr[2][1]-dxdr[0][1]*dxdr[2][2] );
  drdx[0][2] = inv_jac*( dxdr[0][1]*dxdr[1][2]-dxdr[0][2]*dxdr[1][1] );
  drdx[1][0] = inv_jac*( dxdr[1][2]*dxdr[2][0]-dxdr[1][0]*dxdr[2][2] );
  drdx[1][1] = inv_jac*( dxdr[0][0]*dxdr[2][2]-dxdr[0][2]*dxdr[2][0] );
  drdx[1][2] = inv_jac*( dxdr[0][2]*dxdr[1][0]-dxdr[0][0]*dxdr[1][2] );
  drdx[2][0] = inv_jac*( dxdr[1][0]*dxdr[2][1]-dxdr[1][1]*dxdr[2][0] );
  drdx[2][1] = inv_jac*( dxdr[0][1]*dxdr[2][0]-dxdr[0][0]*dxdr[2][1] );
  drdx[2][2] = inv_jac*( dxdr[0][0]*dxdr[1][1]-dxdr[0][1]*dxdr[1][0] );
  
  for(int inode=0;inode<8;inode++){
    dndx[inode][0] = dndr[inode][0]*drdx[0][0] + dndr[inode][1]*drdx[1][0] + dndr[inode][2]*drdx[2][0];
    dndx[inode][1] = dndr[inode][0]*drdx[0][1] + dndr[inode][1]*drdx[1][1] + dndr[inode][2]*drdx[2][1];
    dndx[inode][2] = dndr[inode][0]*drdx[0][2] + dndr[inode][1]*drdx[1][2] + dndr[inode][2]*drdx[2][2];
  }
}


// area of a triangle
static double TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

static double TriArea3D(const double v1[3], const double v2[3], const double v3[3]){
  double x, y, z;
  x = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  y = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return 0.5*sqrt( x*x + y*y + z*z );
}

static double Distance3D(const double p0[3], const double p1[3]){
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

static double Dot3D(const double a[], const double b[]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

static void Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

static void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

// area coordinate inside a triangle
void TriAreaCoord(double vc_p[],
  const double p0[], const double p1[], const double p2[], const double pb[]){

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

// derivative of a shape function of a triangle and constant compornent 
void TriDlDx(double dldx[][2], double const_term[],
  const double p0[], const double p1[], const double p2[])
{
  const double area = TriArea2D(p0, p1, p2);
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
  /*
  assert( fabs( dldx[0][0]+dldx[1][0]+dldx[2][0] ) < 1.0e-15 );
  assert( fabs( dldx[0][1]+dldx[1][1]+dldx[2][1] ) < 1.0e-15 );

  assert( fabs( const_term[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1] - 1.0 ) < 1.0e-10 );
  assert( fabs( const_term[0]+dldx[0][0]*p1[0]+dldx[0][1]*p1[1] ) < 1.0e-10 );
  assert( fabs( const_term[0]+dldx[0][0]*p2[0]+dldx[0][1]*p2[1] ) < 1.0e-10 );

  assert( fabs( const_term[1]+dldx[1][0]*p0[0]+dldx[1][1]*p0[1] ) < 1.0e-10 );
  assert( fabs( const_term[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1] - 1.0 ) < 1.0e-10 );
  assert( fabs( const_term[1]+dldx[1][0]*p2[0]+dldx[1][1]*p2[1] ) < 1.0e-10 );

  assert( fabs( const_term[2]+dldx[2][0]*p0[0]+dldx[2][1]*p0[1] ) < 1.0e-10 );
  assert( fabs( const_term[2]+dldx[2][0]*p1[0]+dldx[2][1]*p1[1] ) < 1.0e-10 );
  assert( fabs( const_term[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1] - 1.0 ) < 1.0e-10 );
  */
}

/*
const static unsigned int NIntTriGauss[3] = { 1, 3, 7 };
const static double TriGauss[3][7][3] =
{
  { // liner
    { 0.3333333333, 0.3333333333, 1.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
  },
  { // quadratic
    { 0.1666666667, 0.1666666667, 0.3333333333 },
    { 0.6666666667, 0.1666666667, 0.3333333333 },
    { 0.1666666667, 0.6666666667, 0.3333333333 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
  },
  { // cubic
    { 0.1012865073, 0.1012865073, 0.1259391805 },
    { 0.7974269854, 0.1012865073, 0.1259391805 },
    { 0.1012865073, 0.7974269854, 0.1259391805 },
    { 0.4701420641, 0.0597158718, 0.1323941527 },
    { 0.4701420641, 0.4701420641, 0.1323941527 },
    { 0.0597158718, 0.4701420641, 0.1323941527 },
    { 0.3333333333, 0.3333333333, 0.225 },
  }
};
 */

void MakeMat_Poisson2D_P1
(const double alpha, const double source,
 const double coords[3][2],
 const double value[3],
 double eres[3],
 double emat[][3])
{
  const int nno = 3;
  const int ndim = 2;
  ////
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;
  for (int i = 0; i<9; ++i){ (&emat[0][0])[i] = 0.0; }
  const double area = TriArea2D(coords[0], coords[1], coords[2]);
  ////
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx, const_term, coords[0], coords[1], coords[2]);
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      emat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
    }
  }
  for (int ino = 0; ino<nno; ino++){
    eres[ino] = source*area*0.33333333333333333;
  }
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      eres[ino] -= emat[ino][jno]*value[jno];
    }
  }
}


void MakeMat_Diffusion2D_P1
(const double alpha, const double source,
 const double dt_timestep, const double gamma_newmark, const double rho,
 const double coords[3][2],
 const double value[3], const double velo[3],
 double eres[3],
 double emat[3][3])
{
  const int nno = 3;
  const int ndim = 2;
  
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;
  for (int i = 0; i<9; ++i){ (&emat[0][0])[i] = 0.0; }
  
  const double area = TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx,const_term,coords[0],coords[1],coords[2]);
  
  ////////////////////////////////////////////////////////////////
  
  double eCmat[nno][nno];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
    }
  }
  double eMmat[nno][nno];
  {
    const double dtmp1 = rho*area*0.08333333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno] = dtmp1;
      }
      eMmat[ino][ino] += dtmp1;
    }
  }
		
  for(int ino=0;ino<nno;ino++){
    eres[ino] = source*area*0.333333333333333333;
  }
  
  ////////////////////////////////////////////////////////////////
  
  {
    const double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno;i++){
      (&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
    }
  }
	for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino]	-= eCmat[ino][jno]*(value[jno]+dt_timestep*velo[jno]) + eMmat[ino][jno]*velo[jno];
    }
  }
}

void MakeMat_LinearSolid2D_Static_P1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y,
 const double disp[3][2],
 const double coords[3][2],
 double eres[3][2],
 double emat[3][3][2][2])
{
  const int nno = 3;
  const int ndim = 2;
  
  const double area = TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim], zero_order_term[nno];
  TriDlDx(dldx, zero_order_term, coords[0],coords[1],coords[2]);
  
  for(int ino=0;ino<nno;ino++){
     for(int jno=0;jno<nno;jno++){
       emat[ino][jno][0][0] = area*(lambda+myu)*dldx[ino][0]*dldx[jno][0];
       emat[ino][jno][0][1] = area*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
       emat[ino][jno][1][0] = area*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
       emat[ino][jno][1][1] = area*(lambda+myu)*dldx[ino][1]*dldx[jno][1];
       const double dtmp1 = (dldx[ino][1]*dldx[jno][1]+dldx[ino][0]*dldx[jno][0])*area*myu;
       emat[ino][jno][0][0] += dtmp1;
       emat[ino][jno][1][1] += dtmp1;
     }
  }
  
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = area*rho*g_x*0.33333333333333333;
    eres[ino][1] = area*rho*g_y*0.33333333333333333;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1];
      eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1];
    }
  }
}


void MakeMat_LinearSolid2D_Dynamic_P1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y,
 const double dt_timestep, const double gamma_newmark,  const double beta_newmark,
 const double disp[3][2], const double velo[3][2], const double acc[3][2],
 const double coords[3][2],
 double eres[3][2],
 double emat[3][3][2][2],
 bool is_initial)
{
  const int nno = 3;
  const int ndim = 2;
  
	const double area = TriArea2D(coords[0],coords[1],coords[2]);
	double dldx[nno][ndim];		// spatial derivative of linear shape function
	{
    double zero_order_term[nno];	// const term of shape function
    TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);
  }
  
  double eKmat[nno][nno][ndim][ndim];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eKmat[ino][jno][0][0] = area*(lambda+myu)*dldx[ino][0]*dldx[jno][0];
      eKmat[ino][jno][0][1] = area*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
      eKmat[ino][jno][1][0] = area*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
      eKmat[ino][jno][1][1] = area*(lambda+myu)*dldx[ino][1]*dldx[jno][1];
      const double dtmp1 = (dldx[ino][1]*dldx[jno][1]+dldx[ino][0]*dldx[jno][0])*area*myu;
      eKmat[ino][jno][0][0] += dtmp1;
      eKmat[ino][jno][1][1] += dtmp1;
    }
  }
  
  double eMmat[nno][nno][ndim][ndim];
  {
    const double dtmp1 = area*rho*0.0833333333333333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno][0][0] = dtmp1;
        eMmat[ino][jno][0][1] = 0.0;
        eMmat[ino][jno][1][0] = 0.0;
        eMmat[ino][jno][1][1] = dtmp1;
      }
      eMmat[ino][ino][0][0] += dtmp1;
      eMmat[ino][ino][1][1] += dtmp1;
    }
  }
  
	// calc external force
	for(int ino=0;ino<nno;ino++){
    eres[ino][0] = area*rho*g_x*0.33333333333333333333333333;
    eres[ino][1] = area*rho*g_y*0.33333333333333333333333333;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  {	// calc coeff matrix for newmark-beta
    double dtmp1 = beta_newmark*dt_timestep*dt_timestep;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
    }
  }
  
  // calc element redisual vector
	for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eKmat[ino][jno][0][0]*disp[jno][0]+eKmat[ino][jno][0][1]*disp[jno][1];
      eres[ino][1] -= eKmat[ino][jno][1][0]*disp[jno][0]+eKmat[ino][jno][1][1]*disp[jno][1];
    }
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eMmat[ino][jno][0][0]*acc[jno][0]+eMmat[ino][jno][0][1]*acc[jno][1];
      eres[ino][1] -= eMmat[ino][jno][1][0]*acc[jno][0]+eMmat[ino][jno][1][1]*acc[jno][1];
    }
  }
	if( is_initial ){
    for(int ino=0;ino<nno;ino++){
     for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= dt_timestep*(eKmat[ino][jno][0][0]*velo[jno][0]+eKmat[ino][jno][0][1]*velo[jno][1]);
        eres[ino][1] -= dt_timestep*(eKmat[ino][jno][1][0]*velo[jno][0]+eKmat[ino][jno][1][1]*velo[jno][1]);
      }
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= 0.5*dt_timestep*dt_timestep*(eKmat[ino][jno][0][0]*acc[jno][0]+eKmat[ino][jno][0][1]*acc[jno][1]);
        eres[ino][1] -= 0.5*dt_timestep*dt_timestep*(eKmat[ino][jno][1][0]*acc[jno][0]+eKmat[ino][jno][1][1]*acc[jno][1]);
      }
    }
  }
}


void MakeMat_Stokes2D_Static_P1P1
(double alpha, double g_x, double g_y,
 const double coords[][2],
 const double velo[3][2], const double press[3],
 double emat_uu[][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3],
 double eres_u[][2], double eres_p[3])
{
  const unsigned int nno = 3;
  const unsigned int ndim = 2;
  
  const double area = TriArea2D(coords[0],coords[1],coords[2]);
  
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  
  ////////////////////////////////////////////////////////////
  
  for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat_uu[0][0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      const double dtmp1 = area*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
      emat_uu[ino][jno][0][0] = dtmp1;
      emat_uu[ino][jno][1][1] = dtmp1;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_up[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_up[ino][jno][0] += area*dldx[ino][0]*0.333333333333333333333333;
      emat_up[ino][jno][1] += area*dldx[ino][1]*0.333333333333333333333333;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_pu[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pu[ino][jno][0] += area*dldx[jno][0]*0.333333333333333333333333;
      emat_pu[ino][jno][1] += area*dldx[jno][1]*0.333333333333333333333333;
    }
  }
  for(unsigned int i=0;i<nno*nno;i++){ *(&emat_pp[0][0]+i) = 0.0; }
  double tau; // relaxation parameter
  {
    const double h = sqrt( area / 3.14 )*2;
    tau = -h*h/alpha*0.1;
//    tau = 0.0;
  }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pp[ino][jno] = area*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
    }
  }
  
  for(unsigned int ino=0;ino<nno;ino++){
    eres_u[ino][0] = area*g_x*0.3333333333333333333;
    eres_u[ino][1] = area*g_y*0.3333333333333333333;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_uu[ino][jno][0][0]*velo[jno][0]+emat_uu[ino][jno][0][1]*velo[jno][1];
      eres_u[ino][1] -= emat_uu[ino][jno][1][0]*velo[jno][0]+emat_uu[ino][jno][1][1]*velo[jno][1];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_up[ino][jno][0]*press[jno];
      eres_u[ino][1] -= emat_up[ino][jno][1]*press[jno];
    }
  }
  for(unsigned int ino=0;ino<nno;ino++){
    eres_p[ino] = 0.0;
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pu[ino][jno][0]*velo[jno][0]+emat_pu[ino][jno][1]*velo[jno][1];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pp[ino][jno]*press[jno];
    }
  }
}

void MakeMat_Stokes2D_Static_P1
(double alpha, double g_x, double g_y,
 const double coords[][2],
 const double velo_press[3][3],
 double emat[3][3][3][3],
 double eres[3][3])
{
  const int nno = 3;
  
  const double velo[3][2] = {
    {velo_press[0][0],velo_press[0][1]},
    {velo_press[1][0],velo_press[1][1]},
    {velo_press[2][0],velo_press[2][1]} };
  const double press[3] = { velo_press[0][2], velo_press[1][2], velo_press[2][2] };
  ////
  double emat_uu[3][3][2][2], emat_up[3][3][2], emat_pu[3][3][2], emat_pp[3][3];
  double eres_u[3][2], eres_p[3];
  MakeMat_Stokes2D_Static_P1P1(alpha, g_x, g_y,
                               coords, velo, press,
                               emat_uu, emat_up, emat_pu, emat_pp,
                               eres_u, eres_p);
  ////
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
      emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
      emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
      emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
      emat[ino][jno][0][2] = emat_up[ino][jno][0];
      emat[ino][jno][1][2] = emat_up[ino][jno][1];
      emat[ino][jno][2][0] = emat_pu[ino][jno][0];
      emat[ino][jno][2][1] = emat_pu[ino][jno][1];
      emat[ino][jno][2][2] = emat_pp[ino][jno];
    }
  }
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = eres_u[ino][0];
    eres[ino][1] = eres_u[ino][1];
    eres[ino][2] = eres_p[ino];
  }
}


void MakeMat_Stokes2D_Dynamic_Newmark_P1P1
(double alpha, double rho, double g_x, double g_y,
 const double dt_timestep, const double gamma_newmark,
 const double coords[3][2],
 const double velo[3][2], const double press[3], const double acc[3][2], const double apress[3],
 double emat_uu[3][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3],
 double eres_u[3][2], double eres_p[3])
{
  //	std::cout << "AddMat_Stokes2D_NonStatic_Newmark_P1P1" << std::endl;
  
  const int nno = 3;
  const int ndim = 2;
  
  ////
  const double area = TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);

  //////////////////////////////////////////////////////////////////////////////
  double eCmat_uu[3][3][2][2];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = area*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
      eCmat_uu[ino][jno][0][0] = dtmp1;
      eCmat_uu[ino][jno][0][1] = 0.0;
      eCmat_uu[ino][jno][1][0] = 0.0;
      eCmat_uu[ino][jno][1][1] = dtmp1;
    }
  }
  
  double eCmat_up[3][3][2];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] = area*dldx[ino][0]*0.333333333333333333333;
      eCmat_up[ino][jno][1] = area*dldx[ino][1]*0.333333333333333333333;
    }
  }
  
  double eCmat_pu[3][3][2];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pu[ino][jno][0] = area*dldx[jno][0]*0.333333333333333333333;
      eCmat_pu[ino][jno][1] = area*dldx[jno][1]*0.333333333333333333333;
    }
  }
  
  double tau;
  {
    const double h = sqrt( area / 3.14 )*2;
    tau = -h*h/alpha*0.1;
  }
  
  double eCmat_pp[3][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pp[ino][jno] = area*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
    }
  }
  
  ////////////////
  double eMmat_uu[3][3][2][2];
  {
    const double dtmp1 = area*rho*0.0833333333333333333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat_uu[ino][jno][0][0] = dtmp1;
        eMmat_uu[ino][jno][0][1] = 0.0;
        eMmat_uu[ino][jno][1][0] = 0.0;
        eMmat_uu[ino][jno][1][1] = dtmp1;
      }
      eMmat_uu[ino][ino][0][0] += dtmp1;
      eMmat_uu[ino][ino][1][1] += dtmp1;
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    eres_u[ino][0] = area*g_x*0.33333333333333333333;
    eres_u[ino][1] = area*g_y*0.33333333333333333333;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  
  ////////////////////////////////
  
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          for(int jdim=0;jdim<ndim;jdim++){
            emat_uu[ino][jno][idim][jdim] = eMmat_uu[ino][jno][idim][jdim]+dtmp1*eCmat_uu[ino][jno][idim][jdim];
          }
        }
      }
    }
    
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          emat_up[ino][jno][idim] = dtmp1*eCmat_up[ino][jno][idim];
          emat_pu[ino][jno][idim] = dtmp1*eCmat_pu[ino][jno][idim];
        }
      }
    }
      
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        emat_pp[ino][jno] = dtmp1*eCmat_pp[ino][jno];
      }
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<ndim;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<ndim;jdim++){
          eres_u[ino][idim] -= eCmat_uu[ino][jno][idim][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim])
          + eMmat_uu[ino][jno][idim][jdim]*acc[jno][jdim];
        }
      }
      for(int jno=0;jno<nno;jno++){
        eres_u[ino][idim] -= eCmat_up[ino][jno][idim]*(press[jno]+dt_timestep*apress[jno]);
      }
    }
  }

  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      for(int jdim=0;jdim<ndim;jdim++){
        eres_p[ino] -= eCmat_pu[ino][jno][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim]);
      }
    }
    for(int jno=0;jno<nno;jno++){
      eres_p[ino] -= eCmat_pp[ino][jno]*(press[jno]+dt_timestep*apress[jno]);
    }
  }
}

void MakeMat_Stokes2D_Dynamic_P1
(double alpha, double rho, double g_x, double g_y,
 const double dt_timestep, const double gamma_newmark,
 const double coords[][2],
 const double velo_press[3][3], const double acc_apress[3][3],
 double emat[3][3][3][3],
 double eres[3][3])
{
  const int nno = 3;

  const double area = TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][2], const_term[nno];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  
  double tau;
  {
    const double h = sqrt( area / 3.14 )*2;
    tau = -h*h/alpha*0.1;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  double eCmat[3][3][3][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = area*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
      eCmat[ino][jno][0][0] = dtmp1;
      eCmat[ino][jno][1][1] = dtmp1;
      eCmat[ino][jno][0][1] = 0.0;
      eCmat[ino][jno][1][0] = 0.0;
      eCmat[ino][jno][2][0] = area*dldx[jno][0]*0.333333333333333333333;
      eCmat[ino][jno][2][1] = area*dldx[jno][1]*0.333333333333333333333;
      eCmat[ino][jno][0][2] = area*dldx[ino][0]*0.333333333333333333333;
      eCmat[ino][jno][1][2] = area*dldx[ino][1]*0.333333333333333333333;
      eCmat[ino][jno][2][2] = area*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
    }
  }
  
  ////////////////
  double eMmat[3][3][3][3];
  {
    const double dtmp1 = area*rho*0.0833333333333333333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno][0][0] = dtmp1;
        eMmat[ino][jno][1][1] = dtmp1;
        eMmat[ino][jno][0][1] = 0.0;
        eMmat[ino][jno][1][0] = 0.0;
        eMmat[ino][jno][0][2] = 0.0;
        eMmat[ino][jno][1][2] = 0.0;
        eMmat[ino][jno][2][0] = 0.0;
        eMmat[ino][jno][2][1] = 0.0;
        eMmat[ino][jno][2][2] = 0.0;
      }
      eMmat[ino][ino][0][0] += dtmp1;
      eMmat[ino][ino][1][1] += dtmp1;
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = area*g_x*0.33333333333333333333;
    eres[ino][1] = area*g_y*0.33333333333333333333;
    eres[ino][2] = 0.0;
  }
  
  ////////////////////////////////
  
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<3*3*3*3;++i){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eCmat[0][0][0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<3;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<3;jdim++){
          eres[ino][idim]
          -= eCmat[ino][jno][idim][jdim]*(velo_press[jno][jdim]+dt_timestep*acc_apress[jno][jdim])
          + eMmat[ino][jno][idim][jdim]*acc_apress[jno][jdim];
        }
      }
    }
  }
}

void MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1
(double rho, double myu, double g_x, double g_y,
 double dt, double gamma,
 const double coords[][2],
 const double velo[][2], const double press[], const double acc[][2], const double apress[],
 double emat_uu[][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3],
 double eres_u[3][2], double eres_p[3])
{
  const int nno = 3;
  const int ndim = 2;

	const double area = TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[3][2], const_term[3];
	TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  const double velo_ave[2] = {
    (velo[0][0]+velo[1][0]+velo[2][0])/3.0,
    (velo[0][1]+velo[1][1]+velo[2][1])/3.0 };

  double tau;
  {  // Calc Stabilization Parameter

    const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
    const double h = sqrt( area / 3.14 )*2;
    const double tau_c = h*0.5/norm_v;
    const double cou_c = norm_v*dt/h;
    if( norm_v*h*rho*1.0e-30 > myu ){ // Re = \infty
      const double dtmp1 = 1/(cou_c*cou_c)+1;
      tau = tau_c / sqrt(dtmp1);
    }
    else if( norm_v*h*rho < myu*1.0e-30 ){ // Re = 0
      tau = h*h*rho*0.5/myu;
    }
    else{
      const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
      const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
      tau = tau_c / sqrt(dtmp1);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  
  double eCmat_uu[3][3][2][2];
	// viscousity
	for(int ino=0;ino<nno;ino++){
	for(int jno=0;jno<nno;jno++){
		eCmat_uu[ino][jno][0][0] = area*myu*dldx[ino][0]*dldx[jno][0];
		eCmat_uu[ino][jno][0][1] = area*myu*dldx[ino][1]*dldx[jno][0];
		eCmat_uu[ino][jno][1][0] = area*myu*dldx[ino][0]*dldx[jno][1];
		eCmat_uu[ino][jno][1][1] = area*myu*dldx[ino][1]*dldx[jno][1];
		const double dtmp1 = area*myu*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		eCmat_uu[ino][jno][0][0] += dtmp1;
		eCmat_uu[ino][jno][1][1] += dtmp1;
	}
	}
	{	// advection
		const double dtmp0[2] = { velo[0][0]+velo[1][0]+velo[2][0], velo[0][1]+velo[1][1]+velo[2][1] };
		for(int jno=0;jno<nno;jno++){
			const double dtmp1 = (dldx[jno][0]*dtmp0[0]+dldx[jno][1]*dtmp0[1]);
			for(int ino=0;ino<nno;ino++){
				double dtmp2 = dtmp1 + (dldx[jno][0]*velo[ino][0]+dldx[jno][1]*velo[ino][1]);
				dtmp2 *= area*rho*0.083333333333333;
				eCmat_uu[ino][jno][0][0] += dtmp2;
				eCmat_uu[ino][jno][1][1] += dtmp2;
			}
		}
	}
  {	// SUPG for advection term
    double tmp_mat[ndim][ndim] = { {0,0}, {0,0} };
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        tmp_mat[0][0] += velo[ino][0]*velo[jno][0];
        tmp_mat[0][1] += velo[ino][0]*velo[jno][1];
        tmp_mat[1][0] += velo[ino][1]*velo[jno][0];
        tmp_mat[1][1] += velo[ino][1]*velo[jno][1];
      }
      tmp_mat[0][0] += velo[ino][0]*velo[ino][0];
      tmp_mat[1][1] += velo[ino][1]*velo[ino][1];
    }
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        double dtmp1 = 0.0;
        dtmp1 += dldx[ino][0]*dldx[jno][0]*tmp_mat[0][0]
        +dldx[ino][0]*dldx[jno][1]*tmp_mat[0][1]
        +dldx[ino][1]*dldx[jno][0]*tmp_mat[1][0]
        +dldx[ino][1]*dldx[jno][1]*tmp_mat[1][1];
        dtmp1 *= tau*rho*area*0.083333333333333333333;
        eCmat_uu[ino][jno][0][0] += dtmp1;
        eCmat_uu[ino][jno][1][1] += dtmp1;
      }
    }
  }

  double eCmat_up[3][3][2], eCmat_pu[3][3][2];
	{   // add pressure gradient, non compression term
		const double dtmp1 = area*0.33333333333333333;
		for(int ino=0;ino<nno;ino++){
		for(int jno=0;jno<nno;jno++){
			eCmat_up[ino][jno][0] = -dtmp1*dldx[ino][0];
			eCmat_up[ino][jno][1] = -dtmp1*dldx[ino][1];
			eCmat_pu[ino][jno][0] =  dtmp1*dldx[jno][0];
			eCmat_pu[ino][jno][1] =  dtmp1*dldx[jno][1];
		}
		}
	}
  // SUPG for pressure gradient
  for(int ino=0;ino<nno;ino++){
    double dtmp1 = (dldx[ino][0]*velo_ave[0]+dldx[ino][1]*velo_ave[1])*tau*area;
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] += dtmp1*dldx[jno][0];
      eCmat_up[ino][jno][1] += dtmp1*dldx[jno][1];
    }
  }
  // PSPG for advection term
  for(int jno=0;jno<nno;jno++){
    const double dtmp1 = (dldx[jno][0]*velo_ave[0]+dldx[jno][1]*velo_ave[1])*tau*area;
    for(int ino=0;ino<nno;ino++){
      eCmat_pu[ino][jno][0] += dtmp1*dldx[ino][0];
      eCmat_pu[ino][jno][1] += dtmp1*dldx[ino][1];
    }
  }
  
  double eCmat_pp[3][3];
  // PSPG for pressure term
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pp[ino][jno] = area*tau/rho*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
    }
  }

	////////////////
  double eMmat_uu[3][3][2][2];
	{	// add inertia term
		const double dtmp1 = area*rho*0.0833333333333333;
		for(int ino=0;ino<nno;ino++){
			for(int jno=0;jno<nno;jno++){
				eMmat_uu[ino][jno][0][0] = dtmp1;
				eMmat_uu[ino][jno][0][1] = 0.0;
        eMmat_uu[ino][jno][1][0] = 0.0;
        eMmat_uu[ino][jno][1][1] = dtmp1;
			}
			eMmat_uu[ino][ino][0][0] += dtmp1;
			eMmat_uu[ino][ino][1][1] += dtmp1;
		}
	}
  // supg for inertia term
  for(int jno=0;jno<nno;jno++){
    double tmp_vec[ndim] = { 0.0, 0.0 };
    for(unsigned int kno=0;kno<nno;kno++){
      tmp_vec[0] += velo[kno][0];
      tmp_vec[1] += velo[kno][1];
    }
    tmp_vec[0] += velo[jno][0];
    tmp_vec[1] += velo[jno][1];
    for(int ino=0;ino<nno;ino++){
      const double dtmp1 = (dldx[ino][0]*tmp_vec[0]+dldx[ino][1]*tmp_vec[1])*rho*tau*area*0.083333333333333;
      eMmat_uu[ino][jno][0][0] += dtmp1;
      eMmat_uu[ino][jno][1][1] += dtmp1;
    }
  }

  // PSPG for inertia term
  double eMmat_pu[3][3][2];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eMmat_pu[ino][jno][0] = -tau*area*0.33333333333333333*dldx[ino][0];
      eMmat_pu[ino][jno][1] = -tau*area*0.33333333333333333*dldx[ino][1];
    }
  }
  
  // eternal force term
	for(int ino=0;ino<nno;ino++){
		eres_u[ino][0] = area*rho*g_x*0.3333333333333333333;
		eres_u[ino][1] = area*rho*g_y*0.3333333333333333333;
	}
	for(int ino=0;ino<nno;ino++){
		eres_p[ino] = 0.0;
	}

/*		
		// SUPG external force
		for(unsigned int ino=0;ino<nno;ino++){
			const double dtmp1 = area*tau*rho*(ave_velo[0]*dldx[ino][0]+ave_velo[1]*dldx[ino][1]);
			eres_u[ino][0] += dtmp1*g_x;
			eres_u[ino][1] += dtmp1*g_y;
		}
*/

/*
		// PSPG external force
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] += area*tau*(dldx[ino][0]*g_x+dldx[ino][1]*g_y);
		}
*/
  //////////////////////////////////////////////////////////////////////
  
  {
    double dtmp1 = gamma*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
    }
    for(int i=0;i<nno*nno;i++){
      (&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
      +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
      +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
      eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
      +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
      +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
    }
  }
	for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
      +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
      +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
    }
  }
}

void MakeMat_NavierStokes2D_Dynamic_P1
(double myu, double rho, double g_x, double g_y,
 const double dt_timestep, const double gamma_newmark,
 const double coords[][2],
 const double velo_press[3][3], const double acc_apress[3][3],
 double emat[3][3][3][3],
 double eres[3][3])
{
  const int nno = 3;
  
  const double velo[3][2] = {
    {velo_press[0][0],velo_press[0][1]},
    {velo_press[1][0],velo_press[1][1]},
    {velo_press[2][0],velo_press[2][1]} };
  const double press[3] = { velo_press[0][2], velo_press[1][2], velo_press[2][2] };
  const double acc[3][2] = {
    {acc_apress[0][0],acc_apress[0][1]},
    {acc_apress[1][0],acc_apress[1][1]},
    {acc_apress[2][0],acc_apress[2][1]} };
  const double apress[3] = { acc_apress[0][2], acc_apress[1][2], acc_apress[2][2] };
  ////
  double emat_uu[3][3][2][2], emat_up[3][3][2], emat_pu[3][3][2], emat_pp[3][3];
  double eres_u[3][2], eres_p[3];
  MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1(rho, myu, g_x, g_y,
                                              dt_timestep, gamma_newmark,
                                              coords,
                                              velo, press, acc, apress,
                                              emat_uu, emat_up, emat_pu, emat_pp,
                                              eres_u, eres_p);
  ////
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
      emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
      emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
      emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
      emat[ino][jno][0][2] = emat_up[ino][jno][0];
      emat[ino][jno][1][2] = emat_up[ino][jno][1];
      emat[ino][jno][2][0] = emat_pu[ino][jno][0];
      emat[ino][jno][2][1] = emat_pu[ino][jno][1];
      emat[ino][jno][2][2] = emat_pp[ino][jno];
    }
  }
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = eres_u[ino][0];
    eres[ino][1] = eres_u[ino][1];
    eres[ino][2] = eres_p[ino];
  }
}


/*
static bool AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1(
		double gamma, double dt, CLinearSystem_Field& ls, 
		double rho, double myu, double g_x, double g_y, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		const unsigned int id_ea )
{
//    std::cout << "NavierStorkes2D_NonStatic_Newmark Triangle 3-point 1st order " << gamma << " " << dt << " " << rho << " " << myu << " " << id_ea << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg es_velo_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg es_velo_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);
	const CElemAry::CElemSeg es_pres_c_va = field_press.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim], eMmat_pu[nno][nno][ndim];
	double  emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eres_u[nno][ndim], eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
        !=  ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo, CORNER,  world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo,CORNER, id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,CORNER, id_field_velo,CORNER,  world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,  world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,CORNER,  world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);//na_velo.GetSeg(id_ns_acc);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);//na_press.GetSeg(id_ns_apress);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no_v[nno];	// 要素節点の全体節点番号				
		es_velo_c_co.GetNodes(ielem,no_v);	// 要素の節点番号を取ってくる
		double coords[nno][ndim];	// 要素節点の座標
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(no_v[ino],coords[ino]);
		}
		// 要素の節点番号を取ってくる
		es_velo_c_va.GetNodes(ielem,no_v);
		double velo[nno][ndim];	// 要素節点の値
		double acc[nno][ndim];	// 要素節点の値
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_velo.GetValue(no_v[ino],velo[ino]);
			ns_acc.GetValue(no_v[ino],acc[ino]);
		}

		unsigned int no_p[nno];	// 要素節点の全体節点番号				
		es_pres_c_va.GetNodes(ielem,no_p);	// 要素の節点番号を取ってくる
		double press[nno];
		double apress[nno];		
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_press.GetValue(no_p[ino],&press[ino]);
			ns_apress.GetValue(no_p[ino],&apress[ino]);
		}

        ////////////////////////////////

        MakeMat_NavierStokes2D_NonStatic_Newmark_P1P1(dt,gamma,
            rho,myu,g_x,g_y,
            coords,   velo,acc,   press,apress,
            eres_u,eres_p, 
            eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp, eMmat_uu,eMmat_pu);

        ////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
			eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
			              +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
			              +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
		}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
			}
		}

		// 要素剛性行列の全体剛性行列へのマージ
		mat_uu.Mearge(nno,no_v,nno,no_v,	4,&emat_uu[0][0][0][0]);
		mat_up.Mearge(nno,no_v,nno,no_p,	2,&emat_up[0][0][0]   );
		mat_pu.Mearge(nno,no_p,nno,no_v,	2,&emat_pu[0][0][0]   );
		mat_pp.Mearge(nno,no_p,nno,no_p,	1,&emat_pp[0][0]      );
		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( no_v[ino],0,eres_u[ino][0]);
			res_u.AddValue( no_v[ino],1,eres_u[ino][1]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_p.AddValue( no_p[ino],0,eres_p[ino]);
		}
	}
	return true;
}


static bool AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1_Combined(
		double gamma, double dt, CLinearSystem_Field& ls, 
		double rho, double myu, double g_x, double g_y, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		const unsigned int id_ea )
{
//    std::cout << "NavierStorkes2D_NonStatic_Newmark Triangle 3-point 1st order Combined" << gamma << " " << dt << " " << rho << " " << myu << " " << id_ea << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg es_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim], eMmat_pu[nno][nno][ndim];
	double  emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eres_u[nno][ndim], eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		 == ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(  id_field_velo, CORNER,world);
	CVector_Blk&    res_u  = ls.GetResidual(id_field_velo, CORNER,world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);//na_velo.GetSeg(id_ns_acc);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);//na_press.GetSeg(id_ns_apress);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int noes[nno];	// 要素節点の全体節点番号
		
		// 要素の節点番号を取ってくる
		es_c_co.GetNodes(ielem,noes);
		double coords[nno][ndim];	// 要素節点の座標
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
		}
        // 要素の節点番号を取ってくる(id_esは流速と圧力で同じはず)
		es_c_va.GetNodes(ielem,noes);
		double velo[nno][ndim];	// 要素節点の値
		double acc[nno][ndim];	// 要素節点の値
		double press[nno];
		double apress[nno];		
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_velo.GetValue(noes[ino],velo[ino]);
			ns_acc.GetValue(noes[ino],acc[ino]);
			ns_press.GetValue(noes[ino],&press[ino]);
			ns_apress.GetValue(noes[ino],&apress[ino]);
		}


        MakeMat_NavierStokes2D_NonStatic_Newmark_P1P1(dt,gamma,
            rho,myu,g_x,g_y,
            coords,   velo,acc,   press,apress,
            eres_u,eres_p, 
            eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp, eMmat_uu,eMmat_pu);


		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
			eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
			              +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
			              +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
		}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
			}
		}

		double emat[nno][nno][3][3];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
			emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
			emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
			emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
			emat[ino][jno][0][2] = emat_up[ino][jno][0];
			emat[ino][jno][1][2] = emat_up[ino][jno][1];
			emat[ino][jno][2][0] = emat_pu[ino][jno][0];
			emat[ino][jno][2][1] = emat_pu[ino][jno][1];
			emat[ino][jno][2][2] = emat_pp[ino][jno];
		}
		}
		mat_uu.Mearge(nno,noes,nno,noes,	9,&emat[0][0][0][0]);
		// merge residual vector
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( noes[ino],0,eres_u[ino][0]);
			res_u.AddValue( noes[ino],1,eres_u[ino][1]);
			res_u.AddValue( noes[ino],2,eres_p[ino]);
		}
    }
	return true;
}


bool Fem::Eqn::AddLinSys_NavierStokes2D_NonStatic_Newmark(
		double dt, double gamma, CLinearSystem_Field& ls,
		double rho, double alpha, double g_x, double g_y,
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( field_velo.GetFieldType() != VECTOR2 ) return false;
	if( field_press.GetFieldType() != SCALAR ) return false;

	if( id_ea != 0 ){
		bool res;
		if( field_velo.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_velo,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_press,CORNER,world) )
            {
			    res = AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1_Combined(
				    gamma,dt,  ls,
				    rho,alpha, g_x,g_y,
				    id_field_velo,id_field_press,world,
                    id_ea);
            }
            else{
			    res = AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1(
				    gamma,dt,  ls,
				    rho,alpha, g_x,g_y,
				    id_field_velo,id_field_press,world,
                    id_ea);
            }
		}
		else{
			assert(0); 
			res = false;
		}
		return res;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_velo.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			// 再帰文
			bool res = AddLinSys_NavierStokes2D_NonStatic_Newmark(
				dt, gamma, ls,
				rho, alpha, g_x, g_y,
				id_field_velo, id_field_press, world,
				id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}

*/

/*
static bool AddLinearSystem_Poisson2D_P1(
  double alpha, double source,
  Fem::Ls::CLinearSystem_Field& ls,
  const unsigned int id_field_val, const CFieldWorld& world,
  const unsigned int id_ea)
{
  //	std::cout << "Poisson2D Triangle 3-point 1st order" << std::endl;

  assert(world.IsIdEA(id_ea));
  const CElemAry& ea = world.GetEA(id_ea);
  assert(ea.ElemType()==TRI);

  if (!world.IsIdField(id_field_val)) return false;
  const CField& field_val = world.GetField(id_field_val);

  const CElemAry::CElemSeg& es_c_va = field_val.GetElemSeg(id_ea, CORNER, true, world);
  const CElemAry::CElemSeg& es_c_co = field_val.GetElemSeg(id_ea, CORNER, false, world);

  const unsigned int nno = 3;
  const unsigned int ndim = 2;

  unsigned int no_c[nno];	// 要素節点の全体節点番号

  double value_c[nno];		// 要素節点の値
  double coord_c[nno][ndim];	// 要素節点の座標

  double emat[nno][nno];	// 要素剛性行列
  double eres_c[nno];	// 要素節点等価内力、外力、残差ベクトル

  CMatDia_BlkCrs& mat_cc = ls.GetMatrix(id_field_val, CORNER, world);
  CVector_Blk&    res_c = ls.GetResidual(id_field_val, CORNER, world);

  const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER, true, world);
  const CNodeAry::CNodeSeg& ns_c_co = field_val.GetNodeSeg(CORNER, false, world);

  for (unsigned int ielem = 0; ielem<ea.Size(); ielem++){
    // 要素配列から要素セグメントの節点番号を取り出す
    es_c_co.GetNodes(ielem, no_c);
    for (unsigned int inoes = 0; inoes<nno; inoes++){
      ns_c_co.GetValue(no_c[inoes], coord_c[inoes]);
    }
    // 節点の値を取って来る
    es_c_va.GetNodes(ielem, no_c);
    for (unsigned int inoes = 0; inoes<nno; inoes++){
      ns_c_val.GetValue(no_c[inoes], &value_c[inoes]);
    }
    // 面積を求める
    const double area = TriArea3D(coord_c[0], coord_c[1], coord_c[2]);
    // 形状関数の微分を求める
    double dldx[nno][ndim];	// 形状関数のxy微分
    double const_term[nno];	// 形状関数の定数項
    TriDlDx(dldx, const_term, coord_c[0], coord_c[1], coord_c[2]);
    // 要素剛性行列を作る
    for (unsigned int ino = 0; ino<nno; ino++){
      for (unsigned int jno = 0; jno<nno; jno++){
        emat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
      }
    }
    // 要素節点等価外力ベクトルを求める
    for (unsigned int ino = 0; ino<nno; ino++){
      eres_c[ino] = source*area*0.33333333333333333;
    }
    // 要素節点等価内力ベクトルを求める
    for (unsigned int ino = 0; ino<nno; ino++){
      for (unsigned int jno = 0; jno<nno; jno++){
        eres_c[ino] -= emat[ino][jno]*value_c[jno];
      }
    }
    // 要素剛性行列にマージする
    mat_cc.Mearge(nno, no_c, nno, no_c, 1, &emat[0][0]);
    // 残差ベクトルにマージする
    for (unsigned int inoes = 0; inoes<nno; inoes++){
      res_c.AddValue(no_c[inoes], 0, eres_c[inoes]);
    }
  }
  return true;
}
*/

// compute energy and its 1st and 2nd derivative for cloth bending
void WdWddW_Bend
(double& W,  // (out) strain energy
 double dW[4][3], // (out) 1st derivative of energy
 double ddW[4][4][3][3], // (out) 2nd derivative of energy
 ////
 const double C[4][3], // (in) undeformed triangle vertex positions
 const double c[4][3], // (in) deformed triangle vertex positions
 double stiff)
{
  const double A0 = TriArea3D(C[0],C[2],C[3]);
  const double A1 = TriArea3D(C[1],C[3],C[2]);
  const double L0 = Distance3D(C[2],C[3]);
  const double H0 = A0*2.0/L0;
  const double H1 = A1*2.0/L0;
  const double e23[3] = { C[3][0]-C[2][0], C[3][1]-C[2][1], C[3][2]-C[2][2] };
  const double e02[3] = { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] };
  const double e03[3] = { C[3][0]-C[0][0], C[3][1]-C[0][1], C[3][2]-C[0][2] };
  const double e12[3] = { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] };
  const double e13[3] = { C[3][0]-C[1][0], C[3][1]-C[1][1], C[3][2]-C[1][2] };
  double cot023, cot032;
  {
    const double r2 = -Dot3D(e02,e23);
    const double r3 = +Dot3D(e03,e23);
    cot023 = r2/H0;
    cot032 = r3/H0;
  }
  double cot123, cot132;
  {
    const double r2 = -Dot3D(e12,e23);
    const double r3 = +Dot3D(e13,e23);
    cot123 = r2/H1;
    cot132 = r3/H1;
  }
  const double tmp0 = stiff/((A0+A1)*L0*L0);
  const double K[4] = { -cot023-cot032, -cot123-cot132, cot032+cot132, cot023+cot123 };
  
  // compute 2nd derivative of energy
  for(int i=0;i<4*4*3*3;i++){ (&ddW[0][0][0][0])[i] = 0; }
  for(int ino=0;ino<4;ino++){
    for(int jno=0;jno<4;jno++){
      const double tmp = K[ino]*K[jno]*tmp0;
      ddW[ino][jno][0][0] = tmp;
      ddW[ino][jno][1][1] = tmp;
      ddW[ino][jno][2][2] = tmp;
    }
  }
  // compute 1st derivative of energy
  for(int ino=0;ino<4;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = 0;
      for(int jno=0;jno<4;jno++){
        for(int jdim=0;jdim<3;jdim++){
          dW[ino][idim] += ddW[ino][jno][idim][jdim]*c[jno][jdim];
        }
      }
      W += dW[ino][idim]*c[ino][idim];
    }
  }
}

void MakePositiveDefinite_Sim22(const double s2[3],double s3[3])
{
  const double b = (s2[0]+s2[1])*0.5;
  const double d = (s2[0]-s2[1])*(s2[0]-s2[1])*0.25 + s2[2]*s2[2];
  const double e = sqrt(d);
  if( b-e > 1.0e-20 ){
    s3[0] = s2[0];
    s3[1] = s2[1];
    s3[2] = s2[2];
    return;
  }
  if( b+e < 0 ){
    s3[0] = 0;
    s3[1] = 0;
    s3[2] = 0;
    return;
  }
  const double l = b+e;
  double t0[2] = { s2[0]-l, s2[2]   };
  double t1[2] = { s2[2],   s2[1]-l };
  //  std::cout << t0[0]*t1[1]-t0[1]*t1[0] << std::endl;
  const double sqlen_t0 = t0[0]*t0[0]+t0[1]*t0[1];
  const double sqlen_t1 = t1[0]*t1[0]+t1[1]*t1[1];
  if( sqlen_t0 > sqlen_t1 ){
    if( sqlen_t0 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t0 = 1.0/sqrt(sqlen_t0);
    t0[0] *= invlen_t0;
    t0[1] *= invlen_t0;
    s3[0] = l*t0[0]*t0[0];
    s3[1] = l*t0[1]*t0[1];
    s3[2] = l*t0[0]*t0[1];
  }
  else{
    if( sqlen_t1 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t1 = 1.0/sqrt(sqlen_t1);
    t1[0] *= invlen_t1;
    t1[1] *= invlen_t1;
    s3[0] = l*t1[0]*t1[0];
    s3[1] = l*t1[1]*t1[1];
    s3[2] = l*t1[0]*t1[1];
  }
  return;
}


void WdWddW_CST
(double& W, // (out) energy
 double dW[3][3], // (out) 1st derivative of energy
 double ddW[3][3][3][3], // (out) 2nd derivative of energy
 ////
 const double C[3][3], // (in) undeformed triangle vertex positions
 const double c[3][3], // (in) deformed triangle vertex positions
 const double lambda, // (in) Lame's 1st parameter
 const double myu)     // (in) Lame's 2nd parameter
{
  double Gd[3][3] = { // undeformed edge vector
    { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
    { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
  
  double Gu[2][3]; // inverse of Gd
  {
    Cross3D(Gu[0], Gd[1], Gd[2]);
    const double invtmp1 = 1.0/Dot3D(Gu[0],Gd[0]);
    Gu[0][0] *= invtmp1;	Gu[0][1] *= invtmp1;	Gu[0][2] *= invtmp1;
    ////
    Cross3D(Gu[1], Gd[2], Gd[0]);
    const double invtmp2 = 1.0/Dot3D(Gu[1],Gd[1]);
    Gu[1][0] *= invtmp2;	Gu[1][1] *= invtmp2;	Gu[1][2] *= invtmp2;
  }
  
  const double gd[2][3] = { // deformed edge vector
    { c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
    { c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5*( Dot3D(gd[0],gd[0]) - Dot3D(Gd[0],Gd[0]) ),
    0.5*( Dot3D(gd[1],gd[1]) - Dot3D(Gd[1],Gd[1]) ),
    1.0*( Dot3D(gd[0],gd[1]) - Dot3D(Gd[0],Gd[1]) ) };
  const double GuGu2[3] = { Dot3D(Gu[0],Gu[0]), Dot3D(Gu[1],Gu[1]), Dot3D(Gu[1],Gu[0]) };
  const double Cons2[3][3] = { // constitutive tensor
    { lambda*GuGu2[0]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[0]),
      lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[0]*GuGu2[2] + 2*myu*(GuGu2[0]*GuGu2[2]) },
    { lambda*GuGu2[1]*GuGu2[0] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[1]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[1]),
      lambda*GuGu2[1]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[1]) },
    { lambda*GuGu2[2]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[2]),
      lambda*GuGu2[2]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[1]),
      lambda*GuGu2[2]*GuGu2[2] + 1*myu*(GuGu2[0]*GuGu2[1] + GuGu2[2]*GuGu2[2]) } };
  const double S2[3] = {  // 2nd Piola-Kirchhoff stress
    Cons2[0][0]*E2[0] + Cons2[0][1]*E2[1] + Cons2[0][2]*E2[2],
    Cons2[1][0]*E2[0] + Cons2[1][1]*E2[1] + Cons2[1][2]*E2[2],
    Cons2[2][0]*E2[0] + Cons2[2][1]*E2[1] + Cons2[2][2]*E2[2] };
  
  // compute energy
  W = 0.5*Area*(E2[0]*S2[0] + E2[1]*S2[1] + E2[2]*S2[2]);
  
  // compute 1st derivative
  const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };
  for(int ino=0;ino<3;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = Area*
      (+S2[0]*gd[0][idim]*dNdr[ino][0]
       +S2[2]*gd[0][idim]*dNdr[ino][1]
       +S2[2]*gd[1][idim]*dNdr[ino][0]
       +S2[1]*gd[1][idim]*dNdr[ino][1]);
    }
  }
  
  double S3[3] = { S2[0], S2[1], S2[2] };
  MakePositiveDefinite_Sim22(S2,S3);
  
  // compute second derivative
  for(int ino=0;ino<3;ino++){
    for(int jno=0;jno<3;jno++){
      for(int idim=0;idim<3;idim++){
        for(int jdim=0;jdim<3;jdim++){
          double dtmp0 = 0;
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];
          ddW[ino][jno][idim][jdim] = dtmp0*Area;
        }
      }
      const double dtmp1 = Area*
      (+S3[0]*dNdr[ino][0]*dNdr[jno][0]
       +S3[2]*dNdr[ino][0]*dNdr[jno][1]
       +S3[2]*dNdr[ino][1]*dNdr[jno][0]
       +S3[1]*dNdr[ino][1]*dNdr[jno][1]);
      ddW[ino][jno][0][0] += dtmp1;
      ddW[ino][jno][1][1] += dtmp1;
      ddW[ino][jno][2][2] += dtmp1;
    }
  }
}



// compute energy and its 1st and 2nd derivative for contact against object
void WdWddW_Contact
(double& W,  // (out) energy
 double dW[3], // (out) 1st derivative of energy
 double ddW[3][3], // (out) 2nd derivative of energy
 ////
 const double c[3], // (in) deformed triangle vertex positions
 double stiff_contact,
 double contact_clearance,
 const CInput_Contact& input )
{
  double n[3];
  double pd = input.penetrationNormal(n[0],n[1],n[2], c[0],c[1],c[2]);
  pd += contact_clearance;
  if( pd  < 0 ){
    W = 0;
    dW[0] = 0;  dW[1] = 0;  dW[2] = 0;
    ddW[0][0] = 0;  ddW[0][1] = 0;  ddW[0][2] = 0;
    ddW[1][0] = 0;  ddW[1][1] = 0;  ddW[1][2] = 0;
    ddW[2][0] = 0;  ddW[2][1] = 0;  ddW[2][2] = 0;
    return;
  }
  W = 0.5*stiff_contact*pd*pd;
  
  dW[0] = -stiff_contact*pd*n[0];
  dW[1] = -stiff_contact*pd*n[1];
  dW[2] = -stiff_contact*pd*n[2];
  
  ddW[0][0] = stiff_contact*n[0]*n[0];
  ddW[0][1] = stiff_contact*n[0]*n[1];
  ddW[0][2] = stiff_contact*n[0]*n[2];
  ddW[1][0] = stiff_contact*n[1]*n[0];
  ddW[1][1] = stiff_contact*n[1]*n[1];
  ddW[1][2] = stiff_contact*n[1]*n[2];
  ddW[2][0] = stiff_contact*n[2]*n[0];
  ddW[2][1] = stiff_contact*n[2]*n[1];
  ddW[2][2] = stiff_contact*n[2]*n[2];
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


const static unsigned int NIntTetGauss[4] = {
  1, 4, 5, 16
};
const static double TetGauss[4][16][4] =
{
  {	// order-1    1point
    { 0.25, 0.25, 0.25, 1.0 },
  },
  {	// order-2    4point
    { 0.585410196624968, 0.138196601125015, 0.138196601125015, 0.25 },
    { 0.138196601125015, 0.585410196624968, 0.138196601125015, 0.25 },
    { 0.138196601125015, 0.138196601125015, 0.585410196624968, 0.25 },
    { 0.138196601125015, 0.138196601125015, 0.138196601125015, 0.25 },
  },
  {	// order-3    5point
    { 0.25, 0.25, 0.25, -0.8 },
    { 0.5, 0.1666666666666667, 0.1666666666666667, 0.45 },
    { 0.1666666666666667, 0.5, 0.1666666666666667, 0.45 },
    { 0.1666666666666667, 0.1666666666666667, 0.5, 0.45 },
    { 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.45 },
  },
  {	// order-4    16point
    { 0.7716429020672371, 0.07611903264425430, 0.07611903264425430, 0.05037379410012282 },
    { 0.07611903264425430, 0.7716429020672371, 0.07611903264425430, 0.05037379410012282 },
    { 0.07611903264425430, 0.07611903264425430, 0.7716429020672371, 0.05037379410012282 },
    { 0.07611903264425430, 0.07611903264425430, 0.07611903264425430, 0.05037379410012282 },

    { 0.1197005277978019, 0.4042339134672644, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.1197005277978019, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.4042339134672644, 0.1197005277978019, 0.06654206863329239 },

    { 0.07183164526766925, 0.4042339134672644, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.07183164526766925, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.4042339134672644, 0.07183164526766925, 0.06654206863329239 },

    { 0.1197005277978019, 0.07183164526766925, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.1197005277978019, 0.07183164526766925, 0.06654206863329239 },
    { 0.07183164526766925, 0.4042339134672644, 0.1197005277978019, 0.06654206863329239 },

    { 0.07183164526766925, 0.1197005277978019, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.07183164526766925, 0.1197005277978019, 0.06654206863329239 },
    { 0.1197005277978019, 0.4042339134672644, 0.07183164526766925, 0.06654206863329239 },
  }
};


static double TetVolume3D
(const double v1[3],
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
static inline void TetDlDx(double dldx[][3], double a[],
  const double p0[], const double p1[], const double p2[], const double p3[])
{
  const double vol = TetVolume3D(p0, p1, p2, p3);
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

  //	std::cout << dldx[0][0]+dldx[1][0]+dldx[2][0]+dldx[3][0] << std::endl;
  //	std::cout << dldx[0][1]+dldx[1][1]+dldx[2][1]+dldx[3][1] << std::endl;
  //	std::cout << dldx[0][2]+dldx[1][2]+dldx[2][2]+dldx[3][2] << std::endl;

  //	std::cout << a[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1]+dldx[0][2]*p0[2] << std::endl;
  //	std::cout << a[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1]+dldx[1][2]*p1[2] << std::endl;
  //	std::cout << a[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1]+dldx[2][2]*p2[2] << std::endl;
  //	std::cout << a[3]+dldx[3][0]*p3[0]+dldx[3][1]*p3[1]+dldx[3][2]*p3[2] << std::endl;
}

void MakeMat_Poisson3D_P1
(const double alpha, const double source,
 const double coords[4][3],
 const double value[4],
 double eres[4],
 double emat[][4])
{
  const int nno = 4;
  const int ndim = 3;
  ////
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;  eres[3] = 0;
  for (int i = 0; i<16; ++i){ (&emat[0][0])[i] = 0.0; }
  const double area = TetVolume3D(coords[0], coords[1], coords[2], coords[3]);
  ////
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term, coords[0], coords[1], coords[2], coords[3]);
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      emat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2]);
    }
  }
  for (int ino = 0; ino<nno; ino++){
    eres[ino] = source*area*0.25;
  }
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      eres[ino] -= emat[ino][jno]*value[jno];
    }
  }
}

void MakeMat_Diffusion3D_P1
(const double alpha, const double source,
 const double dt_timestep, const double gamma_newmark, const double rho,
 const double coords[4][3],
 const double value[4], const double velo[4],
 double eres[4],
 double emat[4][4])
{
  const int nno = 4;
  const int ndim = 3;
  
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;  eres[3] = 0;
  for (int i=0; i<16; ++i){ (&emat[0][0])[i] = 0.0; }
  
  const double vol = TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx,const_term,coords[0],coords[1],coords[2],coords[3]);
  
  ////////////////////////////////////////////////////////////////
  
  double eCmat[nno][nno];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat[ino][jno] = alpha*vol*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2]);
    }
  }
  double eMmat[nno][nno];
  {
    const double dtmp1 = rho*vol*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno] = dtmp1;
      }
      eMmat[ino][ino] += dtmp1;
    }
  }
		
  for(int ino=0;ino<nno;ino++){
    eres[ino] = source*vol*0.25;
  }
  
  ////////////////////////////////////////////////////////////////
  
  {
    const double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno;i++){
      (&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino]	-= eCmat[ino][jno]*(value[jno]+dt_timestep*velo[jno]) + eMmat[ino][jno]*velo[jno];
    }
  }
}


void stress_LinearSolid_TetP2
(double stress[3][3],
const double l0, const double l1, const double l2, const double l3,
const double vol, const double lambda, const double myu,
const double g_x, const double g_y, const double g_z,
const double dldx[4][3],
const double disp[10][3])
{
  /*
  double N[10] = {
    l0*(2*l0-1),
    l1*(2*l1-1),
    l2*(2*l2-1),
    l3*(2*l3-1),
    4*l0*l1,
    4*l1*l2,
    4*l0*l2,
    4*l0*l3,
    4*l1*l3,
    4*l2*l3,
  };
   */
  double dNdx[10][3];
  for (unsigned int i = 0; i<3; i++){
    dNdx[0][i] = (4*l0-1)*dldx[0][i];
    dNdx[1][i] = (4*l1-1)*dldx[1][i];
    dNdx[2][i] = (4*l2-1)*dldx[2][i];
    dNdx[3][i] = (4*l3-1)*dldx[3][i];
    dNdx[4][i] = 4*dldx[0][i]*l1+4*l0*dldx[1][i];
    dNdx[5][i] = 4*dldx[1][i]*l2+4*l1*dldx[2][i];
    dNdx[6][i] = 4*dldx[0][i]*l2+4*l0*dldx[2][i];
    dNdx[7][i] = 4*dldx[0][i]*l3+4*l0*dldx[3][i];
    dNdx[8][i] = 4*dldx[1][i]*l3+4*l1*dldx[3][i];
    dNdx[9][i] = 4*dldx[2][i]*l3+4*l2*dldx[3][i];
  }
  double dudx[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
  for (unsigned int ino = 0; ino<10; ino++){
    for (unsigned int i = 0; i<3; i++){
      for (unsigned int j = 0; j<3; j++){
        dudx[i][j] += disp[ino][i]*dNdx[ino][j];
      }
    }
  }
  double strain[3][3];
  for (unsigned int i = 0; i<3; i++){
    for (unsigned int j = 0; j<3; j++){
      strain[i][j] = 0.5*(dudx[i][j]+dudx[j][i]);
    }
  }
  {
    for (unsigned int i = 0; i<3; i++){
      for (unsigned int j = 0; j<3; j++){
        stress[i][j] = 2*myu*strain[i][j];
      }
      stress[i][i] += lambda*(strain[0][0]+strain[1][1]+strain[2][2]);
    }
  }
}

void matRes_LinearSolid_TetP2
(double emat[10][10][3][3],
double eres[10][3],
const double vol, const double lambda, const double myu,
const double g_x, const double g_y, const double g_z,
const double rho,
const double dldx[4][3],
const double disp[10][3])
{
  for (unsigned int i = 0; i<10*10*3*3; i++){ (&emat[0][0][0][0])[i] = 0; }
  for (unsigned int i = 0; i<10*3; i++){ (&eres[0][0])[i] = 0; }
  unsigned int nOrder = 2;
  unsigned int nInt = NIntTetGauss[nOrder];
  for (unsigned int iint = 0; iint<nInt; iint++){
    double l0 = TetGauss[nOrder][iint][0];
    double l1 = TetGauss[nOrder][iint][1];
    double l2 = TetGauss[nOrder][iint][2];
    double l3 = (1-l0-l1-l2);
    double w = TetGauss[nOrder][iint][3];
    double N[10] = {
      l0*(2*l0-1),
      l1*(2*l1-1),
      l2*(2*l2-1),
      l3*(2*l3-1),
      4*l0*l1,
      4*l1*l2,
      4*l0*l2,
      4*l0*l3,
      4*l1*l3,
      4*l2*l3,
    };
    double dNdx[10][3];
    for (unsigned int i = 0; i<3; i++){
      dNdx[0][i] = (4*l0-1)*dldx[0][i];
      dNdx[1][i] = (4*l1-1)*dldx[1][i];
      dNdx[2][i] = (4*l2-1)*dldx[2][i];
      dNdx[3][i] = (4*l3-1)*dldx[3][i];
      dNdx[4][i] = 4*dldx[0][i]*l1+4*l0*dldx[1][i];
      dNdx[5][i] = 4*dldx[1][i]*l2+4*l1*dldx[2][i];
      dNdx[6][i] = 4*dldx[0][i]*l2+4*l0*dldx[2][i];
      dNdx[7][i] = 4*dldx[0][i]*l3+4*l0*dldx[3][i];
      dNdx[8][i] = 4*dldx[1][i]*l3+4*l1*dldx[3][i];
      dNdx[9][i] = 4*dldx[2][i]*l3+4*l2*dldx[3][i];
    }
    /*
    {
    double tN[4] = {0,0,0,0};
    for(unsigned int ino=0;ino<10;ino++){
    tN[0] += N[ino];
    tN[1] += dNdx[ino][0];
    tN[2] += dNdx[ino][1];
    tN[3] += dNdx[ino][2];
    }
    std::cout << tN[0] << "   " << tN[1] << " " << tN[2] << " " << tN[3] << std::endl;
    }
    */
    for (unsigned int ino = 0; ino<10; ino++){
      for (unsigned int jno = 0; jno<10; jno++){
        emat[ino][jno][0][0] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][0]+myu*dNdx[jno][0]*dNdx[ino][0]);
        emat[ino][jno][0][1] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][1]+myu*dNdx[jno][0]*dNdx[ino][1]);
        emat[ino][jno][0][2] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][2]+myu*dNdx[jno][0]*dNdx[ino][2]);
        emat[ino][jno][1][0] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][0]+myu*dNdx[jno][1]*dNdx[ino][0]);
        emat[ino][jno][1][1] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][1]+myu*dNdx[jno][1]*dNdx[ino][1]);
        emat[ino][jno][1][2] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][2]+myu*dNdx[jno][1]*dNdx[ino][2]);
        emat[ino][jno][2][0] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][0]+myu*dNdx[jno][2]*dNdx[ino][0]);
        emat[ino][jno][2][1] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][1]+myu*dNdx[jno][2]*dNdx[ino][1]);
        emat[ino][jno][2][2] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][2]+myu*dNdx[jno][2]*dNdx[ino][2]);
        const double dtmp1 = dNdx[ino][0]*dNdx[jno][0]+dNdx[ino][1]*dNdx[jno][1]+dNdx[ino][2]*dNdx[jno][2];
        emat[ino][jno][0][0] += w*vol*myu*dtmp1;
        emat[ino][jno][1][1] += w*vol*myu*dtmp1;
        emat[ino][jno][2][2] += w*vol*myu*dtmp1;
      }
    }
    for (unsigned int ino = 0; ino<10; ino++){
      eres[ino][0] += w*vol*rho*g_x*N[ino];
      eres[ino][1] += w*vol*rho*g_y*N[ino];
      eres[ino][2] += w*vol*rho*g_z*N[ino];
      for (unsigned int jno = 0; jno<10; jno++){
        eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
        eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
        eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
      }
    }
  }
}


void MakeMat_LinearSolid3D_Static_P1
(const double myu, const double lambda,
const double rho, const double g_x, const double g_y, const double g_z, 
const double coords[4][3],
const double disp[4][3],
////
double emat[4][4][3][3],
double eres[4][3])
{
  const double vol = TetVolume3D(coords[0], coords[1], coords[2], coords[3]);
  double dldx[4][3];
  {
    double const_term[4];    
    TetDlDx(dldx, const_term, coords[0], coords[1], coords[2], coords[3]);
  }
  /////////////////////////////////////////
  for (int ino = 0; ino<4; ino++){
    for (int jno = 0; jno<4; jno++){
      emat[ino][jno][0][0] = vol*(lambda*dldx[ino][0]*dldx[jno][0]+myu*dldx[jno][0]*dldx[ino][0]);
      emat[ino][jno][0][1] = vol*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
      emat[ino][jno][0][2] = vol*(lambda*dldx[ino][0]*dldx[jno][2]+myu*dldx[jno][0]*dldx[ino][2]);
      emat[ino][jno][1][0] = vol*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
      emat[ino][jno][1][1] = vol*(lambda*dldx[ino][1]*dldx[jno][1]+myu*dldx[jno][1]*dldx[ino][1]);
      emat[ino][jno][1][2] = vol*(lambda*dldx[ino][1]*dldx[jno][2]+myu*dldx[jno][1]*dldx[ino][2]);
      emat[ino][jno][2][0] = vol*(lambda*dldx[ino][2]*dldx[jno][0]+myu*dldx[jno][2]*dldx[ino][0]);
      emat[ino][jno][2][1] = vol*(lambda*dldx[ino][2]*dldx[jno][1]+myu*dldx[jno][2]*dldx[ino][1]);
      emat[ino][jno][2][2] = vol*(lambda*dldx[ino][2]*dldx[jno][2]+myu*dldx[jno][2]*dldx[ino][2]);
      const double dtmp1 = dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2];
      emat[ino][jno][0][0] += vol*myu*dtmp1;
      emat[ino][jno][1][1] += vol*myu*dtmp1;
      emat[ino][jno][2][2] += vol*myu*dtmp1;
    }
  }
  for (int ino = 0; ino<4; ino++){
    eres[ino][0] = vol*rho*g_x*0.25;
    eres[ino][1] = vol*rho*g_y*0.25;
    eres[ino][2] = vol*rho*g_z*0.25;
    for (int jno = 0; jno<4; jno++){
      eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
      eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
      eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
    }
  }
}

void MakeMat_LinearSolid3D_Static_Q1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y, const double g_z,
 const double coords[8][3],
 const double disp[8][3],
 ////
 double emat[8][8][3][3],
 double eres[8][3])
{
  const int nDegInt = 2;
  const int nInt = NIntLineGauss[nDegInt];
  const double (*Gauss)[2] = LineGauss[nDegInt];
  
  for(unsigned int i=0;i<8*8*3*3;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
  for(unsigned int i=0;i<    8*3;i++){ *( &eres[0][0]      +i) = 0.0; }
  
  double vol = 0.0;
  for(int ir1=0;ir1<nInt;ir1++){
  for(int ir2=0;ir2<nInt;ir2++){
  for(int ir3=0;ir3<nInt;ir3++){
    const double r1 = Gauss[ir1][0];
    const double r2 = Gauss[ir2][0];
    const double r3 = Gauss[ir3][0];
    double detjac, detwei, dndx[8][3], an[8];
    ShapeFunc_Hex8(r1,r2,r3,coords,detjac,dndx,an);
    detwei = detjac*Gauss[ir1][1]*Gauss[ir2][1]*Gauss[ir3][1];
    vol += detwei;
    for(int ino=0;ino<8;ino++){
    for(int jno=0;jno<8;jno++){
      double dtmp1 = 0.0;
      for(int idim=0;idim<3;idim++){
        for(int jdim=0;jdim<3;jdim++){
//          emat[ino][jno][idim][jdim] += detwei*( lambda*dndx[ino][idim]*dndx[jno][jdim]
//                                                +myu*dndx[jno][idim]*dndx[ino][jdim] );
        }
        dtmp1 += dndx[ino][idim]*dndx[jno][idim];
      }
      for(int idim=0;idim<3;idim++){
        emat[ino][jno][idim][idim] += detwei*myu*dtmp1;
      }
    }
    }
    for(int ino=0;ino<8;ino++){
      eres[ino][0] += detwei*rho*g_x*an[ino];
      eres[ino][1] += detwei*rho*g_y*an[ino];
      eres[ino][2] += detwei*rho*g_z*an[ino];
    }
  }
  }
  }
  for (int ino = 0; ino<8; ino++){
  for (int jno = 0; jno<8; jno++){
    eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
    eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
    eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
  }
  }
}

void MakeMat_LinearSolid3D_Dynamic_P1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y, const double g_z,
 const double dt, const double gamma_newmark,  const double beta_newmark,
 const double disp[4][3], const double velo[4][3], const double acc[4][3],
 const double coords[4][3],
 double eres[4][3],
 double emat[4][4][3][3],
 bool is_initial)
{
  const int nno = 4;
  const int ndim = 3;
  
  const double vol = TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim];		// spatial derivative of linear shape function
  {
    double zero_order_term[nno];	// const term of shape function
    TetDlDx(dldx, zero_order_term,   coords[0],coords[1],coords[2],coords[3]);
  }
  
  double eKmat[nno][nno][ndim][ndim];
  for (int ino = 0; ino<4; ino++){
    for (int jno = 0; jno<4; jno++){
      eKmat[ino][jno][0][0] = vol*(lambda*dldx[ino][0]*dldx[jno][0]+myu*dldx[jno][0]*dldx[ino][0]);
      eKmat[ino][jno][0][1] = vol*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
      eKmat[ino][jno][0][2] = vol*(lambda*dldx[ino][0]*dldx[jno][2]+myu*dldx[jno][0]*dldx[ino][2]);
      eKmat[ino][jno][1][0] = vol*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
      eKmat[ino][jno][1][1] = vol*(lambda*dldx[ino][1]*dldx[jno][1]+myu*dldx[jno][1]*dldx[ino][1]);
      eKmat[ino][jno][1][2] = vol*(lambda*dldx[ino][1]*dldx[jno][2]+myu*dldx[jno][1]*dldx[ino][2]);
      eKmat[ino][jno][2][0] = vol*(lambda*dldx[ino][2]*dldx[jno][0]+myu*dldx[jno][2]*dldx[ino][0]);
      eKmat[ino][jno][2][1] = vol*(lambda*dldx[ino][2]*dldx[jno][1]+myu*dldx[jno][2]*dldx[ino][1]);
      eKmat[ino][jno][2][2] = vol*(lambda*dldx[ino][2]*dldx[jno][2]+myu*dldx[jno][2]*dldx[ino][2]);
      const double dtmp1 = dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2];
      eKmat[ino][jno][0][0] += vol*myu*dtmp1;
      eKmat[ino][jno][1][1] += vol*myu*dtmp1;
      eKmat[ino][jno][2][2] += vol*myu*dtmp1;
    }
  }
  
  double eMmat[nno][nno][ndim][ndim];
  {
    const double dtmp1 = vol*rho*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno][0][0] = dtmp1;
        eMmat[ino][jno][0][1] = 0.0;
        eMmat[ino][jno][0][2] = 0.0;
        eMmat[ino][jno][1][0] = 0.0;
        eMmat[ino][jno][1][1] = dtmp1;
        eMmat[ino][jno][1][2] = 0.0;
        eMmat[ino][jno][2][0] = 0.0;
        eMmat[ino][jno][2][1] = 0.0;
        eMmat[ino][jno][2][2] = dtmp1;
      }
      eMmat[ino][ino][0][0] += dtmp1;
      eMmat[ino][ino][1][1] += dtmp1;
      eMmat[ino][ino][2][2] += dtmp1;
    }
  }
  
  // calc external force
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = vol*rho*g_x*0.25;
    eres[ino][1] = vol*rho*g_y*0.25;
    eres[ino][2] = vol*rho*g_z*0.25;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  {	// calc coeff matrix for newmark-beta
    double dtmp1 = beta_newmark*dt*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
    }
  }
  
  // calc element redisual vector
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eKmat[ino][jno][0][0]*disp[jno][0]+eKmat[ino][jno][0][1]*disp[jno][1]+eKmat[ino][jno][0][2]*disp[jno][2];
      eres[ino][1] -= eKmat[ino][jno][1][0]*disp[jno][0]+eKmat[ino][jno][1][1]*disp[jno][1]+eKmat[ino][jno][1][2]*disp[jno][2];
      eres[ino][2] -= eKmat[ino][jno][2][0]*disp[jno][0]+eKmat[ino][jno][2][1]*disp[jno][1]+eKmat[ino][jno][2][2]*disp[jno][2];
    }
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eMmat[ino][jno][0][0]*acc[jno][0]+eMmat[ino][jno][0][1]*acc[jno][1]+eMmat[ino][jno][0][2]*acc[jno][2];
      eres[ino][1] -= eMmat[ino][jno][1][0]*acc[jno][0]+eMmat[ino][jno][1][1]*acc[jno][1]+eMmat[ino][jno][1][2]*acc[jno][2];
      eres[ino][2] -= eMmat[ino][jno][2][0]*acc[jno][0]+eMmat[ino][jno][2][1]*acc[jno][1]+eMmat[ino][jno][2][2]*acc[jno][2];
    }
  }
  if( is_initial ){
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= dt*(eKmat[ino][jno][0][0]*velo[jno][0]+eKmat[ino][jno][0][1]*velo[jno][1]+eKmat[ino][jno][0][2]*velo[jno][2]);
        eres[ino][1] -= dt*(eKmat[ino][jno][1][0]*velo[jno][0]+eKmat[ino][jno][1][1]*velo[jno][1]+eKmat[ino][jno][1][2]*velo[jno][2]);
        eres[ino][2] -= dt*(eKmat[ino][jno][2][0]*velo[jno][0]+eKmat[ino][jno][2][1]*velo[jno][1]+eKmat[ino][jno][2][2]*velo[jno][2]);
      }
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= 0.5*dt*dt*(eKmat[ino][jno][0][0]*acc[jno][0]+eKmat[ino][jno][0][1]*acc[jno][1]+eKmat[ino][jno][0][2]*acc[jno][2]);
        eres[ino][1] -= 0.5*dt*dt*(eKmat[ino][jno][1][0]*acc[jno][0]+eKmat[ino][jno][1][1]*acc[jno][1]+eKmat[ino][jno][1][2]*acc[jno][2]);
        eres[ino][2] -= 0.5*dt*dt*(eKmat[ino][jno][2][0]*acc[jno][0]+eKmat[ino][jno][2][1]*acc[jno][1]+eKmat[ino][jno][2][2]*acc[jno][2]);
      }
    }
  }
}

void MakeMat_Stokes3D_Static_P1P1
(double alpha, double g_x, double g_y, double g_z,
 const double coords[4][3],
 const double velo[4][3], const double press[4],
 double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
 double eres_u[4][3], double eres_p[4])
{
  const unsigned int nno = 4;
  const unsigned int ndim = 3;
  
  const double vol = TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);
  
  ////////////////////////////////////////////////////////////
  
  for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat_uu[0][0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      const double dtmp1 = vol*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      emat_uu[ino][jno][0][0] = dtmp1;
      emat_uu[ino][jno][1][1] = dtmp1;
      emat_uu[ino][jno][2][2] = dtmp1;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_up[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_up[ino][jno][0] += vol*dldx[ino][0]*0.25;
      emat_up[ino][jno][1] += vol*dldx[ino][1]*0.25;
      emat_up[ino][jno][2] += vol*dldx[ino][2]*0.25;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_pu[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pu[ino][jno][0] += vol*dldx[jno][0]*0.25;
      emat_pu[ino][jno][1] += vol*dldx[jno][1]*0.25;
      emat_pu[ino][jno][2] += vol*dldx[jno][2]*0.25;
    }
  }
  for(unsigned int i=0;i<nno*nno;i++){ *(&emat_pp[0][0]+i) = 0.0; }
  double tau; // relaxation parameter
  {
    const double h = pow(vol/3.14, 0.3333333333)*2;
    tau = -h*h/alpha*0.1;
  }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pp[ino][jno] = vol*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  for(unsigned int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*g_x*0.25;
    eres_u[ino][1] = vol*g_y*0.25;
    eres_u[ino][2] = vol*g_z*0.25;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  eres_p[3] = 0;
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_uu[ino][jno][0][0]*velo[jno][0]+emat_uu[ino][jno][0][1]*velo[jno][1]+emat_uu[ino][jno][0][2]*velo[jno][2];
      eres_u[ino][1] -= emat_uu[ino][jno][1][0]*velo[jno][0]+emat_uu[ino][jno][1][1]*velo[jno][1]+emat_uu[ino][jno][1][2]*velo[jno][2];
      eres_u[ino][2] -= emat_uu[ino][jno][2][0]*velo[jno][0]+emat_uu[ino][jno][2][1]*velo[jno][1]+emat_uu[ino][jno][2][2]*velo[jno][2];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_up[ino][jno][0]*press[jno];
      eres_u[ino][1] -= emat_up[ino][jno][1]*press[jno];
      eres_u[ino][2] -= emat_up[ino][jno][2]*press[jno];
    }
  }
  for(unsigned int ino=0;ino<nno;ino++){
    eres_p[ino] = 0.0;
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pu[ino][jno][0]*velo[jno][0]+emat_pu[ino][jno][1]*velo[jno][1]+emat_pu[ino][jno][2]*velo[jno][2];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pp[ino][jno]*press[jno];
    }
  }
}



void MakeMat_Stokes3D_Static_P1
(double alpha, double g_x, double g_y, double g_z,
 const double coords[4][3],
 const double velo_press[4][4],
 double emat[4][4][4][4],
 double eres[4][4])
{
  const int nno = 4;
  
  const double velo[4][3] = {
    {velo_press[0][0],velo_press[0][1],velo_press[0][2]},
    {velo_press[1][0],velo_press[1][1],velo_press[1][2]},
    {velo_press[2][0],velo_press[2][1],velo_press[2][2]},
    {velo_press[3][0],velo_press[3][1],velo_press[3][2]} };
  const double press[4] = { velo_press[0][3], velo_press[1][3], velo_press[2][3], velo_press[3][3] };
  ////
  double emat_uu[4][4][3][3], emat_up[4][4][3], emat_pu[4][4][3], emat_pp[4][4];
  double eres_u[4][3], eres_p[4];
  MakeMat_Stokes3D_Static_P1P1(alpha, g_x, g_y, g_z,
                               coords, velo, press,
                               emat_uu, emat_up, emat_pu, emat_pp,
                               eres_u, eres_p);
  ////
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
      emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
      emat[ino][jno][0][2] = emat_uu[ino][jno][0][2];
      ////
      emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
      emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
      emat[ino][jno][1][2] = emat_uu[ino][jno][1][2];
      ////
      emat[ino][jno][2][0] = emat_uu[ino][jno][2][0];
      emat[ino][jno][2][1] = emat_uu[ino][jno][2][1];
      emat[ino][jno][2][2] = emat_uu[ino][jno][2][2];
      ////
      emat[ino][jno][0][3] = emat_up[ino][jno][0];
      emat[ino][jno][1][3] = emat_up[ino][jno][1];
      emat[ino][jno][2][3] = emat_up[ino][jno][2];
      ////
      emat[ino][jno][3][0] = emat_pu[ino][jno][0];
      emat[ino][jno][3][1] = emat_pu[ino][jno][1];
      emat[ino][jno][3][2] = emat_pu[ino][jno][2];
      ////
      emat[ino][jno][3][3] = emat_pp[ino][jno];
    }
  }
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = eres_u[ino][0];
    eres[ino][1] = eres_u[ino][1];
    eres[ino][2] = eres_u[ino][2];
    eres[ino][3] = eres_p[ino];
  }
}



void MakeMat_Stokes3D_Dynamic_Newmark_P1P1
(double alpha, double rho, double g_x, double g_y, double g_z,
 const double dt_timestep, const double gamma_newmark,
 const double coords[4][3],
 const double velo[4][3], const double press[4], const double acc[4][3], const double apress[4],
 double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
 double eres_u[4][3], double eres_p[4])
{
  //	std::cout << "AddMat_Stokes2D_NonStatic_Newmark_P1P1" << std::endl;
  
  const int nno = 4;
  const int ndim = 3;
  ////
  const double vol = TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0],coords[1],coords[2],coords[3]);
  
  //////////////////////////////////////////////////////////////////////////////
  double eCmat_uu[4][4][3][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = vol*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      eCmat_uu[ino][jno][0][0] = dtmp1;
      eCmat_uu[ino][jno][0][1] = 0.0;
      eCmat_uu[ino][jno][0][2] = 0.0;
      ////
      eCmat_uu[ino][jno][1][0] = 0.0;
      eCmat_uu[ino][jno][1][1] = dtmp1;
      eCmat_uu[ino][jno][1][2] = 0.0;
      ////
      eCmat_uu[ino][jno][2][0] = 0.0;
      eCmat_uu[ino][jno][2][1] = 0.0;
      eCmat_uu[ino][jno][2][2] = dtmp1;
    }
  }
  
  double eCmat_up[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] = vol*dldx[ino][0]*0.25;
      eCmat_up[ino][jno][1] = vol*dldx[ino][1]*0.25;
      eCmat_up[ino][jno][2] = vol*dldx[ino][2]*0.25;
    }
  }
  
  double eCmat_pu[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pu[ino][jno][0] = vol*dldx[jno][0]*0.25;
      eCmat_pu[ino][jno][1] = vol*dldx[jno][1]*0.25;
      eCmat_pu[ino][jno][2] = vol*dldx[jno][2]*0.25;
    }
  }
  
  double tau;
  {
    const double h = pow( vol / 3.14, 0.3333333 )*2;
    tau = -h*h/alpha*0.1;
  }
  
  double eCmat_pp[4][4];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pp[ino][jno] = vol*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  ////////////////
  double eMmat_uu[4][4][3][3];
  {
    const double dtmp1 = vol*rho*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat_uu[ino][jno][0][0] = dtmp1;
        eMmat_uu[ino][jno][0][1] = 0.0;
        eMmat_uu[ino][jno][0][2] = 0.0;
        ////
        eMmat_uu[ino][jno][1][0] = 0.0;
        eMmat_uu[ino][jno][1][1] = dtmp1;
        eMmat_uu[ino][jno][1][2] = 0.0;
        ////
        eMmat_uu[ino][jno][2][0] = 0.0;
        eMmat_uu[ino][jno][2][1] = 0.0;
        eMmat_uu[ino][jno][2][2] = dtmp1;
      }
      eMmat_uu[ino][ino][0][0] += dtmp1;
      eMmat_uu[ino][ino][1][1] += dtmp1;
      eMmat_uu[ino][ino][2][2] += dtmp1;
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*g_x*0.25;
    eres_u[ino][1] = vol*g_y*0.25;
    eres_u[ino][2] = vol*g_z*0.25;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  eres_p[3] = 0;
  
  ////////////////////////////////
  
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          for(int jdim=0;jdim<ndim;jdim++){
            emat_uu[ino][jno][idim][jdim] = eMmat_uu[ino][jno][idim][jdim]+dtmp1*eCmat_uu[ino][jno][idim][jdim];
          }
        }
      }
    }
    
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          emat_up[ino][jno][idim] = dtmp1*eCmat_up[ino][jno][idim];
          emat_pu[ino][jno][idim] = dtmp1*eCmat_pu[ino][jno][idim];
        }
      }
    }
    
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        emat_pp[ino][jno] = dtmp1*eCmat_pp[ino][jno];
      }
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<ndim;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<ndim;jdim++){
          eres_u[ino][idim] -= eCmat_uu[ino][jno][idim][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim])
          + eMmat_uu[ino][jno][idim][jdim]*acc[jno][jdim];
        }
      }
      for(int jno=0;jno<nno;jno++){
        eres_u[ino][idim] -= eCmat_up[ino][jno][idim]*(press[jno]+dt_timestep*apress[jno]);
      }
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      for(int jdim=0;jdim<ndim;jdim++){
        eres_p[ino] -= eCmat_pu[ino][jno][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim]);
      }
    }
    for(int jno=0;jno<nno;jno++){
      eres_p[ino] -= eCmat_pp[ino][jno]*(press[jno]+dt_timestep*apress[jno]);
    }
  }
}

void MakeMat_Stokes3D_Dynamic_P1
(double alpha, double rho, double g_x, double g_y, double g_z,
 const double dt_timestep, const double gamma_newmark,
 const double coords[4][3],
 const double velo_press[4][4], const double acc_apress[4][4],
 double emat[4][4][4][4],
 double eres[4][4])
{
  const int nno = 4;
  const int ndim = 4;
  ////
  const double vol = TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][3], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);
  ////
  double tau;
  {
    const double h = pow( vol / 3.14, 0.333333333 )*2;
    tau = -h*h/alpha*0.1;
  }
  //////////////////////////////////////////////////////////////////////////////
  double eCmat[4][4][4][4];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = vol*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      eCmat[ino][jno][0][0] = dtmp1;
      eCmat[ino][jno][1][1] = dtmp1;
      eCmat[ino][jno][2][2] = dtmp1;
      eCmat[ino][jno][0][1] = 0.0;
      eCmat[ino][jno][0][2] = 0.0;
      eCmat[ino][jno][1][0] = 0.0;
      eCmat[ino][jno][1][2] = 0.0;
      eCmat[ino][jno][2][0] = 0.0;
      eCmat[ino][jno][2][1] = 0.0;
      eCmat[ino][jno][3][0] = vol*dldx[jno][0]*0.25;
      eCmat[ino][jno][3][1] = vol*dldx[jno][1]*0.25;
      eCmat[ino][jno][3][2] = vol*dldx[jno][2]*0.25;
      eCmat[ino][jno][0][3] = vol*dldx[ino][0]*0.25;
      eCmat[ino][jno][1][3] = vol*dldx[ino][1]*0.25;
      eCmat[ino][jno][2][3] = vol*dldx[ino][2]*0.25;
      eCmat[ino][jno][3][3] = vol*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  ////////////////
  double eMmat[4][4][4][4];
  {
    const double dtmp1 = vol*rho*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno][0][0] = dtmp1;
        eMmat[ino][jno][0][1] = 0.0;
        eMmat[ino][jno][0][2] = 0.0;
        eMmat[ino][jno][0][3] = 0.0;
        ////
        eMmat[ino][jno][1][0] = 0.0;
        eMmat[ino][jno][1][1] = dtmp1;
        eMmat[ino][jno][1][2] = 0.0;
        eMmat[ino][jno][1][3] = 0.0;
        ////
        eMmat[ino][jno][2][0] = 0.0;
        eMmat[ino][jno][2][1] = 0.0;
        eMmat[ino][jno][2][2] = dtmp1;
        eMmat[ino][jno][2][3] = 0.0;
        /////
        eMmat[ino][jno][3][0] = 0.0;
        eMmat[ino][jno][3][1] = 0.0;
        eMmat[ino][jno][3][2] = 0.0;
        eMmat[ino][jno][3][3] = 0.0;
      }
      eMmat[ino][ino][0][0] += dtmp1;
      eMmat[ino][ino][1][1] += dtmp1;
      eMmat[ino][ino][2][2] += dtmp1;
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = vol*g_x*0.25;
    eres[ino][1] = vol*g_y*0.25;
    eres[ino][2] = vol*g_z*0.25;
    eres[ino][3] = 0.0;
  }
  
  ////////////////////////////////
  
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno*ndim*ndim;++i){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eCmat[0][0][0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<ndim;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<ndim;jdim++){
          eres[ino][idim]
          -= eCmat[ino][jno][idim][jdim]*(velo_press[jno][jdim]+dt_timestep*acc_apress[jno][jdim])
          + eMmat[ino][jno][idim][jdim]*acc_apress[jno][jdim];
        }
      }
    }
  }
}


void MakeMat_NavierStokes3D_Dynamic_Newmark_P1P1
(double rho, double myu, double g_x, double g_y, double g_z,
 double dt, double gamma,
 const double coords[4][3],
 const double velo[4][3], const double press[4], const double acc[4][3], const double apress[4],
 double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
 double eres_u[4][3], double eres_p[4])
{
  const int nno = 4;
  const int ndim = 3;
  ////
  const double vol = TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);

  const double velo_ave[3] = {
    (velo[0][0]+velo[1][0]+velo[2][0]+velo[3][0])*0.25,
    (velo[0][1]+velo[1][1]+velo[2][1]+velo[3][1])*0.25,
    (velo[0][2]+velo[1][2]+velo[2][2]+velo[3][2])*0.25 };

  ////
  double tau;
  {  // Calc Stabilization Parameter

    const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]+velo_ave[2]*velo_ave[2]);
    double tmp = 0.0;
    for(int ino=0;ino<nno;ino++){
      double tmp1 = dldx[ino][0]*velo_ave[0]+dldx[ino][1]*velo_ave[1]+dldx[ino][2]*velo_ave[2];
      tmp += fabs(tmp1);
    }
    double h_ugn = 0;
    if( fabs(tmp) > 1.0e-10 ){
      h_ugn = 2.0 * norm_v / tmp;
    }
    double tau_sugn1 = 0.0;
    if( fabs(norm_v) > 1.0e-10 ){
      tau_sugn1 = h_ugn * 0.5 / norm_v;
    }
    double tau_sugn2 = dt * 0.50;
    double tau_sugn3 = h_ugn * h_ugn * 0.25 / myu;
    double inv1 = 0.0;
    if( fabs(tau_sugn1) > 1.0e-10 ){
      inv1 = 1.0 / tau_sugn1;
    }
    double inv2 = 0.0;
    if( fabs(tau_sugn2) > 1.0e-10 ){
      inv2 = 1.0 / tau_sugn2;
    }
    double inv3 = 0.0;
    if( fabs(tau_sugn3) > 1.0e-10 ){
      inv3 = 1.0 / tau_sugn3;
    }
    double tau_supg = sqrt( inv1*inv1 + inv2*inv2 + inv3*inv3 );
    tau_supg = 1.0 / tau_supg;
//    double tau_pspg = tau_supg;
    tau = tau_supg*0.25;
//    std::cout << tau << "   " << tau_sugn1 << " " << tau_sugn2 << " " << tau_sugn3 << " " << inv1 << " " << inv2 << " " << inv3 << std::endl;
    /*
  ////
    const double h = pow(vol/3.14, 0.333333333333)*2;
    const double tau_c = h*0.5/norm_v;
    const double cou_c = norm_v*dt/h;
    if( norm_v*h*rho*1.0e-30 > myu ){ // Re = \infty
      const double dtmp1 = 1/(cou_c*cou_c)+1;
      tau = tau_c / sqrt(dtmp1);
    }
    else if( norm_v*h*rho < myu*1.0e-30 ){ // Re = 0
      tau = h*h*rho*0.5/myu;
    }
    else{
      const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
      const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
      tau = tau_c / sqrt(dtmp1);
    }
     */
  }

  /*
  double ave_velo[3];
  { // average velocity
    ave_velo[0] = 0.0;
    ave_velo[1] = 0.0;
    ave_velo[2] = 0.0;
    for(int ino=0;ino<nno;ino++){
      ave_velo[0] += velo[ino][0];
      ave_velo[1] += velo[ino][1];
      ave_velo[2] += velo[ino][2];
    }
    ave_velo[0] *= 0.25;
    ave_velo[1] *= 0.25;
    ave_velo[2] *= 0.25;
  }
  */
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  double eCmat_uu[4][4][3][3];
  // viscousity
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_uu[ino][jno][0][0] = vol*myu*dldx[ino][0]*dldx[jno][0];
      eCmat_uu[ino][jno][0][1] = vol*myu*dldx[ino][1]*dldx[jno][0];
      eCmat_uu[ino][jno][0][2] = vol*myu*dldx[ino][2]*dldx[jno][0];
      ///
      eCmat_uu[ino][jno][1][0] = vol*myu*dldx[ino][0]*dldx[jno][1];
      eCmat_uu[ino][jno][1][1] = vol*myu*dldx[ino][1]*dldx[jno][1];
      eCmat_uu[ino][jno][1][2] = vol*myu*dldx[ino][2]*dldx[jno][1];
      ///
      eCmat_uu[ino][jno][2][0] = vol*myu*dldx[ino][0]*dldx[jno][2];
      eCmat_uu[ino][jno][2][1] = vol*myu*dldx[ino][1]*dldx[jno][2];
      eCmat_uu[ino][jno][2][2] = vol*myu*dldx[ino][2]*dldx[jno][2];
      ////
      const double dtmp1 = vol*myu*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      eCmat_uu[ino][jno][0][0] += dtmp1;
      eCmat_uu[ino][jno][1][1] += dtmp1;
      eCmat_uu[ino][jno][2][2] += dtmp1;
    }
  }
  {	// advection
    const double dtmp0[3] = {
      velo[0][0]+velo[1][0]+velo[2][0]+velo[3][0],
      velo[0][1]+velo[1][1]+velo[2][1]+velo[3][1],
      velo[0][2]+velo[1][2]+velo[2][2]+velo[3][2] };
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = (dldx[jno][0]*dtmp0[0]+dldx[jno][1]*dtmp0[1]+dldx[jno][2]*dtmp0[2]);
      for(int ino=0;ino<nno;ino++){
        double dtmp2 = dtmp1 + (dldx[jno][0]*velo[ino][0]+dldx[jno][1]*velo[ino][1]+dldx[jno][2]*velo[ino][2]);
        dtmp2 *= vol*rho*0.05;
        eCmat_uu[ino][jno][0][0] += dtmp2;
        eCmat_uu[ino][jno][1][1] += dtmp2;
        eCmat_uu[ino][jno][2][2] += dtmp2;
      }
    }
  }
  {	// SUPG for advection term
    double tmp_mat[ndim][ndim] = { {0,0,0}, {0,0,0}, {0,0,0} };
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        tmp_mat[0][0] += velo[ino][0]*velo[jno][0];
        tmp_mat[0][1] += velo[ino][0]*velo[jno][1];
        tmp_mat[0][2] += velo[ino][0]*velo[jno][2];
        ////
        tmp_mat[1][0] += velo[ino][1]*velo[jno][0];
        tmp_mat[1][1] += velo[ino][1]*velo[jno][1];
        tmp_mat[1][2] += velo[ino][1]*velo[jno][2];
        ////
        tmp_mat[2][0] += velo[ino][2]*velo[jno][0];
        tmp_mat[2][1] += velo[ino][2]*velo[jno][1];
        tmp_mat[2][2] += velo[ino][2]*velo[jno][2];
      }
      tmp_mat[0][0] += velo[ino][0]*velo[ino][0];
      tmp_mat[1][1] += velo[ino][1]*velo[ino][1];
      tmp_mat[2][2] += velo[ino][2]*velo[ino][2];
    }
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        double dtmp1 = 0.0;
        dtmp1 +=
         dldx[ino][0]*dldx[jno][0]*tmp_mat[0][0]
        +dldx[ino][0]*dldx[jno][1]*tmp_mat[0][1]
        +dldx[ino][0]*dldx[jno][2]*tmp_mat[0][2]
        +dldx[ino][1]*dldx[jno][0]*tmp_mat[1][0]
        +dldx[ino][1]*dldx[jno][1]*tmp_mat[1][1]
        +dldx[ino][1]*dldx[jno][2]*tmp_mat[1][2]
        +dldx[ino][2]*dldx[jno][0]*tmp_mat[2][0]
        +dldx[ino][2]*dldx[jno][1]*tmp_mat[2][1]
        +dldx[ino][2]*dldx[jno][2]*tmp_mat[2][2];
        dtmp1 *= tau*rho*vol*0.05;
        eCmat_uu[ino][jno][0][0] += dtmp1;
        eCmat_uu[ino][jno][1][1] += dtmp1;
        eCmat_uu[ino][jno][2][2] += dtmp1;
      }
    }
  }
  
  double eCmat_up[4][4][3], eCmat_pu[4][4][3];
  {   // add pressure gradient, non compression term
    const double dtmp1 = vol*0.25;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eCmat_up[ino][jno][0] = -dtmp1*dldx[ino][0];
        eCmat_up[ino][jno][1] = -dtmp1*dldx[ino][1];
        eCmat_up[ino][jno][2] = -dtmp1*dldx[ino][2];
        eCmat_pu[ino][jno][0] = +dtmp1*dldx[jno][0];
        eCmat_pu[ino][jno][1] = +dtmp1*dldx[jno][1];
        eCmat_pu[ino][jno][2] = +dtmp1*dldx[jno][2];
      }
    }
  }
  // SUPG for pressure gradient
  for(int ino=0;ino<nno;ino++){
    const double dtmp1 = (dldx[ino][0]*velo_ave[0]+dldx[ino][1]*velo_ave[1]+dldx[ino][2]*velo_ave[2])*tau*vol;
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] += dtmp1*dldx[jno][0];
      eCmat_up[ino][jno][1] += dtmp1*dldx[jno][1];
      eCmat_up[ino][jno][2] += dtmp1*dldx[jno][2];
    }
  }
  // PSPG for advection term
  for(int jno=0;jno<nno;jno++){
    const double dtmp1 = (dldx[jno][0]*velo_ave[0]+dldx[jno][1]*velo_ave[1]+dldx[jno][2]*velo_ave[2])*tau*vol;
    for(int ino=0;ino<nno;ino++){
      eCmat_pu[ino][jno][0] += dtmp1*dldx[ino][0];
      eCmat_pu[ino][jno][1] += dtmp1*dldx[ino][1];
      eCmat_pu[ino][jno][2] += dtmp1*dldx[ino][2];
    }
  }
  
  double eCmat_pp[4][4];
  // PSPG for pressure term
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pp[ino][jno] = vol*tau/rho*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  ////////////////
  double eMmat_uu[4][4][3][3];
  {	// add inertia term
    const double dtmp1 = vol*rho*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat_uu[ino][jno][0][0] = dtmp1;
        eMmat_uu[ino][jno][0][1] = 0.0;
        eMmat_uu[ino][jno][0][2] = 0.0;
        ////
        eMmat_uu[ino][jno][1][0] = 0.0;
        eMmat_uu[ino][jno][1][1] = dtmp1;
        eMmat_uu[ino][jno][1][2] = 0.0;
        ////
        eMmat_uu[ino][jno][2][0] = 0.0;
        eMmat_uu[ino][jno][2][1] = 0.0;
        eMmat_uu[ino][jno][2][2] = dtmp1;
      }
      eMmat_uu[ino][ino][0][0] += dtmp1;
      eMmat_uu[ino][ino][1][1] += dtmp1;
      eMmat_uu[ino][ino][2][2] += dtmp1;
    }
  }
  // supg for inertia term
  for(int jno=0;jno<nno;jno++){
    double tmp_vec[ndim] = { 0.0, 0.0, 0.0 };
    for(unsigned int kno=0;kno<nno;kno++){
      tmp_vec[0] += velo[kno][0];
      tmp_vec[1] += velo[kno][1];
      tmp_vec[2] += velo[kno][2];
    }
    tmp_vec[0] += velo[jno][0];
    tmp_vec[1] += velo[jno][1];
    tmp_vec[2] += velo[jno][2];
    for(int ino=0;ino<nno;ino++){
      const double dtmp1 = (dldx[ino][0]*tmp_vec[0]+dldx[ino][1]*tmp_vec[1]+dldx[ino][2]*tmp_vec[2])*rho*tau*vol*0.05;
      eMmat_uu[ino][jno][0][0] += dtmp1;
      eMmat_uu[ino][jno][1][1] += dtmp1;
      eMmat_uu[ino][jno][2][2] += dtmp1;
    }
  }
  
  // PSPG for inertia term
  double eMmat_pu[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eMmat_pu[ino][jno][0] = -tau*vol*0.25*dldx[ino][0];
      eMmat_pu[ino][jno][1] = -tau*vol*0.25*dldx[ino][1];
      eMmat_pu[ino][jno][2] = -tau*vol*0.25*dldx[ino][2];
    }
  }
  
  // eternal force term
  for(int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*rho*g_x*0.25;
    eres_u[ino][1] = vol*rho*g_y*0.25;
    eres_u[ino][2] = vol*rho*g_z*0.25;
  }
  for(int ino=0;ino<nno;ino++){
    eres_p[ino] = 0.0;
  }
  
   // SUPG external force
//   for(unsigned int ino=0;ino<nno;ino++){
//			const double dtmp1 = area*tau*rho*(ave_velo[0]*dldx[ino][0]+ave_velo[1]*dldx[ino][1]);
//			eres_u[ino][0] += dtmp1*g_x;
//			eres_u[ino][1] += dtmp1*g_y;
//   }
  
  
   // PSPG external force
//   for(unsigned int ino=0;ino<nno;ino++){
//			eres_p[ino] += area*tau*(dldx[ino][0]*g_x+dldx[ino][1]*g_y);
//   }
  
  //////////////////////////////////////////////////////////////////////
  
  {
    double dtmp1 = gamma*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
    }
    for(int i=0;i<nno*nno*ndim;i++){
      (&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
    }
    for(int i=0;i<nno*nno;i++){
      (&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_u[ino][0]
      -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_uu[ino][jno][0][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_uu[ino][jno][0][0]*acc[jno][0]
      +  eMmat_uu[ino][jno][0][1]*acc[jno][1]
      +  eMmat_uu[ino][jno][0][2]*acc[jno][2]
      +  eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
      eres_u[ino][1]
      -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_uu[ino][jno][1][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_uu[ino][jno][1][0]*acc[jno][0]
      +  eMmat_uu[ino][jno][1][1]*acc[jno][1]
      +  eMmat_uu[ino][jno][1][2]*acc[jno][2]
      +  eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
      eres_u[ino][2]
      -= eCmat_uu[ino][jno][2][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_uu[ino][jno][2][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_uu[ino][jno][2][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_uu[ino][jno][2][0]*acc[jno][0]
      +  eMmat_uu[ino][jno][2][1]*acc[jno][1]
      +  eMmat_uu[ino][jno][2][2]*acc[jno][2]
      +  eCmat_up[ino][jno][2]*(press[jno]+dt*apress[jno]);
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres_p[ino]
      -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0])
      +  eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1])
      +  eCmat_pu[ino][jno][2]*(velo[jno][2]+dt*acc[jno][2])
      +  eMmat_pu[ino][jno][0]*acc[jno][0]
      +  eMmat_pu[ino][jno][1]*acc[jno][1]
      +  eMmat_pu[ino][jno][2]*acc[jno][2]
      +  eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
    }
  }
}


void MakeMat_NavierStokes3D_Dynamic_P1
(double myu, double rho, double g_x, double g_y, double g_z,
 const double dt_timestep, const double gamma_newmark,
 const double coords[4][3],
 const double velo_press[4][4], const double acc_apress[4][4],
 double emat[4][4][4][4],
 double eres[4][4])
{
  const int nno = 4;
  const double velo[4][3] = {
    {velo_press[0][0],velo_press[0][1],velo_press[0][2]},
    {velo_press[1][0],velo_press[1][1],velo_press[1][2]},
    {velo_press[2][0],velo_press[2][1],velo_press[2][2]},
    {velo_press[3][0],velo_press[3][1],velo_press[3][2]} };
  const double press[4] = { velo_press[0][3], velo_press[1][3], velo_press[2][3], velo_press[3][3] };
  const double acc[4][3] = {
    {acc_apress[0][0],acc_apress[0][1],acc_apress[0][2]},
    {acc_apress[1][0],acc_apress[1][1],acc_apress[1][2]},
    {acc_apress[2][0],acc_apress[2][1],acc_apress[2][2]},
    {acc_apress[3][0],acc_apress[3][1],acc_apress[3][2]} };
  const double apress[4] = { acc_apress[0][3], acc_apress[1][3], acc_apress[2][3], acc_apress[3][3] };
  ////
  double emat_uu[4][4][3][3], emat_up[4][4][3], emat_pu[4][4][3], emat_pp[4][4];
  double eres_u[4][3], eres_p[4];
  MakeMat_NavierStokes3D_Dynamic_Newmark_P1P1(rho, myu, g_x, g_y,g_z,
                                              dt_timestep, gamma_newmark,
                                              coords,
                                              velo, press, acc, apress,
                                              emat_uu, emat_up, emat_pu, emat_pp,
                                              eres_u, eres_p);
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
      emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
      emat[ino][jno][0][2] = emat_uu[ino][jno][0][2];
      ////
      emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
      emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
      emat[ino][jno][1][2] = emat_uu[ino][jno][1][2];
      ////
      emat[ino][jno][2][0] = emat_uu[ino][jno][2][0];
      emat[ino][jno][2][1] = emat_uu[ino][jno][2][1];
      emat[ino][jno][2][2] = emat_uu[ino][jno][2][2];
      ////
      emat[ino][jno][0][3] = emat_up[ino][jno][0];
      emat[ino][jno][1][3] = emat_up[ino][jno][1];
      emat[ino][jno][2][3] = emat_up[ino][jno][2];
      ////
      emat[ino][jno][3][0] = emat_pu[ino][jno][0];
      emat[ino][jno][3][1] = emat_pu[ino][jno][1];
      emat[ino][jno][3][2] = emat_pu[ino][jno][2];
      ////
      emat[ino][jno][3][3] = emat_pp[ino][jno];
    }
  }
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = eres_u[ino][0];
    eres[ino][1] = eres_u[ino][1];
    eres[ino][2] = eres_u[ino][2];
    eres[ino][3] = eres_p[ino];
  }
}




/////////////////////////////





