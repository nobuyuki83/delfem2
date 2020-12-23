/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/femem2.h"
#include <math.h>
#include <assert.h>
#include <complex>

// --------------------------------------------------------

namespace delfem2 {
namespace femem2 {


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

// area of a triangle
DFM2_INLINE double TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

// area coordinate inside a triangle
DFM2_INLINE void TriAreaCoord(double vc_p[],
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

}
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
  const double area = femem2::TriArea2D(p0, p1, p2);
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



DFM2_INLINE void delfem2::EMat_Poisson_Tri2D
(double eres[3], double emat[3][3],
 const double alpha, const double source,
 const double coords[3][2],
 const double value[3])
{
  const int nno = 3;
  const int ndim = 2;
  ////
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;
  for (int i = 0; i<9; ++i){ (&emat[0][0])[i] = 0.0; }
  const double area = femem2::TriArea2D(coords[0], coords[1], coords[2]);
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


DFM2_INLINE void delfem2::EMat_Helmholtz_Tri2D(
    std::complex<double> eres[3],
    std::complex<double> emat[3][3],
    const double wave_length,
    const double coords[3][2],
    const std::complex<double> value[3])
{
  const int nno = 3;
  const int ndim = 2;
  const double area = femem2::TriArea2D(coords[0],coords[1],coords[2]);
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


DFM2_INLINE void delfem2::EMat_Diffusion_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha,
    const double source,
    const double dt_timestep,
    const double gamma_newmark,
    const double rho,
    const double coords[3][2],
    const double value[3],
    const double velo[3])
{
  const int nno = 3;
  const int ndim = 2;
  
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;
  for (int i = 0; i<9; ++i){ (&emat[0][0])[i] = 0.0; }
  
  const double area = femem2::TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx,const_term,coords[0],coords[1],coords[2]);
  
  // --------------------------
  
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

DFM2_INLINE void delfem2::EMat_SolidStaticLinear_Tri2D(
    double eres[3][2],
    double emat[3][3][2][2],
    const double myu,
    const double lambda,
    const double rho,
    const double g_x,
    const double g_y,
    const double disp[3][2],
    const double coords[3][2])
{
  const int nno = 3;
  const int ndim = 2;
  
  const double area = femem2::TriArea2D(coords[0],coords[1],coords[2]);
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
  
  // ----------------------------------------------
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1];
      eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1];
    }
  }
}


DFM2_INLINE void delfem2::EMat_SolidDynamicLinear_Tri2D(
    double eres[3][2],
    double emat[3][3][2][2],
    const double myu, const double lambda,
    const double rho, const double g_x, const double g_y,
    const double dt_timestep, const double gamma_newmark,  const double beta_newmark,
    const double disp[3][2], const double velo[3][2], const double acc[3][2],
    const double coords[3][2],
    bool is_initial)
{
  const int nno = 3;
  const int ndim = 2;
  
	const double area = femem2::TriArea2D(coords[0],coords[1],coords[2]);
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

// above: solid linear
// ----------------------------------------
// below: stokes 2D

DFM2_INLINE void delfem2::MakeMat_Stokes2D_Static_P1P1(
    double alpha, double g_x, double g_y,
    const double coords[][2],
    const double velo[3][2], const double press[3],
    double emat_uu[][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3],
    double eres_u[][2], double eres_p[3])
{
  const unsigned int nno = 3;
  const unsigned int ndim = 2;
  
  const double area = femem2::TriArea2D(coords[0],coords[1],coords[2]);
  
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  
  // -------------------------
  
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

DFM2_INLINE void delfem2::EMat_Stokes2D_Static_P1(
    double alpha, double g_x, double g_y,
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


DFM2_INLINE void delfem2::EMat_Stokes2D_Dynamic_P1(
    double alpha, double rho, double g_x, double g_y,
    const double dt_timestep,
    const double gamma_newmark,
    const double coords[][2],
    const double velo_press[3][3],
    const double acc_apress[3][3],
    double emat[3][3][3][3],
    double eres[3][3])
{
  const int nno = 3;

  const double area = femem2::TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][2], const_term[nno];
  TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);
  
  double tau;
  {
    const double h = sqrt( area / 3.14 )*2;
    tau = -h*h/alpha*0.1;
  }
  
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

// above: stokes equation
// --------------------------
// below: navier-stokes equation

DFM2_INLINE void delfem2::MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1(
    double rho,
    double myu,
    double g_x,
    double g_y,
    double dt,
    double gamma,
    const double coords[][2],
    const double velo[][2],
    const double press[],
    const double acc[][2],
    const double apress[],
    double emat_uu[][3][2][2],
    double emat_up[][3][2],
    double emat_pu[][3][2],
    double emat_pp[][3],
    double eres_u[3][2],
    double eres_p[3])
{
  const int nno = 3;
  const int ndim = 2;
  
  const double area = femem2::TriArea2D(coords[0],coords[1],coords[2]);
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
      const double re_c = 0.5*norm_v*h*rho/myu;  // 0.5*norm_v*h*rho/myu;
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
  {  // advection
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
  {  // SUPG for advection term
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
  
  // ------------------------------
  double eMmat_uu[3][3][2][2];
  {  // add inertia term
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

DFM2_INLINE void delfem2::EMat_NavierStokes2D_Dynamic_P1(
    double myu, double rho, double g_x, double g_y,
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
  //
  double emat_uu[3][3][2][2], emat_up[3][3][2], emat_pu[3][3][2], emat_pp[3][3];
  double eres_u[3][2], eres_p[3];
  MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1(rho, myu, g_x, g_y,
                                              dt_timestep, gamma_newmark,
                                              coords,
                                              velo, press, acc, apress,
                                              emat_uu, emat_up, emat_pu, emat_pp,
                                              eres_u, eres_p);
  //
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




