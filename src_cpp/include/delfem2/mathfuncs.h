/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef MATHFUNCS_H
#define MATHFUNCS_H

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>

#ifndef COMPLEX
typedef std::complex<double> COMPLEX;
#endif

// n <= 9
// size(Y) = (n+1)*(n+1)
// Y_l^m = Y[l*(l+1)+m]
void makeArray_SphericalHarmonics(double* a, int norder, double x, double y, double z);

void makeArray_CoeffSphericalHankelFirst(int n, COMPLEX* h, double x);


// size(h) = (n+1)
inline void makeArray_SphericalHankelFirst(int n, COMPLEX* h, double x){
  makeArray_CoeffSphericalHankelFirst(n, h, x);
  const double y = 1.0/x;
  COMPLEX eixy = y*exp(COMPLEX(0, x));
  for(int i=0;i<n+1;i++){
    h[i] *= eixy;
  }
}

inline void makeArray_Legendre(int n, double* P, double x)
{
  P[0] = 1;
  if (n==0) return; 

  P[1] = x;
  if (n==1) return;

  const double x2 = x*x;
  P[2] = 0.5*(3*x2-1);
  if (n==2) return;

  const double x3 = x*x2;
  P[3] = 0.5*x*(5*x2-3);
  if (n==3) return;

  const double x4 = x2*x2;
  P[4] = 0.125*(35*x4-30*x2+3);
  if (n==4) return;

//  const double x5 = x3*x2;
  P[5] = 0.125*x*(63*x4-70*x2+15);
  if (n==5) return;

  const double x6 = x3*x3;
  P[6] = 0.0625*(231*x6-315*x4+105*x2-5);
  if (n==6) return;

//  const double x7 = x3*x4;
  P[7] = 0.0625*x*(429*x6-693*x4+315*x2-35);
//  P[7] = 0.0625*x*(((13*x2-21)*11*x2+105)*3*x2-35);
  if (n==7) return;

  const double x8 = x4*x4;
  P[8] = 0.0078125*(6435*x8-12012*x6+6930*x4-1260*x2+35);
//  P[8] = 0.0078125*((2145*x6-4004*x4+2310*x2-420)*3*x2+35);
//  P[8] = 0.0078125*((((15*x2-28)*13*x2+210)*11*x2-420)*3*x2+35);
  if (n==8) return;
}

// size(Y) = (n+1)*(n+1)
// P_l^m = P[l*(l+1)+m]
inline void makeArrayAssociatedLegendre
(int n, double* P, double x)
{
  double y = 1.0-x*x;
  if (y<0.0){ y = 0.0; }
  y = sqrt(y);
  ////
  { // 0
    P[0] = 1;
  }
  if (n==0) return;
  { // 1
    P[2] = x;

    P[3] = -y;
    P[1] = -P[3];
  }
  if (n==1) return;
  { // 2    
    P[6] = 0.5*(3*x*x-1);

    P[7] = -3*x*y;
    P[5] = -P[7];

    P[8] = 3*y*y;
    P[4] = +P[8];
  }
  if (n==2) return;
  { // 3
    P[12] = 0.5*x*(5*x*x-3);

    P[13] = -1.5*(5*x*x-1)*y;
    P[11] = -P[13];

    P[14] = 15.0*x*y*y;
    P[10] = +P[14];

    P[15] = -15.0*y*y*y;
    P[9] = -P[15];
  }
  if (n==3) return;
}


////////////////////////////////////////////////////////////////////////////////////
/*
#include <vector>
double evaluateSHRealDecomposition
(double x, double y, double z,
 const std::vector<double>& aDecmpSHReal);

COMPLEX evaluateSHDecomposition
(double phase,
 double x, double y, double z,
 const std::vector<COMPLEX>& aDecmpSH);

void decomposeSH
(std::vector<COMPLEX>& aDecmpSH,
 double (*value)(double,double,double) );

void decomposeSHReal_Radom
(std::vector<double>& aDecmpSHReal,
 double (*value)(double,double,double));

void decomposeSHReal_Sphere
(std::vector<double>& aDecmpSHReal,
 double (*value)(double,double,double));

double sphericalConvolution
(int nla, int nlo,
 double (*value1)(double,double,double),
 double (*value2)(double,double,double) );
*/

#endif
