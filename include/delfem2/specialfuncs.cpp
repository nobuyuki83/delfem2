/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include "delfem2/specialfuncs.h"

//namespace dfm2 = delfem2;

//const double pi = 3.1415926535;
//const double invpi = 1.0/pi;

namespace delfem2 {
namespace specialfuncs {

DFM2_INLINE double normalize(double &x, double &z, double &y) {
  double r = sqrt(x * x + y * y + z * z);
  double invr = 1.0 / r;
  x *= invr;
  y *= invr;
  z *= invr;
  return r;
}

}
}

// size(Y) = (n+1)*(n+1)
// Y_l^m = Y[l*(l+1)+m]
// https://ciks.cbt.nist.gov/~garbocz/paper134/node10.html
void delfem2::makeArray_SphericalHarmonics
(double* Y,
 int n,
 double x, double y, double z)
{
  const double pi = 3.1415926535;
  const double invpi = 1.0/pi;
  std::complex<double> ep = std::complex<double>(+x,y);
  ////
  { // 0
    Y[ 0] =+0.5*sqrt(invpi);
  }
  if( n == 0 ) return;
  double r1 = real(ep);
  double i1 = imag(ep);
  { // 1
    double v1=-0.5*sqrt(1.5*invpi);
    Y[ 1] =v1*i1;
    Y[ 2] =+0.5*sqrt(3.0*invpi)*z;
    Y[ 3] =v1*r1;
  }
  if( n == 1 ) return;
  std::complex<double> ep2 = ep*ep;
  double r2 = real(ep2);
  double i2 = imag(ep2);
  double z2 = z*z;
  { // 2
    double v1=-0.50*sqrt(7.5*invpi)*z;
    double v2=+0.25*sqrt(7.5*invpi);
    Y[ 4] =v2*i2;
    Y[ 5] =v1*i1;
    Y[ 6] =+0.25*sqrt(5.0*invpi)*(2*z*z-x*x-y*y);
    Y[ 7] =v1*r1;
    Y[ 8] =v2*r2;
  }
  if( n == 2 ) return;
  double r3 = real(ep2*ep);
  double i3 = imag(ep2*ep);
  { // 3
    double v1=-0.125*sqrt(21*invpi)*(4*z*z-x*x-y*y);
    double v2=+0.250*sqrt(52.5*invpi)*z;
    double v3=-0.125*sqrt(35*invpi);
    Y[ 9] = v3*i3;
    Y[10] = v2*i2;
    Y[11] = v1*i1;
    Y[12] =+0.250*sqrt(7*invpi)*z*(-3*x*x-3*y*y+2*z*z);
    Y[13] = v1*r1;
    Y[14] = v2*r2;
    Y[15] = v3*r3;
  }
  if( n == 3 ) return;
  std::complex<double> ep4 = ep2*ep2;
  double r4 = real(ep4);
  double i4 = imag(ep4);
  double z4 = z2*z2;
  { // 4
    double v1=-3.0/8.00*sqrt(5.0*invpi)*z*(7.0*z2-3);
    double v2=+3.0/8.00*sqrt(5*0.5*invpi)*(7.0*z2-1);
    double v3=-3.0/8.00*sqrt(35*invpi)*z;
    double v4=+3.0/16.0*sqrt(35*0.5*invpi);
    Y[16] = v4*i4;
    Y[17] = v3*i3;
    Y[18] = v2*i2;
    Y[19] = v1*i1;
    Y[20]=+3.0/16.0*sqrt(invpi)*(35*z4-30*z2+3);
    Y[21] = v1*r1;
    Y[22] = v2*r2;
    Y[23] = v3*r3;
    Y[24] = v4*r4;
  }
  if( n == 4 ) return;
  double r5 = real(ep4*ep);
  double i5 = imag(ep4*ep);
  { // 5
    double v1=-1.0/16.0*sqrt(82.5*invpi)*(21*z4-14*z2+1);
    double v2=+1.0/8.00*sqrt(577.5*invpi)*z*(3*z2-1);
    double v3=-1.0/32.0*sqrt(385*invpi)*(9*z2-1);
    double v4=+3.0/16.0*sqrt(192.5*invpi)*z;
    double v5=-3.0/32.0*sqrt(77*invpi);
    Y[25] = v5*i5;
    Y[26] = v4*i4;
    Y[27] = v3*i3;
    Y[28] = v2*i2;
    Y[29] = v1*i1;
    Y[30]=+1.0/16.0*sqrt(11*invpi)*z*(63*z4-70*z2+15);
    Y[31] = v1*r1;
    Y[32] = v2*r2;
    Y[33] = v3*r3;
    Y[34] = v4*r4;
    Y[35] = v5*r5;
  }
  if( n == 5 ) return;
  double r6 = real(ep4*ep2);
  double i6 = imag(ep4*ep2);
  { // 6
    double v1=-1.0/16.0*sqrt(273*0.5*invpi)*z*(33*z4-30*z2+5);
    double v2=+1.0/64.0*sqrt(1365*invpi)*(33*z4-18*z2+1);
    double v3=-1.0/32.0*sqrt(1365*invpi)*z*(11*z2-3);
    double v4=+3.0/32.0*sqrt(91*0.5*invpi)*(11*z2-1);
    double v5=-3.0/32.0*sqrt(1001*invpi)*z;
    double v6=+1.0/64.0*sqrt(3003*invpi);
    Y[36] = v6*i6;
    Y[37] = v5*i5;
    Y[38] = v4*i4;
    Y[39] = v3*i3;
    Y[40] = v2*i2;
    Y[41] = v1*i1;
    Y[42]=+1.0/32.0*sqrt(13*invpi)*(231*z4*z2-315*z4+105*z2-5);
    Y[43] = v1*r1;
    Y[44] = v2*r2;
    Y[45] = v3*r3;
    Y[46] = v4*r4;
    Y[47] = v5*r5;
    Y[48] = v6*r6;
  }
  if( n == 6 ) return;
  double r7 = real(ep4*ep2*ep);
  double i7 = imag(ep4*ep2*ep);
  { // 7
    double v1=-1.0/64.0*sqrt(105*0.5*invpi)*(429*z4*z2-495*z4+135*z2-5);
    double v2=+3.0/64.0*sqrt(35*invpi)*(143*z4*z-110*z2*z+15*z);
    double v3=-3.0/64.0*sqrt(35*0.5*invpi)*(143*z4-66*z2+3);
    double v4=+3.0/32.0*sqrt(385*0.5*invpi)*(13*z2*z-3*z);
    double v5=-3.0/64.0*sqrt(385*0.5*invpi)*(13*z2-1);
    double v6=+3.0/64.0*sqrt(5005*invpi)*z;
    double v7=-3.0/128.0*sqrt(1430*invpi);
    Y[49] = v7*i7;
    Y[50] = v6*i6;
    Y[51] = v5*i5;
    Y[52] = v4*i4;
    Y[53] = v3*i3;
    Y[54] = v2*i2;
    Y[55] = v1*i1;
    Y[56]=+1.0/32.0*sqrt(15*invpi)*(429*z4*z2*z-693*z4*z+315*z2*z-35*z);
    Y[57] = v1*r1;
    Y[58] = v2*r2;
    Y[59] = v3*r3;
    Y[60] = v4*r4;
    Y[61] = v5*r5;
    Y[62] = v6*r6;
    Y[63] = v7*r7;
  }
  if( n == 7 ) return;
  std::complex<double> ep8 = ep4*ep4;
  double r8 = real(ep8);
  double i8 = imag(ep8);
  double z8 = z4*z4;  
  {  // 8
    double v1=-3.0/64.00*sqrt(17*0.5*invpi)*(715*z4*z2*z-1001*z4*z+385*z2*z-35*z);
    double v2=+3.0/128.0*sqrt(595*invpi)*(143*z4*z2-143*z4+33*z2-1);
    double v3=-1.0/64.00*sqrt(19635*0.5*invpi)*(39*z4*z-26*z2*z+3*z);
    double v4=+3.0/128.0*sqrt(1309*0.5*invpi)*(65*z4-26*z2+1);
    double v5=-3.0/64.00*sqrt(17017*0.5*invpi)*(5*z2*z-z);
    double v6=+1.0/128.0*sqrt(7293*invpi)*(15*z2-1);
    double v7=-3.0/64.00*sqrt(12155*0.5*invpi)*z;
    double v8=+3.0/256.0*sqrt(12155*0.5*invpi);
    Y[64] = v8*i8;
    Y[65] = v7*i7;
    Y[66] = v6*i6;
    Y[67] = v5*i5;
    Y[68] = v4*i4;
    Y[69] = v3*i3;
    Y[70] = v2*i2;
    Y[71] = v1*i1;
    Y[72]=+1.0/256.0*sqrt(17*invpi)*(6435*z8-12012*z4*z2+6930*z4-1260*z2+35);
    Y[73] = v1*r1;
    Y[74] = v2*r2;
    Y[75] = v3*r3;
    Y[76] = v4*r4;
    Y[77] = v5*r5;
    Y[78] = v6*r6;
    Y[79] = v7*r7;
    Y[80] = v8*r8;
  }
  if( n == 8 ) return;
  double r9 = real(ep8*ep);
  double i9 = imag(ep8*ep);
  { // 9
    double v1=-3.0/256*sqrt(95*0.5*invpi)*(2431*z8-4004*z4*z2+2002*z4-308*z2+7); //1
    double v2=+3.0/128*sqrt(1045*invpi)*z*(221*z4*z2-273*z4+91*z2-7); //2
    double v3=-1.0/256*sqrt(21945*invpi)*(221*z4*z2-195*z4+39*z2-1); //3
    double v4=+3.0/256*sqrt(95095*2*invpi)*z*(17*z4-10*z2+1); //4 *
    double v5=-3.0/256*sqrt(2717*invpi)*(85*z4-30*z2+1); //5
    double v6=+1.0/128*sqrt(40755*invpi)*z*(17*z2-3); //6
    double v7=-3.0/512*sqrt(13585*invpi)*(17*z2-1); //7
    double v8=-3.0/512*sqrt(230945*2*invpi)*z; //7
    double v9=-1.0/512*sqrt(230945*invpi); //7
    Y[81] = v9*i9;
    Y[82] = v8*i8;
    Y[83] = v7*i7;
    Y[84] = v6*i6;
    Y[85] = v5*i5;
    Y[86] = v4*i4;
    Y[87] = v3*i3;
    Y[88] = v2*i2;
    Y[89] = v1*i1;
    Y[90]=+1.0/256*sqrt(19*invpi)*z*(12155*z8-25740*z4*z2+18018*z4-4620*z2+315); // 0
    Y[91] = v1*r1;
    Y[92] = v2*r2;
    Y[93] = v3*r3;
    Y[94] = v4*r4;
    Y[95] = v5*r5;
    Y[96] = v6*r6;
    Y[97] = v7*r7;
    Y[98] = v8*r8;
    Y[99] = v9*r9;
  }
  if (n==9) return;
}



void delfem2::makeArray_CoeffSphericalHankelFirst(int n, std::complex<double>* h, double x)
{
  std::complex<double> i = std::complex<double>(0, 1);
  const double y = 1.0/x;
  //  COMPLEX eixy = y*exp(COMPLEX(0, x));
  
  h[0] = -i;
  if (n==0) return;
  
  h[1] = -(1.0+i*y);
  if (n==1) return;
  
  const double y2 = y*y;
  {
    double rl = 1.0-3*y2;
    double ig = 3*y;
    h[2] = (i*rl-ig);
  }
  if (n==2) return;
  
  const double y3 = y2*y;
  {
    double rl = 1  -15*y2;
    double ig = 6*y-15*y3;
    h[3] = (rl+i*ig);
  }
  if (n==3) return;
  
  const double y4 = y2*y2;
  {
    double rl = 1   - 45*y2+105*y4;
    double ig = 10*y-105*y3;
    h[4] = -(i*rl-ig);
  }
  if (n==4) return;
  
  const double y5 = y2*y3;
  {
    double rl = 1   -105*y2+945*y4;
    double ig = 15*y-420*y3+945*y5;
    h[5] = -(rl+i*ig);
  }
  if (n==5) return;
  
  const double y6 = y3*y3;
  {
    double rl = 1.0-210*y2+4725*y4-10395*y6;
    double ig = 21*y-1260*y3+10395*y5;
    //    double ig = 21*y*(1-60*y2+495*y4);
    //    double ig = 21*y*(1-15*y2*(4-33*y2));
    //    h[6] = i*eixy*(rl+i*ig);
    h[6] = (i*rl-ig);
  }
  if (n==6) return;
  
  const double y7 = y3*y4;
  {
    double rl = 1.0 - 378*y2+17325*y4-135135*y6;
    double ig = 28*y-3150*y3+62370*y5-135135*y7;
    h[7] = (rl+i*ig);
  }
  if (n==7) return;
  
  const double y8 = y4*y4;
  {
    double rl = 1.0 - 630*y2+ 51975*y4- 945945*y6+2027025*y8;
    double ig = 36*y-6930*y3+270270*y5-2027025*y7;
    h[8] = -(i*rl-ig);
  }
  if (n==8) return;
}


/*
void randomUnitVector(double r[3], double s){
  for(int itr=0;itr<100;itr++){
    r[0] = s*(2*(double)rand()/(RAND_MAX+1.0)-1);
    r[1] = s*(2*(double)rand()/(RAND_MAX+1.0)-1);
    r[2] = s*(2*(double)rand()/(RAND_MAX+1.0)-1);
    double l = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    if( (l <= 1 || itr==99) && l > 1.0e-5 ){
      r[0] /= l;
      r[1] /= l;
      r[2] /= l;
      return;
    }
  }
  r[0] = s;
  r[1] = 0;
  r[2] = 0;
}

COMPLEX evaluateSHDecomposition
(double phase,
 double x, double y, double z,
 const std::vector<COMPLEX>& aDecmpSH)
{
  COMPLEX c(0,0);
  for(int l=0;l<nSphericalHarmonics;l++){
    for(int m=-l;m<l+1;m++){
      int i = l*l+l+m;
      c += aDecmpSH[i]*sphericalHarmonics(l,m,x,y,z);
    }
  }
  return c;
}

double evaluateSHRealDecomposition
(double x, double y, double z,
 const std::vector<double>& aDecmpSHReal)
{
  double d = 0;
  for(int l=0;l<nSphericalHarmonics;l++){
    for(int m=0;m<l+1;m++){
      int i = l*(l+1)/2+m;
      d += aDecmpSHReal[i]*sphericalHarmonicsReal(l,m,x,y,z);
    }
  }
  return d;
}

void decomposeSH
(std::vector<COMPLEX>& aDecmpSH,
 double (*value)(double,double,double) )
{
  std::cout << "MakeDecomposition" << std::endl;
  aDecmpSH.resize(nSphericalHarmonics*nSphericalHarmonics);
  for(int ish=0;ish<nSphericalHarmonics*nSphericalHarmonics;ish++){
    aDecmpSH[ish] = COMPLEX(0,0);
  }
  const int n = 10000;
  for(int i=0;i<n;i++){
    double r[3]; randomUnitVector(r,1);
    for(int l=0;l<nSphericalHarmonics;l++){
      for(int m=-l;m<l+1;m++){
        int i = l*l+l+m;
        COMPLEX c = sphericalHarmonics(l,m, r[0],r[1],r[2]);
        double g = value(r[0], r[1], r[2]);
        aDecmpSH[i] += c*g;
      }
    }
  }
  for(int ish=0;ish<nSphericalHarmonics*nSphericalHarmonics;ish++){
    aDecmpSH[ish] *= 4*3.14159/n;
  }
}

void decomposeSHReal_Radom
(std::vector<double>& aDecmpSHReal,
 double (*value)(double,double,double))
{
  std::cout << "MakeDecomposition" << std::endl;
  const int nsh = nSphericalHarmonicsReal*(nSphericalHarmonicsReal+1)/2;
  aDecmpSHReal.resize(nsh);
  for(int ish=0;ish<nsh;ish++){ aDecmpSHReal[ish] = 0; }
  const int n = 100000;
  for(int i=0;i<n;i++){
    double r[3]; randomUnitVector(r,1);
    for(int l=0;l<nSphericalHarmonicsReal;l++){
      for(int m=0;m<l+1;m++){
        int i = l*(l+1)/2+m;
        double d = sphericalHarmonicsReal(l,m, r[0],r[1],r[2]);
        if( m != 0){ d *= 2; }
        double g = value(r[0], r[1], r[2]);
        aDecmpSHReal[i] += d*g;
      }
    }
  }
  for(int ish=0;ish<nsh;ish++){
    aDecmpSHReal[ish] *= 4*3.14159/n;
  }
}

double sphericalConvolution
(int nla, int nlo,
 double (*value1)(double,double,double),
 double (*value2)(double,double,double) )
{
  double dla = 3.14159*1/nla;
  double dlo = 3.14159*2/nlo;
  double v = 0.0;
  for(int ila=0;ila<nla+1;ila++){
    double y0 = cos(dla*ila);
    double r0 = sin(dla*ila);
    for(int ilo=0;ilo<nlo;ilo++){
      double x0 = sin(dlo*ilo)*r0;
      double z0 = cos(dlo*ilo)*r0;
      double aj = r0*dlo*dla;
      if( ila == 0 || ila == nla ){
        aj = 0.25*dla*dla/nlo;
        r0 = dla*0.25;
        x0 = sin(dlo*ilo)*r0;
        z0 = cos(dlo*ilo)*r0;
        if( ilo == 0 ){ y0 = +cos(dla*0.25); }
        else{           y0 = -cos(dla*0.25); }
        double invl = 1.0/sqrt(x0*x0+y0*y0+z0*z0);
        x0 *= invl;
        y0 *= invl;
        z0 *= invl;
      }
      v += aj*value1(x0,y0,z0)*value2(x0,y0,z0);
    }
  }
  return v;
}

void decomposeSHReal_Sphere
(std::vector<double>& aDecmpSHReal,
 double (*value)(double,double,double))
{
  const int nsh = nSphericalHarmonicsReal*(nSphericalHarmonicsReal+1)/2;
  aDecmpSHReal.resize(nsh);
  for(int ish=0;ish<nsh;ish++){ aDecmpSHReal[ish] = 0; }
  ////
  int nla =  64;
  int nlo =  32;
  for(int l=0;l<nSphericalHarmonicsReal;l++){
    for(int m=0;m<l+1;m++){
      int ish = l*(l+1)/2+m;
      double v = sphericalConvolution(nla, nlo, arraySphericalHarmonicsReal[ish],value);
      if( m != 0 ){ v *= 2; }
      aDecmpSHReal[ish] += v;
    }
  }
}
 */
