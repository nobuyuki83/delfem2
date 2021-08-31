/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_NOISE_H
#define DFM2_NOISE_H

#include <cassert>
#include <vector>

namespace delfem2 {

void Shuffle(std::vector<int>& X)
{
  const int N = (int)X.size();
  for(int i=0;i<N;i++){
    int j = rand()%N;
    int k = X[i];
    X[i] = X[j];
    X[j] = k;
  }
}

double grad2
(double x, double y,
 int igrad,
 const std::vector<double>& aGrad)
{
  const int ng = static_cast<int>(aGrad.size()/2);
  const int ig0 = igrad%ng;
  return aGrad[ig0*2+0]*x + aGrad[ig0*2+1]*y;
}

double grad3
(double x, double y, double z,
 int igrad,
 const std::vector<double>& aGrad)
{
  const int ng = static_cast<int>(aGrad.size()/3);
  const int ig0 = igrad%ng;
  return aGrad[ig0*3+0]*x + aGrad[ig0*3+1]*y + aGrad[ig0*3+2]*z;
}

template <typename REAL>
REAL fade(REAL t) {
  //  return t;
  return t*t*t*(t*(t*6 - 15) + 10);
}

double noise_perlin_2d
(double x, double y, int np,
 const std::vector<double>& aGrad,
 const std::vector<int>& aP)
{
  const int ix0 = (int)x % np; // ix [0,255]
  const int iy0 = (int)y % np; // iy [0,255]
  const int ix1 = (ix0+1) % np;
  const int iy1 = (iy0+1) % np;
  const double ux = x - (int)x; // [0,1]
  const double uy = y - (int)y; // [0,1]
  const int ig_aa = aP[aP[ix0]+iy0];
  const int ig_ba = aP[aP[ix1]+iy0];
  const int ig_ab = aP[aP[ix0]+iy1];
  const int ig_bb = aP[aP[ix1]+iy1];
  double gaa = grad2(  ux,  uy, ig_aa,aGrad);
  double gba = grad2(1-ux,  uy, ig_ba,aGrad);
  double gab = grad2(  ux,1-uy, ig_ab,aGrad);
  double gbb = grad2(1-ux,1-uy, ig_bb,aGrad);
  double fux = fade(ux);
  double fuy = fade(uy);
  //  std::cout << fux << " " << fuy << std::endl;
  double vxa = (1-fux)*gaa + fux*gba;
  double vxb = (1-fux)*gab + fux*gbb;
  return (1-fuy)*vxa + fuy*vxb;
}

double noise_perlin_3d(
	double x, 
	double y, 
	double z, 
	int np,
	const std::vector<double>& aGrad,
	const std::vector<int>& aP)
{
  assert( aP.size() > 2 );
  const int ix0 = (int)x % np; // ix [0,255]
  const int iy0 = (int)y % np; // iy [0,255]
  const int iz0 = (int)z % np; // iy [0,255]
  const int ix1 = (ix0+1) % np;
  const int iy1 = (iy0+1) % np;
  const int iz1 = (iz0+1) % np;
  const double ux = x - (int)x; // [0,1]
  const double uy = y - (int)y; // [0,1]
  const double uz = z - (int)z; // [0,1]
  const int ig_aaa = aP[aP[aP[ix0]+iy0]+iz0];
  const int ig_baa = aP[aP[aP[ix1]+iy0]+iz0];
  const int ig_aba = aP[aP[aP[ix0]+iy1]+iz0];
  const int ig_bba = aP[aP[aP[ix1]+iy1]+iz0];
  const int ig_aab = aP[aP[aP[ix0]+iy0]+iz1];
  const int ig_bab = aP[aP[aP[ix1]+iy0]+iz1];
  const int ig_abb = aP[aP[aP[ix0]+iy1]+iz1];
  const int ig_bbb = aP[aP[aP[ix1]+iy1]+iz1];
  double gaaa = grad3(  ux,  uy,  uz, ig_aaa,aGrad);
  double gbaa = grad3(1-ux,  uy,  uz, ig_baa,aGrad);
  double gaba = grad3(  ux,1-uy,  uz, ig_aba,aGrad);
  double gbba = grad3(1-ux,1-uy,  uz, ig_bba,aGrad);
  double gaab = grad3(  ux,  uy,1-uz, ig_aab,aGrad);
  double gbab = grad3(1-ux,  uy,1-uz, ig_bab,aGrad);
  double gabb = grad3(  ux,1-uy,1-uz, ig_abb,aGrad);
  double gbbb = grad3(1-ux,1-uy,1-uz, ig_bbb,aGrad);
  double fux = fade(ux);
  double fuy = fade(uy);
  double fuz = fade(uz);
  double vxaa = (1-fux)*gaaa + fux*gbaa;
  double vxba = (1-fux)*gaba + fux*gbba;
  double vxab = (1-fux)*gaab + fux*gbab;
  double vxbb = (1-fux)*gabb + fux*gbbb;
  double vxya = (1-fuy)*vxaa + fuy*vxba;
  double vxyb = (1-fuy)*vxab + fuy*vxbb;
  double vxyz = (1-fuz)*vxya + fuz*vxyb;
  return vxyz;
}

double noise_perlin_2d_oct
(double x, double y,
 int nrep,
 int noct, double persistence,
 const std::vector<double>& aGrad,
 const std::vector<int>& aP)
{
  double val = 0.;
  double mag = 1.;
  double freq = 1.0;
  for(int ioct=0;ioct<noct;++ioct){
    val += mag*noise_perlin_2d(
		x*freq,y*freq, 
		static_cast<int>(nrep*freq), 
		aGrad,aP);
    freq *= 2.0;
    mag *= persistence;
  }
  return val;
}

double noise_perlin_3d_oct(
	double x, 
	double y, 
	double z,
	int nrep,
	int noct, 
	double persistence,
	const std::vector<double>& aGrad,
	const std::vector<int>& aP)
{
  double val = 0;
  double  mag = 1;
  double freq = 1;
  for(int ioct=0;ioct<noct;++ioct){
    val += mag*noise_perlin_3d(
		x*freq,y*freq,z*freq,
		static_cast<int>(nrep*freq), 
		aGrad,aP);
    freq *= 2.0;
    mag *= persistence;
  }
  return val;
}

void ComputePerlin(
    std::vector<double> &aV,
    unsigned int nW,
    unsigned int nH,
    int nrep,
    int noct,
    double persistance) {
  std::vector<int> aP;
  aP.resize(256);
  for (int i = 0; i < 256; ++i) { aP[i] = i; }
  Shuffle(aP);
  aP.resize(512);
  for (int i = 0; i < 256; ++i) { aP[256 + i] = i; }

  std::vector<double> aGrad = {-1, -1, -1, +1, +1, -1, +1, +1};

  aV.resize(nH * nW);
  for (unsigned int ih = 0; ih < nH; ++ih) {
    for (unsigned int iw = 0; iw < nW; ++iw) {
      double x = (double) iw / (nW) * nrep;
      double y = (double) ih / (nH) * nrep;
      aV[ih * nW + iw] = noise_perlin_2d_oct(x, y,
          nrep, noct, persistance,
          aGrad, aP);
    }
  }
}
  
}

#endif
