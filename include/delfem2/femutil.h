/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMUTIL_H
#define DFM2_FEMUTIL_H

#include "delfem2/dfm2_inline.h"
#include <cassert>
#include <climits>

namespace delfem2 {

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

const static unsigned int NIntTetGauss[4] = {
    1, 4, 5, 16
};

const static double TetGauss[4][16][4] = {
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


namespace femutil {

// area of a triangle
DFM2_INLINE double TriArea2D(
    const double p0[],
    const double p1[],
    const double p2[]);

// area coordinate inside a triangle
DFM2_INLINE void TriAreaCoord(
    double vc_p[],
    const double p0[],
    const double p1[],
    const double p2[],
    const double pb[]);

DFM2_INLINE double Distance3D(
    const double p0[3],
    const double p1[3]);

DFM2_INLINE double TriArea3D(
    const double v1[3],
    const double v2[3],
    const double v3[3]);

DFM2_INLINE void UnitNormalAreaTri3D(
    double n[3],
    double& a,
    const double v1[3],
    const double v2[3],
    const double v3[3]);

DFM2_INLINE double Dot3D(
    const double a[],
    const double b[]);

DFM2_INLINE void Cross3D(
    double r[3],
    const double v1[3],
    const double v2[3]);

DFM2_INLINE double TetVolume3D(
    const double v1[3],
    const double v2[3],
    const double v3[3],
    const double v4[3]);

DFM2_INLINE void MatVec3(
    double y[3],
    const double m[9],
    const double x[3]);

DFM2_INLINE void MatMat3(
    double* C,
    const double* A,
    const double* B);

DFM2_INLINE void MatMatTrans3(
    double* C,
    const double* A,
    const double* B);


}

// =======================================================================

// derivative of a shape function of a triangle and constant compornent 
DFM2_INLINE void TriDlDx(
    double dldx[][2],
    double const_term[],
    const double p0[],
    const double p1[],
    const double p2[]);

DFM2_INLINE void TetDlDx(
    double dldx[][3],
    double a[],
    const double p0[],
    const double p1[],
    const double p2[],
    const double p3[]);

DFM2_INLINE void ShapeFunc_Hex8
    (const double& r0, const double& r1,	const double& r2,
     const double coords[][3],
     double& detjac,
     double dndx[][3],
     double an[] );

template <int nno, int ndim>
DFM2_INLINE void FetchData(
    double val_to[nno][ndim],
    const unsigned int* aIP,
    const double* val_from,
    int nstride=-1)
{
  if( nstride == -1 ){ nstride = ndim; }
  assert( nstride >= ndim );
  for(int ino=0;ino<nno;++ino){
    unsigned int ip = aIP[ino];
    for(int idim=0;idim<ndim;++idim){
      val_to[ino][idim] = val_from[ip*nstride+idim];
    }
  }
}

DFM2_INLINE void ddW_MassConsistentVal3D_Tet3D(
    double* eMmat,
    double rho, double vol,
    bool is_add,
    unsigned int nstride);

} // namespace delfem2


#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femutil.cpp"
#endif

  
#endif /* fem_ematrix_h */
