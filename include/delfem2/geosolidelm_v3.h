 /*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_GEOSOLIDELM_V3_H
#define DFM2_GEOSOLIDELM_V3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template <typename T>
T Volume_Tet3(
    const T v1[3],
    const T v2[3],
    const T v3[3],
    const T v4[3]);

template <typename T>
double Height(
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3,
    const CVec3<T>& v4);

// ---------------------------------------------------------

template <typename T>
bool barycentricCoord_Origin_Tet(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3);
  
template <typename T>
bool barycentricCoord_Origin_Pyramid(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4);
  
template <typename T>
bool barycentricCoord_Origin_Wedge(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5);

template <typename T>
CVec3<T> positionBarycentricCoord_Pyramid(
    double r0,
    double r1,
    double r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4);
  
template <typename T>
CVec3<T> positionBarycentricCoord_Wedge(
    double r0,
    double r1,
    double r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5);
  
template <typename T>
void iteration_barycentricCoord_Origin_Solid(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& q, // q=positionBarycentricCoord_Wedge
    const CVec3<T>& dpdr0,
    const CVec3<T>& dpdr1,
    const CVec3<T>& dpdr2,
    double damp=1.0);

  
// -----------------------------------------------

template <typename T>
T Volume_OrgTet(
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 );
  
template <typename T>
T Volume_Tet(
    const CVec3<T>& v0,
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 );

template <typename T>
double Volume_Pyramid(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4);

template <typename T>
double Volume_Wedge(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5);

// ---------------------------------------------

template <typename T>
double SqareLongestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

template <typename T>
double LongestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

template <typename T>
double SqareShortestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

template <typename T>
double ShortestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

// ----------------------------------------------------------
// here starts std::vector<CVector3>
  
template <typename T>
double Volume_Tet(
    int iv1,
    int iv2,
    int iv3,
    int iv4,
    const std::vector<CVec3<T> >& aPoint);

 template <typename T>
 double SolidAngleTri(
     const CVec3<T>& v1,
     const CVec3<T>& v2,
     const CVec3<T>& v3);


} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geosolidelm_v3.cpp"
#endif


#endif // VEC3_H
