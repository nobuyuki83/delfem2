/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_GEO_CCD_H
#define DFM2_GEO_CCD_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template<typename T>
bool FindCoplanerInterp(
    double &r,
    const CVec3<T> &s0, const CVec3<T> &s1, const CVec3<T> &s2, const CVec3<T> &s3,
    const CVec3<T> &e0, const CVec3<T> &e1, const CVec3<T> &e2, const CVec3<T> &e3);

DFM2_INLINE double Nearest_LineSeg_LineSeg_CCD_Iteration(
  double p[3],
  const CVec3d &p0s,
  const CVec3d &p0e,
  const CVec3d &p1s,
  const CVec3d &p1e,
  const CVec3d &q0s,
  const CVec3d &q0e,
  const CVec3d &q1s,
  const CVec3d &q1e,
  unsigned int nitr);

template<typename T>
bool IsContact_FV_CCD2(
    int ino0, int ino1, int ino2, int ino3,
    const CVec3<T> &p0, const CVec3<T> &p1, const CVec3<T> &p2, const CVec3<T> &p3,
    const CVec3<T> &q0, const CVec3<T> &q1, const CVec3<T> &q2, const CVec3<T> &q3);


} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_ccd.cpp"
#endif

#endif // DFM2_GEO_CCD_H
