 /*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_GEODELAUNAY3_V3_H
#define DFM2_GEODELAUNAY3_V3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {
  
/**
 * @function check if Delaunay condition satisfied
 * @return
 * 0 : p3 is inside circum circle on the p0,p1,p2
 * 1 :       on
 * 2 :       outsdie
 */
template <typename T>
int DetDelaunay(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3);

template <typename T>
double SquareCircumradius(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3);

template <typename T>
CVec3<T> CircumCenter(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3);

template <typename T>
double Circumradius(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geodelaunay3_v3.cpp"
#endif


#endif // VEC3_H
