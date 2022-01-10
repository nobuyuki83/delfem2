/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/fem_discreteshell.h"

#include "delfem2/mat2.h"
#include "delfem2/geo3_v23m34q.h"

// =======================================


template<typename T>
DFM2_INLINE void delfem2::CdC_DiscreteShell(
    double &C,
    CVec3<T> dC[4],
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3) {
  const CVec3<T> v02 = p2 - p0;
  const CVec3<T> v03 = p3 - p0;
  const CVec3<T> v12 = p2 - p1;
  const CVec3<T> v13 = p3 - p1;
  const CVec3<T> v23 = p3 - p2;
  // ---
  const CVec3<T> A = v02.cross(v03);
  const CVec3<T> B = v13.cross(v12);
  const double lA = A.norm();
  const double lB = B.norm();
  const CVec3<T> a = A / lA;
  const CVec3<T> b = B / lB;
  const double ab = a.dot(b);
  //  C = acos(ab);
  C = ab - 1;
  const double sab = 1.0;//-1.0/sin(C);
  const CVec3<T> tmpBA = (b - a * (a.dot(b))) * (sab / lA);
  const CVec3<T> tmpAB = (a - b * (b.dot(a))) * (sab / lB);
  dC[0] = tmpBA.cross(v23);
  dC[1] = v23.cross(tmpAB);
  dC[2] = v03.cross(tmpBA) + tmpAB.cross(v13);
  dC[3] = tmpBA.cross(v02) + v12.cross(tmpAB);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CdC_DiscreteShell(
    double &,
    CVec3d dC[4],
    const CVec3d &,
    const CVec3d &,
    const CVec3d &,
    const CVec3d &);
#endif
