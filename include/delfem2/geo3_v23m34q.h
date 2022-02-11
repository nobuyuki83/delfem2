/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @discussion this file may be split into v2 and v3m34q.
 * 2D application and 3D application typically don't have the same dependency
 */

#ifndef DFM2_V23M34Q_H
#define DFM2_V23M34Q_H

#include <cstdio>

#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2 {

// ---------------------------------------------------
// matrix 3

template<typename T>
CVec3<T> MatVec(
    const CMat3<T> &m,
    const CVec3<T> &vec0);

template<typename T>
CVec3<T> MatVecTrans(
    const CMat3<T> &m,
    const CVec3<T> &vec0);

DFM2_INLINE void SetDiag(
    CMat3d &m,
    const CVec3d &d);

CVec3d GetSpinVector(
    const CMat3d &m);

CVec3d GetCartesianRotationVector(
    const CMat3d &m);

void Mat4_MatTransl(
    double m[16],
    const CMat3d &mat,
    const CVec3d &trans);

DFM2_INLINE void Mat4_ScaleMatTransl(
    double m[16],
    double scale,
    const CMat3d &mat,
    const CVec3d &trans);

// mat3
// -----------------------------------------------------------
// quaternion

template<typename REAL>
DFM2_INLINE CVec3<REAL> operator*(
    const CQuat<REAL> &v,
    const CVec3<REAL> &m);

DFM2_INLINE CQuatd Quat_CartesianAngle(
    const CVec3d &p);

DFM2_INLINE bool Distortion_MappingTriangleFrom2To3Dim(
    double thresA,
    double thresE,
    unsigned int it0,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ,
    const std::vector<double> &aTexP);

} //


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo3_v23m34q.cpp"
#endif

#endif /* DFM2_V23M34Q_H */
