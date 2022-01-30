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

// -------------------------------------------------------
// vec2


// --------------------------------------------------------
// vec3

DFM2_INLINE CVec2d screenXYProjection(
    const CVec3d &v,
    const float *mMV,
    const float *mPj);


// -------------
// vec3 and vec2



// ---------------------------------------------------
// matrix 3

template <typename T>
CVec3<T> operator*(
    const CVec3<T> &v,
    const CMat3<T> &m);

template <typename T>
CVec3<T> operator*(
    const CMat3<T> &m,
    const CVec3<T> &v);

template<typename T>
CVec3<T> MatVec(
    const CMat3<T> &m,
    const CVec3<T> &vec0);

template<typename T>
CVec3<T> MatVecTrans(
    const CMat3<T> &m,
    const CVec3<T> &vec0);

void SetProjection(CMat3d &m, const CVec3d &vec0);

DFM2_INLINE void SetDiag(
    CMat3d &m,
    const CVec3d &d);

DFM2_INLINE void SetRotMatrix_Cartesian(
    CMat3d &m,
    const CVec3d &v);

void SetSpinTensor(
    CMat3d &m,
    const CVec3d &vec0);

void SetOuterProduct(
    CMat3d &m,
    const CVec3d &vec0,
    const CVec3d &vec1);

CVec3d GetSpinVector(
    const CMat3d &m);

CVec3d GetCartesianRotationVector(
    const CMat3d &m);

template <typename VEC, typename REAL = typename VEC::Scalar>
std::array<REAL,9> Mat3_From3Bases(
    const VEC &vec0,
    const VEC &vec1,
    const VEC &vec2)
{
  return  {
    vec0.x, vec1.x, vec2.x,
    vec0.y, vec1.y, vec2.y,
    vec0.z, vec1.z, vec2.z};
}

CMat3d Mat3_FromCartesianRotationVector(
    const CVec3d &vec0);

/**
 * @brief output outer product Vec0 * Vec1^T
 */
CMat3d Mat3_OuterProduct(
    const CVec3d &vec0,
    const CVec3d &vec1);

CMat3d Mat3_Spin(
    const CVec3d &vec0);

void Mat4_MatTransl(
    double m[16],
    const CMat3d &mat,
    const CVec3d &trans);

DFM2_INLINE void Mat4_ScaleMatTransl(
    double m[16],
    double scale,
    const CMat3d &mat,
    const CVec3d &trans);


// ----------------------
// below: inertia tensor

CMat3d Mat3_IrotTri(
    const CVec3d &d0,
    const CVec3d &d1,
    const CVec3d &d2);

/**
 * @brief moment of inertia triangle pyramid with vtx (origin,d0,d1,d2) volume_density = 1
 */
CMat3d Mat3_IrotTriSolid(
    const CVec3d &d0,
    const CVec3d &d1,
    const CVec3d &d2);

CMat3d Mat3_IrotLineSeg(
    const CVec3d &d0,
    const CVec3d &d1);

CMat3d Mat3_IrotPoint(
    const CVec3d &d0);

// above: inertia tensor
// ----------------------

CMat3d Mirror(
    const CVec3d &n);

/**
 * @brief matrix for two cross products
 * @details Ma = v^(v^a)
 */
CMat3d Mat3_CrossCross(
    const CVec3d &v);

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
