/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_V23M3Q_H
#define DFM2_V23M3Q_H

#include <stdio.h>
#include "vec2.h"
#include "vec3.h"
#include "mat3.h"
#include "quat.h"

namespace delfem2 {

CVec3d operator* (const CVec3d& v, const CMat3d& m);
CVec3d operator* (const CMat3d& m, const CVec3d& v);
  
CVec3d MatVec(const CMat3d& m, const CVec3d& vec0);
CVec3d MatVecTrans(const CMat3d& m, const CVec3d& vec0);

// ---------------------------------------------

CVec2d screenXYProjection(const CVec3d& v,
                          const float* mMV,
                          const float* mPj);

void SetProjection(CMat3d& m, const CVec3d& vec0);
void SetDiag(CMat3d& m, const CVec3d& d);
void SetRotMatrix_Cartesian(CMat3d& m, const CVec3d& v);
void SetSpinTensor(CMat3d& m, const CVec3d& vec0);
void SetOuterProduct(CMat3d& m,
                     const CVec3d& vec0,
                     const CVec3d& vec1 );
CVec3d GetSpinVector(const CMat3d& m);
CVec3d GetCartesianRotationVector(const CMat3d& m);

  
// --------------------------------------------

CMat3d Mat3(const CVec3d& vec0,
              const CVec3d& vec1,
              const CVec3d& vec2);
CMat3d Mat3(const CVec3d& vec0);
CMat3d Mat3(const CVec3d& vec0,
              const CVec3d& vec1);
CMat3d RotMatrix_Cartesian(const CVec3d& v);
CMat3d Mat3_RotCartesian(const CVec3d& vec0);

/**
 * @brief output outer product Vec0 * Vec1^T
 */
CMat3d Mat3_OuterProduct(
    const CVec3d& vec0,
    const CVec3d& vec1 );
CMat3d Mat3_Spin(const CVec3d& vec0);
CMat3d Mat3_ParallelTransport(const CVec3d& p0,
                                const CVec3d& p1,
                                const CVec3d& q0,
                                const CVec3d& q1);

/**
 * @brief 3x3 Rotation matrix to rotate V into v with minimum rotation angle
 * @param V (in) rotation from
 * @param v (in) rotation to
 */
template <typename REAL>
CMat3<REAL> Mat3_MinimumRotation(
      const CVec3<REAL>& V,
      const CVec3<REAL>& v);

// --------------------------------------
// below: inertia tensor

CMat3d Mat3_IrotTri(const CVec3d& d0,
                      const CVec3d& d1,
                      const CVec3d& d2);

/**
 * @brief moment of inertia triangle pyramid with vtx (origin,d0,d1,d2) volume_density = 1
 */
CMat3d Mat3_IrotTriSolid(const CVec3d& d0,
                           const CVec3d& d1,
                           const CVec3d& d2);
CMat3d Mat3_IrotLineSeg(const CVec3d& d0,
                          const CVec3d& d1);
CMat3d Mat3_IrotPoint(const CVec3d& d0);

// above: inertia tensor
// ---------------------------------------

CMat3d Mirror(const CVec3d& n);

void Mat4_MatTransl(double m[16],
                    const CMat3d& mat, const CVec3d& trans);
void Mat4_ScaleMatTransl(double m[16],
                         double scale, const CMat3d& mat, const CVec3d& trans);

int PickHandlerRotation_PosQuat(const CVec3d& src, const CVec3d& dir,
                                const CVec3d& pos, const double quat[4], double rad,
                                double tol);
int PickHandlerRotation_Mat4(const CVec3d& src, const CVec3d& dir,
                             const double mat[16], double rad,
                             double tol);
bool DragHandlerRot_PosQuat(double quat[4], int ielem,
                            const CVec2d& sp0, const CVec2d& sp1,
                            const CVec3d& pos,
                            const float mMV[16], const float mPj[16]);
bool DragHandlerRot_Mat4(double quat[4], int ielem,
                         const CVec2d& sp0, const CVec2d& sp1, double mat[16],
                         const float mMV[16], const float mPj[16]);
CVec3d drag_AxisHandler(const CVec2d& sp0,
                       const CVec2d& sp1,
                       const CVec3d& p,
                       const CVec3d& axis,
                       double len,
                       const float* mMV,
                       const float* mPj);
bool isPickPoint(const CVec2d& sp,
                 const CVec3d& p,
                 const float* mMV,
                 const float* mPj,
                 double pick_tol);
bool isPick_AxisHandler(const CVec2d& sp,
                        const CVec3d& p,
                        const CVec3d& axis,
                        double len,
                        const float* mMV,
                        const float* mPj,
                        double pick_tol);
bool isPickQuad(const CVec3d& p0,
                const CVec3d& p1,
                const CVec3d& p2,
                const CVec3d& p3,
                const CVec2d& sp,
                const CVec3d& pick_dir,
                const float mMV[16],
                const float mPj[16],
                double eps);

bool isPickCircle(const CVec3d& axis,
                  const CVec3d& org,
                  double rad,
                  const CVec3d& src,
                  const CVec3d& dir,
                  double pick_tol);
double DragCircle(const CVec2d& sp0,
                  const CVec2d& sp1,
                  const CVec3d& p,
                  const CVec3d& axis,
                  const float* mMV,
                  const float* mPj);
  
bool isPickCircle(const CVec2d& sp,
                  const CVec3d& p,
                  const CVec3d& axis,
                  double r,
                  const float* mMV,
                  const float* mPj,
                  double pick_tol);  
}

#endif /* vec23mat3quat_hpp */
