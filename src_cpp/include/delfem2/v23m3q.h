/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef V23M3Q_H
#define V23M3Q_H

#include <stdio.h>

#include "vec2.h"
#include "vec3.h"
#include "mat3.h"
#include "quat.h"

CVector3 operator* (const CVector3& v, const CMatrix3& m);
CVector3 operator* (const CMatrix3& m, const CVector3& v);

// ---------------------------------------------

CVector2 screenXYProjection(const CVector3& v,
                            const float* mMV,
                            const float* mPj);


// --------------------------------------------

CMatrix3 Mat3(const CVector3& vec0,
              const CVector3& vec1,
              const CVector3& vec2);
CMatrix3 Mat3_RotCartesian(const CVector3& vec0);

CMatrix3 Mat3_OuterProduct(const CVector3& vec0,
                      const CVector3& vec1 );
CMatrix3 Mat3_Spin(const CVector3& vec0);
CMatrix3 Mat3_ParallelTransport(const CVector3& p0,
                                const CVector3& p1,
                                const CVector3& q0,
                                const CVector3& q1);
CMatrix3 Mat3_MinimumRotation(const CVector3& V,
                              const CVector3& v);
CMatrix3 Mat3_IrotTri(const CVector3& d0,
                      const CVector3& d1,
                      const CVector3& d2);
CMatrix3 Mat3_IrotTriSolid(const CVector3& d0,
                           const CVector3& d1,
                           const CVector3& d2);
CMatrix3 Mat3_IrotLineSeg(const CVector3& d0,
                          const CVector3& d1);
CMatrix3 Mat3_IrotPoint(const CVector3& d0);

void Mat4_MatTransl(double m[16],
                    const CMatrix3& mat, const CVector3& trans);
void Mat4_ScaleMatTransl(double m[16],
                         double scale, const CMatrix3& mat, const CVector3& trans);

int PickHandlerRotation_PosQuat(const CVector3& src, const CVector3& dir,
                                const CVector3& pos, const double quat[4], double rad,
                                double tol);
int PickHandlerRotation_Mat4(const CVector3& src, const CVector3& dir,
                             const double mat[16], double rad,
                             double tol);
bool DragHandlerRot_PosQuat(double quat[4], int ielem,
                            const CVector2& sp0, const CVector2& sp1,
                            const CVector3& pos,
                            const float mMV[16], const float mPj[16]);
bool DragHandlerRot_Mat4(double quat[4], int ielem,
                         const CVector2& sp0, const CVector2& sp1, double mat[16],
                         const float mMV[16], const float mPj[16]);
CVector3 drag_AxisHandler(const CVector2& sp0,
                          const CVector2& sp1,
                          const CVector3& p,
                          const CVector3& axis,
                          double len,
                          const float* mMV,
                          const float* mPj);
bool isPickPoint(const CVector2& sp,
                 const CVector3& p,
                 const float* mMV,
                 const float* mPj,
                 double pick_tol);
bool isPick_AxisHandler(const CVector2& sp,
                        const CVector3& p,
                        const CVector3& axis,
                        double len,
                        const float* mMV,
                        const float* mPj,
                        double pick_tol);
bool isPickQuad(const CVector3& p0,
                const CVector3& p1,
                const CVector3& p2,
                const CVector3& p3,
                const CVector2& sp,
                const CVector3& pick_dir,
                const float mMV[16],
                const float mPj[16],
                double eps);

bool isPickCircle(const CVector3& axis,
                  const CVector3& org,
                  double rad,
                  const CVector3& src,
                  const CVector3& dir,
                  double pick_tol);
double DragCircle(const CVector2& sp0,
                  const CVector2& sp1,
                  const CVector3& p,
                  const CVector3& axis,
                  const float* mMV,
                  const float* mPj);


#endif /* vec23mat3quat_hpp */
