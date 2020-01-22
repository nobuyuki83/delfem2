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

CVec3 operator* (const CVec3& v, const CMatrix3& m);
CVec3 operator* (const CMatrix3& m, const CVec3& v);
  
CVec3 MatVec(const CMatrix3& m, const CVec3& vec0);
CVec3 MatVecTrans(const CMatrix3& m, const CVec3& vec0);

// ---------------------------------------------

CVec2 screenXYProjection(const CVec3& v,
                            const float* mMV,
                            const float* mPj);

void SetProjection(CMatrix3& m, const CVec3& vec0);
void SetDiag(CMatrix3& m, const CVec3& d);
void SetRotMatrix_Cartesian(CMatrix3& m, const CVec3& v);
void SetSpinTensor(CMatrix3& m, const CVec3& vec0);
void SetOuterProduct(CMatrix3& m,
                     const CVec3& vec0,
                     const CVec3& vec1 );
CVec3 GetSpinVector(const CMatrix3& m);
CVec3 GetCartesianRotationVector(const CMatrix3& m);

  
// --------------------------------------------

CMatrix3 Mat3(const CVec3& vec0,
              const CVec3& vec1,
              const CVec3& vec2);
CMatrix3 Mat3(const CVec3& vec0);
CMatrix3 Mat3(const CVec3& vec0,
              const CVec3& vec1);
CMatrix3 RotMatrix_Cartesian(const CVec3& v);
CMatrix3 Mat3_RotCartesian(const CVec3& vec0);
CMatrix3 Mat3_OuterProduct(const CVec3& vec0,
                      const CVec3& vec1 );
CMatrix3 Mat3_Spin(const CVec3& vec0);
CMatrix3 Mat3_ParallelTransport(const CVec3& p0,
                                const CVec3& p1,
                                const CVec3& q0,
                                const CVec3& q1);
CMatrix3 Mat3_MinimumRotation(const CVec3& V,
                              const CVec3& v);
CMatrix3 Mat3_IrotTri(const CVec3& d0,
                      const CVec3& d1,
                      const CVec3& d2);
CMatrix3 Mat3_IrotTriSolid(const CVec3& d0,
                           const CVec3& d1,
                           const CVec3& d2);
CMatrix3 Mat3_IrotLineSeg(const CVec3& d0,
                          const CVec3& d1);
CMatrix3 Mat3_IrotPoint(const CVec3& d0);
CMatrix3 Mirror(const CVec3& n);

void Mat4_MatTransl(double m[16],
                    const CMatrix3& mat, const CVec3& trans);
void Mat4_ScaleMatTransl(double m[16],
                         double scale, const CMatrix3& mat, const CVec3& trans);

int PickHandlerRotation_PosQuat(const CVec3& src, const CVec3& dir,
                                const CVec3& pos, const double quat[4], double rad,
                                double tol);
int PickHandlerRotation_Mat4(const CVec3& src, const CVec3& dir,
                             const double mat[16], double rad,
                             double tol);
bool DragHandlerRot_PosQuat(double quat[4], int ielem,
                            const CVec2& sp0, const CVec2& sp1,
                            const CVec3& pos,
                            const float mMV[16], const float mPj[16]);
bool DragHandlerRot_Mat4(double quat[4], int ielem,
                         const CVec2& sp0, const CVec2& sp1, double mat[16],
                         const float mMV[16], const float mPj[16]);
CVec3 drag_AxisHandler(const CVec2& sp0,
                          const CVec2& sp1,
                          const CVec3& p,
                          const CVec3& axis,
                          double len,
                          const float* mMV,
                          const float* mPj);
bool isPickPoint(const CVec2& sp,
                 const CVec3& p,
                 const float* mMV,
                 const float* mPj,
                 double pick_tol);
bool isPick_AxisHandler(const CVec2& sp,
                        const CVec3& p,
                        const CVec3& axis,
                        double len,
                        const float* mMV,
                        const float* mPj,
                        double pick_tol);
bool isPickQuad(const CVec3& p0,
                const CVec3& p1,
                const CVec3& p2,
                const CVec3& p3,
                const CVec2& sp,
                const CVec3& pick_dir,
                const float mMV[16],
                const float mPj[16],
                double eps);

bool isPickCircle(const CVec3& axis,
                  const CVec3& org,
                  double rad,
                  const CVec3& src,
                  const CVec3& dir,
                  double pick_tol);
double DragCircle(const CVec2& sp0,
                  const CVec2& sp1,
                  const CVec3& p,
                  const CVec3& axis,
                  const float* mMV,
                  const float* mPj);
  
bool isPickCircle(const CVec2& sp,
                  const CVec3& p,
                  const CVec3& axis,
                  double r,
                  const float* mMV,
                  const float* mPj,
                  double pick_tol);
 
void Energy_MIPS(double& E, double dE[3][3], double ddE[3][3][3][3],
                 const double c[3][3],
                 const double C[3][3]);

  
}

#endif /* vec23mat3quat_hpp */
