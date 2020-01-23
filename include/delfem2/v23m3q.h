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

CVec3 operator* (const CVec3& v, const CMat3d& m);
CVec3 operator* (const CMat3d& m, const CVec3& v);
  
CVec3 MatVec(const CMat3d& m, const CVec3& vec0);
CVec3 MatVecTrans(const CMat3d& m, const CVec3& vec0);

// ---------------------------------------------

CVec2d screenXYProjection(const CVec3& v,
                          const float* mMV,
                          const float* mPj);

void SetProjection(CMat3d& m, const CVec3& vec0);
void SetDiag(CMat3d& m, const CVec3& d);
void SetRotMatrix_Cartesian(CMat3d& m, const CVec3& v);
void SetSpinTensor(CMat3d& m, const CVec3& vec0);
void SetOuterProduct(CMat3d& m,
                     const CVec3& vec0,
                     const CVec3& vec1 );
CVec3 GetSpinVector(const CMat3d& m);
CVec3 GetCartesianRotationVector(const CMat3d& m);

  
// --------------------------------------------

CMat3d Mat3(const CVec3& vec0,
              const CVec3& vec1,
              const CVec3& vec2);
CMat3d Mat3(const CVec3& vec0);
CMat3d Mat3(const CVec3& vec0,
              const CVec3& vec1);
CMat3d RotMatrix_Cartesian(const CVec3& v);
CMat3d Mat3_RotCartesian(const CVec3& vec0);
CMat3d Mat3_OuterProduct(const CVec3& vec0,
                      const CVec3& vec1 );
CMat3d Mat3_Spin(const CVec3& vec0);
CMat3d Mat3_ParallelTransport(const CVec3& p0,
                                const CVec3& p1,
                                const CVec3& q0,
                                const CVec3& q1);
CMat3d Mat3_MinimumRotation(const CVec3& V,
                              const CVec3& v);
CMat3d Mat3_IrotTri(const CVec3& d0,
                      const CVec3& d1,
                      const CVec3& d2);
CMat3d Mat3_IrotTriSolid(const CVec3& d0,
                           const CVec3& d1,
                           const CVec3& d2);
CMat3d Mat3_IrotLineSeg(const CVec3& d0,
                          const CVec3& d1);
CMat3d Mat3_IrotPoint(const CVec3& d0);
CMat3d Mirror(const CVec3& n);

void Mat4_MatTransl(double m[16],
                    const CMat3d& mat, const CVec3& trans);
void Mat4_ScaleMatTransl(double m[16],
                         double scale, const CMat3d& mat, const CVec3& trans);

int PickHandlerRotation_PosQuat(const CVec3& src, const CVec3& dir,
                                const CVec3& pos, const double quat[4], double rad,
                                double tol);
int PickHandlerRotation_Mat4(const CVec3& src, const CVec3& dir,
                             const double mat[16], double rad,
                             double tol);
bool DragHandlerRot_PosQuat(double quat[4], int ielem,
                            const CVec2d& sp0, const CVec2d& sp1,
                            const CVec3& pos,
                            const float mMV[16], const float mPj[16]);
bool DragHandlerRot_Mat4(double quat[4], int ielem,
                         const CVec2d& sp0, const CVec2d& sp1, double mat[16],
                         const float mMV[16], const float mPj[16]);
CVec3 drag_AxisHandler(const CVec2d& sp0,
                       const CVec2d& sp1,
                       const CVec3& p,
                       const CVec3& axis,
                       double len,
                       const float* mMV,
                       const float* mPj);
bool isPickPoint(const CVec2d& sp,
                 const CVec3& p,
                 const float* mMV,
                 const float* mPj,
                 double pick_tol);
bool isPick_AxisHandler(const CVec2d& sp,
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
                const CVec2d& sp,
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
double DragCircle(const CVec2d& sp0,
                  const CVec2d& sp1,
                  const CVec3& p,
                  const CVec3& axis,
                  const float* mMV,
                  const float* mPj);
  
bool isPickCircle(const CVec2d& sp,
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
