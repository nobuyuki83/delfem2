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

///////////////////////////////////////////////


///////////////////////////////////////////////

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

void AffineMatrixTrans(double m[16], const CMatrix3& mat, const CVector3& trans);

#endif /* vec23mat3quat_hpp */
