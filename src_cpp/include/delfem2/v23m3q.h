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
CMatrix3 OuterProduct(const CVector3& vec0,
                      const CVector3& vec1 );
CMatrix3 Spin(const CVector3& vec0);
CMatrix3 ParallelTransport(const CVector3& p0,
                           const CVector3& p1,
                           const CVector3& q0,
                           const CVector3& q1);
CMatrix3 MinimumRotation(const CVector3& V,
                         const CVector3& v);
CMatrix3 Irot_Tri(const CVector3& d0,
                  const CVector3& d1,
                  const CVector3& d2);
CMatrix3 Irot_TriSolid(const CVector3& d0,
                       const CVector3& d1,
                       const CVector3& d2);
CMatrix3 Irot_LineSeg(const CVector3& d0,
                      const CVector3& d1);
CMatrix3 Irot_Point(const CVector3& d0);

#endif /* vec23mat3quat_hpp */
