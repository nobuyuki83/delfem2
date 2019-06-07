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

void Energy_MIPS(double& E, double dE[3][3], double ddE[3][3][3][3],
                 const double c[3][3],
                 const double C[3][3]);

void ConstraintProjection_Rigid2D(double* aXYt,
                                  double stiffness,
                                  const int* clstr_ind, int nclstr_ind,
                                  const int* clstr,     int nclstr0,
                                  const double* aXY0,   int nXY0);

void ConstraintProjection_Rigid3D(double* aXYZt,
                                  double stiffness,
                                  const int* clstr_ind, int nclstr_ind,
                                  const int* clstr,     int nclstr0,
                                  const double* aXYZ0,   int nXYZ0);

void PBD_ConstraintProjection_Strain(double C[3],
                              double dCdp[3][9],
                              const double P[3][2], // (in) undeformed triangle vertex positions
                              const double p[3][3]); // (in) deformed triangle vertex positions
void PBD_ConstraintProjection_DistanceTri2D3D(double C[3],
                                              double dCdp[3][9],
                                              const double P[3][2], // (in) undeformed triangle vertex positions
                                              const double p[3][3]); // (in) deformed triangle vertex positions
void PBD_ConstraintProjection_EnergyStVK(double& C,
                                         double dCdp[9],
                                         const double P[3][2], // (in) undeformed triangle vertex positions
                                         const double p[3][3], // (in) deformed triangle vertex positions)
                                         const double lambda,
                                         const double myu);
void PBD_ConstraintProjection_DistanceTet(double C[6],
                                          double dCdp[6][12],
                                          const double P[4][3], // (in) undeformed triangle vertex positions
                                          const double p[4][3]); // (in) deformed triangle vertex positions

void Check_ConstraintProjection_DistanceTri2D3D(const double P[3][2], // (in) undeformed triangle vertex positions
                                                const double p[3][3]); // (in) deformed triangle vertex positions)
void Check_ConstraintProjection_Strain(const double P[3][2], // (in) undeformed triangle vertex positions
                                    const double p[3][3]); // (in) deformed triangle vertex positions)
void Check_ConstraintProjection_EnergyStVK(const double P[3][2], // (in) undeformed triangle vertex positions
                                           const double p[3][3], // (in) deformed triangle vertex positions)
                                           const double lambda,
                                           const double myu);

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
