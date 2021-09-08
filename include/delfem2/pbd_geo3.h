/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_PBD_GEO3_H
#define DFM2_PBD_GEO3_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/mat3.h"
#include "delfem2/vec3.h"

namespace delfem2 {

DFM2_INLINE void PBD_Pre3D(
    std::vector<double> &aXYZt,
    double dt,
    const double gravity[3],
    const std::vector<double> &aXYZ,
    const std::vector<double> &aUVW,
    const std::vector<int> &aBCFlag);

DFM2_INLINE void PBD_Post(
    std::vector<double> &aXYZ,
    std::vector<double> &aUVW,
    double dt,
    const std::vector<double> &aXYZt,
    const std::vector<int> &aBCFlag);

DFM2_INLINE void PBD_Update_Const3_Point3_Dim3(
    std::vector<double> &aXYZt,
    const double m[3],
    const double C[3],
    const double dCdp[3][9],
    const int aIP[3]);

DFM2_INLINE void PBD_Update_Const3(
    double *aXYZt,
    const int np,
    const int ndim,
    const double *m,
    const double *C,
    const double *dCdp,
    const unsigned int *aIP,
    double ratio);

DFM2_INLINE void PBD_ConstProj_Rigid2D(
    double *aXYt,
    double stiffness,
    const unsigned int *clstr_ind, 
	size_t nclstr_ind,
    const unsigned int *clstr, 
	size_t nclstr0,
    const double *aXY0, 
	size_t nXY0);

DFM2_INLINE void PBD_ConstProj_Rigid3D(
    double *aXYZt,
    double stiffness,
    const int *clstr_ind, int nclstr_ind,
    const int *clstr, int nclstr0,
    const double *aXYZ0, int nXYZ0);

/**
 *
 * @param C
 * @param dCdp
 * @param[in] P undeformed triangle vertex positions
 * @param[in] p deformed triangle vertex positions
 */
DFM2_INLINE void PBD_CdC_TriStrain2D3D(
    double C[3],
    double dCdp[3][9],
    const double P[3][2],
    const double p[3][3]);

/**
 *
 * @param C
 * @param dCdp
 * @param[in] P undeformed triangle vertex positions
 * @param[in] p deformed triangle vertex positions
 */
DFM2_INLINE void PBD_ConstraintProjection_DistanceTri2D3D(
    double C[3],
    double dCdp[3][9],
    const double P[3][2],
    const double p[3][3]);

/**
 *
 * @param C
 * @param dCdp
 * @param[in] P undeformed triangle vertex positions
 * @param[in] p deformed triangle vertex positions
 * @param lambda
 * @param myu
 */
DFM2_INLINE void PBD_ConstraintProjection_EnergyStVK(
    double &C,
    double dCdp[9],
    const double P[3][2],
    const double p[3][3],
    const double lambda,
    const double myu);

/**
 *
 * @param C
 * @param dCdp
 * @param[in] P undeformed triangle vertex positions
 * @param[in] p deformed triangle vertex positions
 */
DFM2_INLINE void PBD_ConstraintProjection_DistanceTet(
    double C[6],
    double dCdp[6][12],
    const double P[4][3],
    const double p[4][3]);

DFM2_INLINE void PBD_CdC_QuadBend(
    double C[3],
    double dCdp[3][12],
    const double P[4][3],
    const double p[4][3]);

DFM2_INLINE void PBD_Seam(
    double *aXYZt,
    size_t nXYZ,
    const unsigned int *aLine,
    size_t nline);

template<typename T>
DFM2_INLINE void GetConstConstDiff_Bend(
    double &C,
    CVec3<T> dC[4],
    // -----
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3);
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/pbd_geo3.cpp"
#endif

#endif /* pbd_v23_h */
