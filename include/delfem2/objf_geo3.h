/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_OBJFUNC_v23_H
#define DFM2_OBJFUNC_v23_H

#include "delfem2/dfm2_inline.h"
#include <vector>
#include "delfem2/mat3.h"
#include "delfem2/vec3.h"

namespace delfem2 {

DFM2_INLINE void PBD_Pre3D(
    std::vector<double>& aXYZt,
    double dt,
    const double gravity[3],
    const std::vector<double>& aXYZ,
    const std::vector<double>& aUVW,
    const std::vector<int>& aBCFlag);

DFM2_INLINE void PBD_Post(
    std::vector<double>& aXYZ,
    std::vector<double>& aUVW,
    double dt,
    const std::vector<double>& aXYZt,
    const std::vector<int>& aBCFlag);

DFM2_INLINE void PBD_Update_Const3_Point3_Dim3(
    std::vector<double>& aXYZt,
    const double m[3],
    const double C[3],
    const double dCdp[3][9],
    const int aIP[3]);

DFM2_INLINE void PBD_Update_Const3(
    double* aXYZt,
    const int np,
    const int ndim,
    const double* m,
    const double* C,
    const double* dCdp,
    const unsigned int* aIP,
    double ratio);

DFM2_INLINE void PBD_ConstProj_Rigid2D(
    double* aXYt,
    double stiffness,
    const unsigned int *clstr_ind, unsigned int nclstr_ind,
    const unsigned int *clstr, unsigned int nclstr0,
    const double* aXY0, unsigned int nXY0);

DFM2_INLINE void PBD_ConstProj_Rigid3D(
    double* aXYZt,
    double stiffness,
    const int* clstr_ind, int nclstr_ind,
    const int* clstr,     int nclstr0,
    const double* aXYZ0,   int nXYZ0);

DFM2_INLINE void PBD_CdC_TriStrain2D3D(
    double C[3],
    double dCdp[3][9],
    const double P[3][2], // (in) undeformed triangle vertex positions
    const double p[3][3]); // (in) deformed triangle vertex positions

DFM2_INLINE void PBD_ConstraintProjection_DistanceTri2D3D(
    double C[3],
    double dCdp[3][9],
    const double P[3][2], // (in) undeformed triangle vertex positions
    const double p[3][3]); // (in) deformed triangle vertex positions

DFM2_INLINE void PBD_ConstraintProjection_EnergyStVK(
    double& C,
    double dCdp[9],
    const double P[3][2], // (in) undeformed triangle vertex positions
    const double p[3][3], // (in) deformed triangle vertex positions)
    const double lambda,
    const double myu);

DFM2_INLINE void PBD_ConstraintProjection_DistanceTet(
    double C[6],
    double dCdp[6][12],
    const double P[4][3], // (in) undeformed triangle vertex positions
    const double p[4][3]); // (in) deformed triangle vertex positions

DFM2_INLINE void PBD_CdC_QuadBend(
    double C[3],
    double dCdp[3][12],
    const double P[4][3],
    const double p[4][3]);

DFM2_INLINE void PBD_Seam(
    double* aXYZt,
    unsigned int nXYZ,
    const unsigned int* aLine,
    unsigned int nline);

DFM2_INLINE void WdWddW_MIPS(
    double& E, double dE[3][3], double ddE[3][3][3][3],
    const double c[3][3],
    const double C[3][3]);
  
/**
 * @brief energy W and its derivative dW and second derivative ddW
 * where W = a^T R(dn) b(theta)
 */
DFM2_INLINE void RodFrameTrans
 (CVec3d frm[3],
  const CVec3d& S0,
  const CVec3d& V01,
  const CVec3d& du,
  double dtheta);


/**
 * @param dF_dv (out) how i-th frame axis moves w.r.t the movement of the edge
 * @param dF_dt (out) how i-th frame axis moves w.r.t the rotation along 2nd axis
 * @param l01 (in) length of the edge
 * @param Frm (in) frame axis vectors
 */
DFM2_INLINE void DiffFrameRod
 (CMat3d dF_dv[3],
  CVec3d dF_dt[3],
  //
  double l01,
  const CVec3d Frm[3]);

/**
 * @brief second derivative (ddW) of W, where W = Q^T Frm[iaxis]
 */
DFM2_INLINE void DifDifFrameRod
 (CMat3d& ddW_ddv,
  CVec3d& ddW_dvdt,
  double& ddW_dtt,
  //
  unsigned int iaxis,
  double l01,
  const CVec3d& Q,
  const CVec3d Frm[3]);

DFM2_INLINE double WdWddW_DotFrame
 (CVec3d dV_dP[3],
  double dV_dt[2],
  CMat3d ddV_ddP[3][3],
  CVec3d ddV_dtdP[2][3],
  double ddV_ddt[2][2],
  //
  const CVec3d P[3],
  const CVec3d S[2],
  const double off[3]);

DFM2_INLINE void Darboux_Rod
 (CVec3d& darboux,
  //
  const CVec3d P[3],
  const CVec3d S[2]);

/**
 * @param P[3] (in) point of rod
 * @param S[2] (in) director vectors on the edges
 * @param off[3] (in) Daboux vector in the material frame
 * @param is_eaxct
 */
DFM2_INLINE double WdWddW_Rod
 (CVec3d dW_dP[3],
  double dW_dt[2],
  CMat3d ddW_ddP[3][3],
  CVec3d ddW_dtdP[2][3],
  double ddW_ddt[2][2],
  //
  const CVec3d P[3],
  const CVec3d S[2],
  const CVec3d& darboux0,
  bool is_exact);

DFM2_INLINE double WdWddW_SquareLengthLineseg3D
 (CVec3d dW_dP[2],
  CMat3d ddW_ddP[2][2],
  //
  const CVec3d P[2],
  double L0);

DFM2_INLINE double W_ArapEnergy
 (const std::vector<double>& aXYZ0,
  const std::vector<double>& aXYZ1,
  const std::vector<double>& aQuat1,
  const std::vector<unsigned int>& psup_ind,
  const std::vector<unsigned int>& psup);

DFM2_INLINE void dW_ArapEnergy
 (std::vector<double>& aRes,
  const std::vector<double>& aXYZ0,
  const std::vector<double>& aXYZ1,
  const std::vector<double>& aQuat1,
  const std::vector<unsigned int>& psup_ind,
  const std::vector<unsigned int>& psup);

DFM2_INLINE void ddW_ArapEnergy
 (std::vector<double>& eM,
  const std::vector<unsigned int>& aIP,
  const std::vector<double>& aXYZ0,
  const std::vector<double>& aQuat1);


template <typename T>
DFM2_INLINE void GetConstConstDiff_Bend
 (double& C,
  CVec3<T> dC[4],
  // -----
  const CVec3<T>& p0,
  const CVec3<T>& p1,
  const CVec3<T>& p2,
  const CVec3<T>& p3);

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/objf_geo3.cpp"
#endif

#endif /* pbd_v23_h */
