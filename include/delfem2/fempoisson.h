/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file make element matrix (EMAT) and merge it to the global matri for linear solid equation
 *
 * (2021/05/05) changed header guard name
 * (2020/12/26) TODO: use template to generalize the merge functions
 * (2020/12/25) created. separated from "femem3" and "femem2"
 */

#ifndef DFM2_FEMPOISSON_H
#define DFM2_FEMPOISSON_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/femutil.h"

#ifdef DFM2_STATIC_LIBRARY
// Merge use explicitly use the template so for static library we need to include the template itself.
#  include "delfem2/lsmats.h"
#endif

namespace delfem2 {

// ------------------------------------------------------------------
// below: fem element matrix for 2D mesh

// -[\alpha]\nabla^2[value] = [source]
void EMat_Poisson_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha,
    const double source,
    const double coords[3][2],
    const double value[3]);

void EMat_Poission2_QuadOrth(
    double emat[4][4],
    double lx,
    double ly);

/**
 *
 * @param emat
 * @param lx
 * @param ly
 * @param[in] ngauss ngauss=1 is enough for analytically exact integration
 */
void EMat_Poisson2_QuadOrth_GaussInt(
    double emat[4][4],
    double lx,
    double ly,
    unsigned int ngauss);

DFM2_INLINE void EMat_Poisson_Tet3D(
    double eres[4],
    double emat[4][4],
    const double alpha, const double source,
    const double coords[4][3],
    const double value[4]);

// [\rho][velo] - [\alpha]\nabla^2[value] = [source]
void EMat_Diffusion_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha,
    const double source,
    const double dt,
    const double gamma,
    const double rho,
    const double coords[3][2],
    const double value[3],
    const double velo[3]);

DFM2_INLINE void EMat_Diffusion_Newmark_Tet3D(
    double eres[4],
    double emat[4][4],
    const double alpha, const double source,
    const double dt_timestep, const double gamma_newmark, const double rho,
    const double coords[4][3],
    const double value[4], const double velo[4]);

// --------------------------------------------


template <class MAT>
void MergeLinSys_Poission_MeshTri2D(
    MAT& mat_A,
    double* vec_b,
    const double alpha,
    const double source,
    const double* aXY1,
    size_t np,
    const unsigned int* aTri1,
    size_t nTri,
    const double* aVal)
{
  const size_t nDoF = np;
  //
  std::vector<unsigned int> tmp_buffer(nDoF, UINT_MAX);
  for (unsigned int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData<3,2>(coords, aIP,aXY1);
    const double value[3] = { aVal[i0], aVal[i1], aVal[i2] };
    //
    double eres[3], emat[3][3];
    EMat_Poisson_Tri2D(
        eres,emat,
        alpha, source,
        coords, value);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    Merge<3,3,double>(mat_A,aIP,aIP,emat,tmp_buffer);
  }
}

template <class MAT>
void MergeLinSys_Poission_MeshTet3D(
    MAT& mat_A,
    double* vec_b,
    const double alpha,
    const double source,
    const double* aXYZ,
    size_t nXYZ,
    const unsigned int* aTet,
    size_t nTet,
    const double* aVal)
{
  const size_t np = nXYZ;
  std::vector<unsigned int> tmp_buffer(np, UINT_MAX);
  for (unsigned int itet = 0; itet<nTet; ++itet){
    const unsigned int i0 = aTet[itet*4+0];
    const unsigned int i1 = aTet[itet*4+1];
    const unsigned int i2 = aTet[itet*4+2];
    const unsigned int i3 = aTet[itet*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData<4,3>(coords, aIP,aXYZ);
    const double value[4] = { aVal[i0], aVal[i1], aVal[i2], aVal[i3] };
    //
    double eres[4], emat[4][4];
    EMat_Poisson_Tet3D(
        eres,emat,
        alpha, source,
        coords, value);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    Merge<4,4,double>(mat_A,aIP,aIP,emat,tmp_buffer);
  }
}

template <class MAT>
void MergeLinSys_Diffusion_MeshTri2D(
    MAT& mat_A,
    double* vec_b,
    const double alpha,
    const double rho,
    const double source,
    const double dt_timestep,
    const double gamma_newmark,
    const double* aXY1,
    size_t nXY,
    const unsigned int* aTri1,
    size_t nTri,
    const double* aVal,
    const double* aVelo)
{
  std::vector<unsigned int> tmp_buffer(nXY, UINT_MAX);
  for (unsigned int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData<3,2>(coords, aIP,aXY1);
    const double value[3] = { aVal[ i0], aVal[ i1], aVal[ i2] };
    const double velo[ 3] = { aVelo[i0], aVelo[i1], aVelo[i2] };
    // --
    double eres[3], emat[3][3];
    EMat_Diffusion_Tri2D(
        eres,emat,
        alpha, source,
        dt_timestep, gamma_newmark, rho,
        coords, value, velo);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    Merge<3,3,double>(mat_A,aIP,aIP,emat,tmp_buffer);
  }
}

template <class MAT>
void MergeLinSys_Diffusion_MeshTet3D(
    MAT& mat_A,
    double* vec_b,
    const double alpha,
    const double rho,
    const double source,
    const double dt_timestep,
    const double gamma_newmark,
    const double* aXYZ,
    size_t nXYZ,
    const unsigned int* aTet,
    size_t nTet,
    const double* aVal,
    const double* aVelo)
{
  const size_t np = nXYZ;
  std::vector<unsigned int> tmp_buffer(np, UINT_MAX);
  for (unsigned int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData<4,3>(coords, aIP,aXYZ);
    const double value[4] = { aVal[ i0], aVal[ i1], aVal[ i2], aVal[ i3] };
    const double velo[ 4] = { aVelo[i0], aVelo[i1], aVelo[i2], aVelo[i3] };
    // ---------------------
    double eres[4], emat[4][4];
    EMat_Diffusion_Newmark_Tet3D(
        eres,emat,
        alpha, source,
        dt_timestep, gamma_newmark, rho,
        coords, value, velo);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
//    mat_A.Mearge(4, aIP, 4, aIP, 1, &emat[0][0], tmp_buffer);
    Merge<4,4,double>(mat_A,aIP,aIP,emat,tmp_buffer);
  }
}


} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fempoisson.cpp"
#endif
  
#endif /* fem_ematrix_h */
