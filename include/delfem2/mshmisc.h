/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file should not depends on anything except for  std::vector, std::functional
 */

#ifndef DFM2_MSHMISC_H
#define DFM2_MSHMISC_H

#include <vector>
#include <functional>  // maybe we should separate the functions with this dependency

#include "delfem2/dfm2_inline.h"

// -----------------
// work on points

namespace delfem2 {

DFM2_INLINE void RemoveUnreferencedPoints_MeshElem(
    std::vector<double> &aXYZ1,
    std::vector<unsigned int> &aElem1,
    std::vector<int> &aMap01,
    unsigned int ndim,
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aElem0);

DFM2_INLINE void Quality_MeshTri2D(
    double &max_aspect, 
	double &min_area,
    const double *aXY,
    const unsigned int *aTri,
    unsigned int nTri);

// ----------------------
// set primitive mesh

DFM2_INLINE void SetTopology_ExtrudeTri2Tet(
    unsigned int *aTet,
    int nXY,
    const unsigned int *aTri, 
	int nTri,
    int nlayer);

DFM2_INLINE void ExtrudeTri2Tet(
    int nlayer, double h,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    const std::vector<double> &aXY,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE double SolidAngleTri3D(
    const double v1[3],
    const double v2[3],
    const double v3[3]);

DFM2_INLINE void makeSolidAngle(
    std::vector<double> &aSolidAngle,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aNorm,
    std::vector<int> &elsup_ind,
    std::vector<int> &elsup);

/**
 * @brief Compute mass of the points (lumped mass) for 3D tet mesh
 * @details aMassMatrixLumped need to be allocated before in the size of nXY
 * this is here because both "fem" and "pbd" use this function
 */
DFM2_INLINE void MassPoint_Tet3D(
    double *aMassMatrixLumped,
    double rho,
    const double *aXYZ, 
	size_t nXYZ,
    const unsigned int *aTet, 
	size_t nTet);

/**
 * @brief Compute mass of the points (lumped mass) for 2D triangle mesh
 * @param aMass (out) this need to be allocated before in the size of nXY
 * @details this is here because both "fem" and "pbd" use this function
 */
DFM2_INLINE void MassPoint_Tri2D(
    double *aMassMatrixLumped,
    //
    double rho,
    const double *aXY,
    size_t nXY,
    const unsigned int *aTri,
    size_t nTri);

/**
 * @brief Compute mass of the points (lumped mass) for 3D triangle mesh
 * @param aMass (out) this need to be allocated before in the size of nXY
 * @details this is here because both "fem" and "pbd" use this function
 */
DFM2_INLINE void MassPoint_Tri3D(
    double *aMass,
    double rho,
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aTri,
    size_t nTri);

DFM2_INLINE void LaplacianSmoothing(
    std::vector<double> &vtx_xyz,
    const std::vector<int> &tri_vtx,
    const std::vector<int> &elsup_ind,
    const std::vector<int> &elsup);

// ---------------------------------------------------------

DFM2_INLINE void AddMesh(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri0);

DFM2_INLINE double Area_MeshTri3(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx,
    const std::function<bool(unsigned int)>& flag);


DFM2_INLINE bool Distortion_MappingTriangleFrom2To3Dim(
    double thresA,
    double thresE,
    unsigned int it0,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ,
    const std::vector<double> &aTexP);

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mshmisc.cpp"
#endif

#endif /* DFM2_MSHMISC_H */
