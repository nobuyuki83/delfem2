/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file should not depends on anything except for  std::vector
 */

#ifndef DFM2_MSHMISC_H
#define DFM2_MSHMISC_H

#include <vector>

#include "delfem2/dfm2_inline.h"

// -----------------
// work on points

namespace delfem2 {

DFM2_INLINE void GetCenterWidthGroup(
    double &cx, double &cy, double &cz,
    double &wx, double &wy, double &wz,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElem,
    const int nnoel,
    int igroup,
    const std::vector<int> &aIndGroup);

DFM2_INLINE void GetCenterWidthGroup(
    double &cx, double &cy, double &cz,
    double &wx, double &wy, double &wz,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    int igroup,
    const std::vector<int> &aIndGroup);

DFM2_INLINE void GetCenterWidth3DGroup(
    double cw[6],
    //
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    int igroup,
    const std::vector<int> &aIndGroup);

DFM2_INLINE void CG_Tri(
    double &cgx,
    double &cgy,
    double &cgz,
    int itri,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri);

/**
 * @brief center positions of each triangle and the maximum radius of the triangle
 * @details this funciton is implemented for "float" and double.
 * the aXYZ_c0 will be resized to aTri.size()/3
 */
template<typename T>
DFM2_INLINE T CentsMaxRad_MeshTri3(
    std::vector<T> &aXYZ_c0,
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTri);

template<typename T>
DFM2_INLINE void CG_MeshTri3_Shell(
    T cg[3],
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTri);

template<typename T>
DFM2_INLINE T CG_TriMsh3Flg_Shell(
    T cg[3],
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTri,
    int iflg,
    const std::vector<int> &aFlg);

template<typename T>
DFM2_INLINE void CG_MeshTri3_Solid(
    T cg[3],
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTri);

template<typename T>
DFM2_INLINE void CG_MeshTet3(
    T &v_tot,
    T cg[3],
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTet);

// ---------------------------------

DFM2_INLINE void RemoveUnreferencedPoints_MeshElem(
    std::vector<double> &aXYZ1,
    std::vector<unsigned int> &aElem1,
    std::vector<int> &aMap01,
    unsigned int ndim,
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aElem0);

/**
 * @brief Normal at the vertex of a triangle mesh.
 */
template<typename REAL>
DFM2_INLINE void Normal_MeshTri3D(
    REAL *aNorm,
    const REAL *aXYZ,
    size_t nXYZ,
    const unsigned int *aTri,
    size_t nTri);

/**
 * @brief Normal at the vertex of a quad mesh. Defined for "float" and "double"
 */
template<typename REAL>
DFM2_INLINE void Normal_MeshQuad3(
    std::vector<REAL> &aNorm,
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aQuad);

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
    std::vector<double> &aXYZ,
    const std::vector<int> &aTri,
    const std::vector<int> &elsup_ind,
    const std::vector<int> elsup);

// ---------------------------------------------------------

DFM2_INLINE void AddMesh(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri0);

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mshmisc.cpp"
#endif

#endif /* DFM2_MSHMISC_H */
