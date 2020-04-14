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
#include "delfem2/dfm2_inline.h"
#include <vector>

// -----------------
// work on points

namespace delfem2{

/**
 * @brief update minimum and maximum coordinates
 * @details implemented for "float" and "double"
 */
template<typename T>
DFM2_INLINE void updateMinMaxXYZ
 (T& x_min, T& x_max,
  T& y_min, T& y_max,
  T& z_min, T& z_max,
  T x, T y, T z);

/**
 * @param bb3 (out) bounding box in the order of <minx, miny, minz, maxx, maxy, maxz>
 * @param aXYZ (in) array of 3D coordinates of points
 * @param nXYZ (in) number of points
 * @details implemented for "float" and "double"
 */
template<typename T>
DFM2_INLINE void Min3Max3_Points3(
    T min3[3],
    T max3[3],
    const T* aXYZ,
    const unsigned int nXYZ);

// center & width
template <typename T>
DFM2_INLINE void CenterWidth_Point3
 (T& cx, T& cy, T& cz,
  T& wx, T& wy, T& wz,
  const T* paXYZ, unsigned int nXYZ);
  
template <typename T>
DFM2_INLINE void CenterWidth_Points3
 (T& cx, T& cy, T& cz,
  T& wx, T& wy, T& wz,
  const std::vector<T>& aXYZ);
  
template <typename T>
DFM2_INLINE void CenterWidth_Points3
 (T c[3],
  T w[3],
  const std::vector<T>& aXYZ);
  
DFM2_INLINE void GetCenterWidthGroup
 (double& cx, double& cy, double& cz,
  double& wx, double& wy, double& wz,
  const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aElem,
  const int nnoel,
  int igroup,
  const std::vector<int>& aIndGroup);

DFM2_INLINE void GetCenterWidthGroup
 (double& cx, double& cy, double& cz,
  double& wx, double& wy, double& wz,
  const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  int igroup,
  const std::vector<int>& aIndGroup);

DFM2_INLINE void GetCenterWidth3DGroup
 (double cw[6],
                     //
  const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  int igroup,
  const std::vector<int>& aIndGroup);

// local coordinate
DFM2_INLINE void GetCenterWidthLocal
 (double& lcx, double& lcy, double& lcz,
  double& lwx, double& lwy, double& lwz,
  const std::vector<double>& aXYZ,
  const double lex[3],
  const double ley[3],
  const double lez[3]);


// ------------------------------------------


/**
 * @brief rotate with the Bryant angle (in the  order of XYZ) around the origin.
 * @details the angles are in the radian.
 */
template <typename T>
DFM2_INLINE void Rotate_Points3
 (std::vector<T>& aXYZ,
  T radx, T rady, T radz);

template <typename T>
DFM2_INLINE void Translate_Points3
 (std::vector<T>& aXYZ,
  T tx, T ty, T tz);

template <typename T>
DFM2_INLINE void Translate_Points3
 (T* pXYZs_,
  const unsigned int nnode_,
  T tx, T ty, T tz);

template <typename T>
DFM2_INLINE void Translate_Points2
 (std::vector<T>& aXY,
  T tx, T ty);
  
template <typename T>
DFM2_INLINE void Scale_PointsX
 (std::vector<T>& aXYZ,
                   T s);

template <typename T>
DFM2_INLINE void Scale_Points3
 (T* pXYZs_,
  const unsigned int nnode_,
  T s);


/**
 * @brief uniformly scale the coordinate of opints such that the longest edge of axis-aligned boudning box is the scale.
 * @param length_longest_aabb_edge length of longest edge of axis-aligned bounding box
 */
template <typename REAL>
DFM2_INLINE void Normalize_Points3
 (std::vector<REAL>& aXYZ,
  REAL length_longest_aabb_edge = 1.0);

DFM2_INLINE double Size_Points3D_LongestAABBEdge
 (const std::vector<double>& aXYZ);
  
/**
 * @details implemented for "float" and "double"
 */
template <typename T>
DFM2_INLINE void CG_Point3
 (T cg[3],
  const std::vector<T>& aXYZ);

// points above here
// ----------------------------------------------------------------------------------------------
// mesh from here
  

DFM2_INLINE void CG_Tri
 (double& cgx, double& cgy, double& cgz,
  int itri,
  const std::vector<double>& aXYZ,
  const std::vector<int>& aTri);


/**
 * @brief center positions of each triangle and the maximum radius of the triangle
 * @details this funciton is implemented for "float" and double.
 * the aXYZ_c0 will be resized to aTri.size()/3
 */
template <typename T>
DFM2_INLINE T CentsMaxRad_MeshTri3
 (std::vector<T>& aXYZ_c0,
  const std::vector<T>& aXYZ,
  const std::vector<unsigned int>& aTri);
  
template <typename T>
DFM2_INLINE void CG_MeshTri3_Shell
 (T cg[3],
  const std::vector<T>& aXYZ,
  const std::vector<unsigned int>& aTri);
  
template <typename T>
DFM2_INLINE T CG_TriMsh3Flg_Shell
 (T cg[3],
  const std::vector<T>& aXYZ,
  const std::vector<unsigned int>& aTri,
  int iflg,
  const std::vector<int>& aFlg);
  
template <typename T>
DFM2_INLINE void CG_MeshTri3_Solid
 (T cg[3],
  const std::vector<T>& aXYZ,
  const std::vector<unsigned int>& aTri);

template <typename T>
DFM2_INLINE void CG_MeshTet3
 (T& v_tot,
  T cg[3],
  const std::vector<T>& aXYZ,
  const std::vector<unsigned int>& aTet);
  
// ---------------------------------
  
DFM2_INLINE void RemoveUnreferencedPoints_MeshElem
 (std::vector<double>& aXYZ1,
  std::vector<unsigned int>& aElem1,
  std::vector<int>& aMap01,
  unsigned int ndim,
  const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aElem0);

/**
 * @brief Normal at the vertex of a triangle mesh.
 */
DFM2_INLINE void Normal_MeshTri3D
 (double *aNorm,
  const double *aXYZ, unsigned int nXYZ,
  const unsigned int *aTri, unsigned int nTri);
 
/**
 * @brief Normal at the vertex of a quad mesh. Defined for "float" and "double"
 */
template <typename REAL>
DFM2_INLINE void Normal_MeshQuad3(
    std::vector<REAL>& aNorm,
    const std::vector<REAL>& aXYZ,
    const std::vector<unsigned int>& aQuad);

DFM2_INLINE void Quality_MeshTri2D(
    double &max_aspect, double &min_area,
    const double *aXY,
    const unsigned int *aTri, unsigned int nTri);

// ----------------------
// set primitive mesh

DFM2_INLINE void SetTopology_ExtrudeTri2Tet(
    unsigned int *aTet,
    int nXY,
    const unsigned int *aTri, int nTri,
    int nlayer);

DFM2_INLINE void ExtrudeTri2Tet(
    int nlayer, double h,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    const std::vector<double> &aXY,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE double SolidAngleTri3D
 (const double v1[3],
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
    const double *aXYZ, unsigned int nXYZ,
    const unsigned int *aTet, unsigned int nTet);

/**
 * @brief Compute mass of the points (lumped mass) for 2D triangle mesh
 * @param aMassMatrixLumped (out) this need to be allocated before in the size of nXY
 * @details this is here because both "fem" and "pbd" use this function
 */
DFM2_INLINE void MassPoint_Tri2D(
    double *aMassMatrixLumped,
    //
    double rho,
    const double *aXY, unsigned int nXY,
    const unsigned int *aTri, unsigned int nTri);


DFM2_INLINE void LaplacianSmoothing
 (std::vector<double>& aXYZ,
  const std::vector<int>& aTri,
  const std::vector<int>& elsup_ind,
  const std::vector<int> elsup);

// ---------------------------------------------------------

DFM2_INLINE void SubdivisionPoints_QuadCatmullClark
 (std::vector<double>& aXYZ1,
  //
  const std::vector<unsigned int>& aQuad1,
  const std::vector<int>& aEdgeFace0,
  const std::vector<unsigned int> &psupIndQuad0,
  const std::vector<unsigned int> &psupQuad0,
  const unsigned int* aQuad0, unsigned int nQuad0,
  const double* aXYZ0, unsigned int nXYZ0);
  
DFM2_INLINE void SubdivPoints3_MeshQuad
 (std::vector<double>& aXYZ1,
                            //
  const std::vector<int>& aEdgeFace0,
  const std::vector<unsigned int>& aQuad0,
  const std::vector<double>& aXYZ0);

DFM2_INLINE void SubdivisionPoints_Hex
 (std::vector<double>& aXYZ1,
                          //
  const std::vector<unsigned int> &psupIndHex0,
  const std::vector<unsigned int> &psupHex0,
  const std::vector<unsigned int>& aQuadHex0,
  const unsigned int* aHex0, unsigned int nHex0,
  const double* aXYZ0, unsigned int nXYZ0);

} // delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/mshmisc.cpp"
#endif


#endif
