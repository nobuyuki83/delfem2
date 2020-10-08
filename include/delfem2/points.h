/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file interfaces for functions that handle coordinates of points stored in a flat vector.
 * @brief interfaces for functions that handle coordinates of points stored in a flat vector.
 * @details this file should not depends on anything except for  std::vector
 */

#ifndef DFM2_POINTS_H
#define DFM2_POINTS_H
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
DFM2_INLINE void BoundingBox3_Points3(
    T min3[3],
    T max3[3],
    const T* aXYZ,
    unsigned int nXYZ);

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
  unsigned int nnode_,
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
  unsigned int nnode_,
  T s);


/**
 * @brief uniformly scale the coordinate of opints such that the longest edge of axis-aligned boudning box is the scale.
 * @param length_longest_aabb_edge length of longest edge of axis-aligned bounding box
 */
template <typename REAL>
DFM2_INLINE void Normalize_Points3(
    std::vector<REAL>& aXYZ,
    REAL length_longest_aabb_edge = 1.0);


/**
 * @brief scale each vector to make norm == 1
 * @param aVec coordinates packed in a single std::vector array
 * @param ndim dimension (must be either 2 or 3)
 */
template <typename REAL>
DFM2_INLINE void NormalizeVector_Points(
    REAL* aVec,
    unsigned int np,
    unsigned int ndim);


DFM2_INLINE double Size_Points3D_LongestAABBEdge
 (const std::vector<double>& aXYZ);
  
/**
 * @details implemented for "float" and "double"
 */
template <typename T>
DFM2_INLINE void CG_Point3
 (T cg[3],
  const std::vector<T>& aXYZ);

template <typename T>
DFM2_INLINE void Points_RandomUniform(
    T* aOdir,
    unsigned int npo,
    unsigned int ndim,
    const T* minCoords,
    const T* maxCoords);

DFM2_INLINE void TangentVector_Points3(
    std::vector<double>& aOdir,
    const std::vector<double>& aNorm);



} // delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/points.cpp"
#endif


#endif
