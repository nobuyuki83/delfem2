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

#include <vector>

#include "delfem2/dfm2_inline.h"

// -----------------
// work on points

namespace delfem2 {

// ------------------------------------------


/**
 * @brief rotate with the Bryant angle (in the  order of XYZ) around the origin.
 * @details the angles are in the radian.
 */
template<typename T>
DFM2_INLINE void Rotate_Points3(
    std::vector<T> &vec_xyz,
    T radx,
    T rady,
    T radz);

// ----

template<typename T>
DFM2_INLINE void Translate_Points(
    T *vtx_coords,
    const size_t num_vtx,
    const unsigned int ndim,
    const T *translation);

template<typename T>
DFM2_INLINE void Translate_Points2(
    std::vector<T> &aXY,
    T tx,
    T ty);

template<typename T>
DFM2_INLINE void Translate_Points3(
    std::vector<T> &aXYZ,
    T tx,
    T ty,
    T tz);

// ----

template<typename T>
DFM2_INLINE void Scale_PointsX(
    std::vector<T> &vtx_xyz,
    T scale);

template<typename T>
DFM2_INLINE void Scale_Points(
    T *vtx_coords,
    size_t num_vtx,
    unsigned int ndim,
    T scale);

/**
 * @brief uniformly scale & translate the coordinate of opints
 * specifying the longest edge of AABB and the center of AABB is origin
 * @param length_longest_edge_boundingbox length of longest edge of axis-aligned bounding box
 */
template<typename REAL>
DFM2_INLINE void Normalize_Points3(
    std::vector<REAL> &vtx_xyz,
    REAL length_longest_edge_boundingbox = 1);

/**
 * @brief scale each vector to make norm == 1
 * @param aVec coordinates packed in a single std::vector array
 * @param ndim dimension (must be either 2 or 3)
 */
template<typename REAL>
DFM2_INLINE void NormalizeVector_Points(
    REAL *aVec,
    unsigned int np,
    unsigned int ndim);

/**
 * @details implemented for "float" and "double"
 */
template<typename T>
DFM2_INLINE void CG_Point3(
    T *center_of_gravity,
    const std::vector<T> &vtx_xyz);

DFM2_INLINE double EnergyKinetic(
    const double *aUVW,
    size_t np);

template<typename T>
DFM2_INLINE void Points_RandomUniform(
    T *aOdir,
    size_t npo,
    unsigned int ndim,
    const T *minCoords,
    const T *maxCoords);

DFM2_INLINE void TangentVector_Points3(
    std::vector<double> &aOdir,
    const std::vector<double> &aNorm);

class CKineticDamper {
 public:
  void Damp(std::vector<double> &aUVW) {
    aEnergy.push_back(EnergyKinetic(aUVW.data(), aUVW.size() / 3));
    if (aEnergy.size() > 3) {
      aEnergy.erase(aEnergy.begin());
      const double g0 = aEnergy[1] - aEnergy[0];
      const double g1 = aEnergy[2] - aEnergy[1];
      if (g0 > 0 && g1 < 0) { aUVW.assign(aUVW.size(), 0.0); }
    }
  }
 public:
  std::vector<double> aEnergy;
};

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_points.cpp"
#endif

#endif
