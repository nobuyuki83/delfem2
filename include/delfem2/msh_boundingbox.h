/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file should not depends on anything except for  std::vector, std::functional
 */

#ifndef DFM2_MSH_BOUNDINGBOX_H
#define DFM2_MSH_BOUNDINGBOX_H

#include <vector>
#include <functional>  // maybe we should separate the functions with this dependency

#include "delfem2/dfm2_inline.h"

// -----------------
// work on points

namespace delfem2 {


/**
 * @brief update minimum and maximum coordinates
 * @details implemented for "float" and "double"
 */
template<typename T>
DFM2_INLINE void updateMinMaxXYZ(
    T &x_min, T &x_max,
    T &y_min, T &y_max,
    T &z_min, T &z_max,
    T x, T y, T z);

/**
 * @param bb3 (out) bounding box in the order of <minx, miny, minz, maxx, maxy, maxz>
 * @param vtx_xyz (in) array of 3D coordinates of points
 * @param num_vtx (in) number of points
 * @details implemented for "float" and "double"
 */
template<typename T>
DFM2_INLINE void BoundingBox3_Points3(
    T min3[3],
    T max3[3],
    const T *vtx_xyz,
    const size_t num_vtx);

// center & width
template<typename T>
DFM2_INLINE void CenterWidth_Point3(
    T &cx, T &cy, T &cz,
    T &wx, T &wy, T &wz,
    const T *vtx_xyz,
    const size_t num_vtx);

template<typename T>
DFM2_INLINE void CenterWidth_Points3(
    T &cx, T &cy, T &cz,
    T &wx, T &wy, T &wz,
    const std::vector<T> &vtx_xyz);

template<typename T>
DFM2_INLINE void CenterWidth_Points3(
    T c[3],
    T w[3],
    const std::vector<T> &vtx_xyz);

DFM2_INLINE double Size_Points3D_LongestAABBEdge(
    const std::vector<double> &aXYZ);

// local coordinate
DFM2_INLINE void GetCenterWidthLocal(
    double &lcx, double &lcy, double &lcz,
    double &lwx, double &lwy, double &lwz,
    const std::vector<double> &aXYZ,
    const double lex[3],
    const double ley[3],
    const double lez[3]);    

DFM2_INLINE void GetCenterWidthGroup(
    double &cx, double &cy, double &cz,
    double &wx, double &wy, double &wz,
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &elm_vtx,
    unsigned int nnoel,
    int igroup,
    const std::vector<int> &elm_groupidx);

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

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_boundingbox.cpp"
#endif

#endif /* DFM2_MSHMISC_H */
