/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file should not depends on anything except for  std::vector, std::functional
 */

#ifndef DFM2_MSH_CENTER_OF_GRAVITY_H
#define DFM2_MSH_CENTER_OF_GRAVITY_H

#include <vector>
#include <functional>  // maybe we should separate the functions with this dependency

#include "delfem2/dfm2_inline.h"

// -----------------
// work on points

namespace delfem2 {


/**
 * @details implemented for "float" and "double"
 */
template<typename T>
DFM2_INLINE void CG_Point3(
    T *center_of_gravity,
    const std::vector<T> &vtx_xyz);    


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
    std::vector<T> &tri_centerxyz,
    const std::vector<T> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx);

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
    T center_of_gravity_xyz[3],
    const std::vector<T> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx);

template<typename T>
DFM2_INLINE void CG_MeshTet3(
    T &v_tot,
    T cg[3],
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTet);

template<typename T>
void CG_MeshLine3(
  T &v_tot,
  T cg[3],
  const std::vector<T> &vtx_xyz,
  const std::vector<unsigned int> &line_vtx);

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_center_of_gravity.cpp"
#endif

#endif /* DFM2_MSHMISC_H */
