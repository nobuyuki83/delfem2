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
