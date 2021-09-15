/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief functions for mesh export/import
 */

#ifndef DFM2_MSH_IO_PLY_H
#define DFM2_MSH_IO_PLY_H

#include <cstdio>
#include <vector>
#include <string>
#include <filesystem>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE void Write_Ply(
    const std::string &fname,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri);

DFM2_INLINE void Write_Ply(
    const std::string &fname,
    unsigned int nXYZ, double *paXYZ,
    unsigned int nTri, unsigned int *paTri);

DFM2_INLINE void Read_Ply(
    const std::string &fname,
    int &nnode_,
    double *&pXYZs_,
    int &ntri_,
    unsigned int *&aTri_);

DFM2_INLINE void Read_Ply(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    const std::filesystem::path &file_path);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_io_ply.cpp"
#endif

#endif // DFM2_MSH_IO_PLY_H
