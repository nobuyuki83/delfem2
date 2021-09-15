/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_CAD2_IO_DXF_H
#define DFM2_CAD2_IO_DXF_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/cad2_dtri2.h"

namespace delfem2 {

/**
 * @brief  write the shape of cad into DXF file
 */
bool WriteCAD_DXF(
    const std::string &file_name,
    const CCad2D &cad,
    double scale);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cad2_io_dxf.cpp"
#endif

#endif
