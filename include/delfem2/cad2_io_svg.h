/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_CAD2_IO_SVG_H
#define DFM2_CAD2_IO_SVG_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/cad2_dtri2.h"

namespace delfem2 {

/**
 * @details read an SVG image file and output first path elemnet as a loop of curves.
 * If there is no path element, output first polygon elmenet if they are.
 */
void ReadSVG_LoopEdgeCCad2D(
  std::vector<std::vector<CCad2D_EdgeGeo> > &aaEdge,
  const std::string &fname);

void ReadSVG_Cad2D(
    CCad2D &cad,
    const std::string &fpath,
    double scale);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cad2_io_svg.cpp"
#endif

#endif
