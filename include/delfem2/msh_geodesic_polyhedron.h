/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_MSH_GEODESIC_POLYHEDRON_H
#define DFM2_MSH_GEODESIC_POLYHEDRON_H

#include <string>
#include <vector>
#include <fstream>

void MeshTri3D_GeodesicPolyhedron(
    std::vector<double>& aXYZ1,
    std::vector<unsigned int>& aTri1);

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_geodesic_polyhedron.cpp"
#endif

#endif
