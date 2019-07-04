/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef gl_voxsdf_h
#define gl_voxsdf_h

#include <stdio.h>

#include "delfem2/voxel.h"
#include "delfem2/sdf.h"
#include "delfem2/bv.h"

void Draw_CubeGrid(bool is_picked, int iface_picked,
                   double elen, const CVector3& org,
                   const CCubeGrid& cube);


#endif /* gl_voxsdf_hpp */
