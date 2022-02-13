/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"

#include "delfem2/msh_topology_uniform.h"
#include "delfem2/mshmisc.h"
#include "delfem2/gridvoxel.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(gridvoxel, to_tri) {
  dfm2::CGrid3<int> vg;
  vg.Initialize(1, 1, 1, 0);
  vg.Set(0, 0, 0, 1);
  {
    std::vector<double> aXYZ0;
    std::vector<unsigned int> aQuad0;
    dfm2::MeshQuad3D_VoxelGrid(
        aXYZ0, aQuad0,
        vg.ndivx, vg.ndivy, vg.ndivz,
        vg.aVal);
    EXPECT_EQ(aXYZ0.size(), 8 * 3);
    EXPECT_EQ(aQuad0.size(), 6 * 4);
    {
      std::vector<unsigned int> aTri0;
      delfem2::convert2Tri_Quad(aTri0, aQuad0);
      EXPECT_EQ(aTri0.size(), 12 * 3);
    }
  }
}

TEST(gridvoxel, to_quad) {
  dfm2::CGrid3<int> vg;
  vg.Initialize(2, 2, 1, 0);
  vg.Set(0, 0, 0, 1);
  vg.Set(1, 0, 0, 1);
  vg.Set(1, 1, 0, 1);
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aQuad0;
  dfm2::MeshQuad3D_VoxelGrid(
      aXYZ0, aQuad0,
      vg.ndivx, vg.ndivy, vg.ndivz,
      vg.aVal);
  EXPECT_EQ(aXYZ0.size(), 18 * 3);
  EXPECT_EQ(aQuad0.size(), 14 * 4);
  std::vector<double> aXYZ0a;
  std::vector<unsigned int> aQuad0a;
  std::vector<int> mapInOut;
  dfm2::RemoveUnreferencedPoints_MeshElem(
    aXYZ0a, aQuad0a, mapInOut,
    3, aXYZ0, aQuad0);
  EXPECT_EQ(aXYZ0a.size(), 16 * 3);
  EXPECT_EQ(aQuad0a.size(), 14 * 4);
}
