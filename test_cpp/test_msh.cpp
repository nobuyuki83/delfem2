/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"

#include "delfem2/mshtopo.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/gridvoxel.h"
#include <cstring>
#include <random>

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(mshio,load_obj)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Obj(std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
                    aXYZ, aTri);
  EXPECT_EQ(aTri.size(),1000*3);
}

TEST(meshtopo,quad_subdiv0)
{
  dfm2::CGrid3<int> vg;
  vg.Initialize(1,1,1, 0);
  vg.Set(0,0,0, 1);
  {
    std::vector<double> aXYZ0;
    std::vector<unsigned int> aQuad0;
    dfm2::MeshQuad3D_VoxelGrid(aXYZ0, aQuad0,
                               vg.ndivx, vg.ndivy, vg.ndivz,
                               vg.aVal);
    EXPECT_EQ(aXYZ0.size(),8*3);
    EXPECT_EQ(aQuad0.size(),6*4);
    {
      std::vector<unsigned int> aTri0;
      delfem2::convert2Tri_Quad(aTri0, aQuad0);
      EXPECT_EQ(aTri0.size(),12*3);
    }
  }
}

TEST(meshtopo,quad_subdiv1)
{
  dfm2::CGrid3<int> vg;
  vg.Initialize(2,2,1, 0);
  vg.Set(0,0,0, 1);
  vg.Set(1,0,0, 1);
  vg.Set(1,1,0, 1);
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aQuad0;
  dfm2::MeshQuad3D_VoxelGrid(aXYZ0, aQuad0,
                             vg.ndivx, vg.ndivy, vg.ndivz,
                             vg.aVal);
  EXPECT_EQ(aXYZ0.size(),18*3);
  EXPECT_EQ(aQuad0.size(),14*4);
  std::vector<double> aXYZ0a;
  std::vector<unsigned int> aQuad0a;
  std::vector<int> mapInOut;
  dfm2::RemoveUnreferencedPoints_MeshElem(aXYZ0a,aQuad0a, mapInOut,
                                          3,aXYZ0,aQuad0);

  EXPECT_EQ(aXYZ0a.size(),16*3);
  EXPECT_EQ(aQuad0a.size(),14*4);
}