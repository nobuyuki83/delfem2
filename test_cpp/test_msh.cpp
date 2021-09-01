/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"

#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_ioobj.h"
#include "delfem2/msh_iomisc.h"
#include "delfem2/points.h"
#include "delfem2/slice.h"
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
  delfem2::Read_Obj(
      std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
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


TEST(slice,test1){
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CSliceTriMesh> aCS;
  std::vector< std::set<unsigned int> > ReebGraphCS;
  // ----------------------
  delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
                    aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ,1.0);
  std::vector<unsigned int> aTriSuTri;
  dfm2::ElSuEl_MeshElem(aTriSuTri,
                        aTri.data(), aTri.size()/3, dfm2::MESHELEM_TRI,
                        aXYZ.size()/3);
  // ----------------------
  std::vector<double> aHeight;
  aHeight.push_back(-0.3);
  aHeight.push_back(-0.2);
  aHeight.push_back(-0.1);
  aHeight.push_back(-0.0);
  aHeight.push_back(+0.1);
  aHeight.push_back(+0.2);
  aHeight.push_back(+0.3);
  const double nrm[3] = {0,1,0};
  const double org[3] = {0,0,0};
  std::vector<double> aHeightVtx(aXYZ.size()/3);
  for(unsigned int ip=0;ip<aXYZ.size()/3;++ip){
    double x0 = aXYZ[ip*3+0] - org[0];
    double y0 = aXYZ[ip*3+1] - org[1];
    double z0 = aXYZ[ip*3+2] - org[2];
    aHeightVtx[ip] = nrm[0]*x0 + nrm[1]*y0 + nrm[2]*z0;
  }
  delfem2::Slice_MeshTri3D_Heights(aCS,
                                   aHeight,
                                   aHeightVtx,
                                   aTri,aTriSuTri);
  MakeReebGraph(ReebGraphCS,
                aCS, aTri, aTriSuTri);
  EXPECT_EQ( aCS.size(), ReebGraphCS.size() );
  for(int ics=0;ics<ReebGraphCS.size();++ics){
    for(auto itr = ReebGraphCS[ics].begin();itr!=ReebGraphCS[ics].end();++itr){
      const unsigned int jcs1 = *itr;
      EXPECT_LT( jcs1, aCS.size());
      EXPECT_EQ( abs(aCS[ics].IndHeight() - aCS[jcs1].IndHeight()), 1 );
    }
  }

}
