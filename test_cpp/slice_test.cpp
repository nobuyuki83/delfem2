/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"

#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/points.h"
#include "delfem2/slice.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(slice, test1) {
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CSliceTriMesh> aCS;
  std::vector<std::set<unsigned int> > ReebGraphCS;
  // ----------------------
  delfem2::Read_Ply(
      aXYZ, aTri,
      std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
  delfem2::Normalize_Points3(aXYZ, 1.0);
  std::vector<unsigned int> aTriSuTri;
  dfm2::ElSuEl_MeshElem(
      aTriSuTri,
      aTri.data(), aTri.size() / 3, dfm2::MESHELEM_TRI,
      aXYZ.size() / 3);
  // ----------------------
  std::vector<double> aHeight;
  aHeight.push_back(-0.3);
  aHeight.push_back(-0.2);
  aHeight.push_back(-0.1);
  aHeight.push_back(-0.0);
  aHeight.push_back(+0.1);
  aHeight.push_back(+0.2);
  aHeight.push_back(+0.3);
  const double nrm[3] = {0, 1, 0};
  const double org[3] = {0, 0, 0};
  std::vector<double> aHeightVtx(aXYZ.size() / 3);
  for (unsigned int ip = 0; ip < aXYZ.size() / 3; ++ip) {
    double x0 = aXYZ[ip * 3 + 0] - org[0];
    double y0 = aXYZ[ip * 3 + 1] - org[1];
    double z0 = aXYZ[ip * 3 + 2] - org[2];
    aHeightVtx[ip] = nrm[0] * x0 + nrm[1] * y0 + nrm[2] * z0;
  }
  delfem2::Slice_MeshTri3D_Heights(aCS,
                                   aHeight,
                                   aHeightVtx,
                                   aTri, aTriSuTri);
  MakeReebGraph(ReebGraphCS,
                aCS, aTri, aTriSuTri);
  EXPECT_EQ(aCS.size(), ReebGraphCS.size());
  for (int ics = 0; ics < ReebGraphCS.size(); ++ics) {
    for (auto itr = ReebGraphCS[ics].begin(); itr != ReebGraphCS[ics].end(); ++itr) {
      const unsigned int jcs1 = *itr;
      EXPECT_LT(jcs1, aCS.size());
      EXPECT_EQ(abs(aCS[ics].IndHeight() - aCS[jcs1].IndHeight()), 1);
    }
  }

}
