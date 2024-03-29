/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"
#include "delfem2/vec3.h"
#include "delfem2/geo_tri.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/specialfuncs.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

TEST(mathfunc, sherical_harmonics_orthgonality) {
  std::vector<double> aXYZ, aVal;
  std::vector<unsigned int> aTri;
  delfem2::MeshTri3D_Cube(aXYZ, aTri, 50);
  for (int ip = 0; ip < aXYZ.size() / 3; ip++) {
    delfem2::Normalize3(aXYZ.data()+ip*3);
  }
  const int norder = 9;
  const int N = (norder + 1) * (norder + 1);
  double area_sum = 0.0;
  std::vector<double> A(N * N, 0.0);
  for (int it = 0; it < aTri.size() / 3; ++it) {
    const unsigned int i0 = aTri[it * 3 + 0];
    const unsigned int i1 = aTri[it * 3 + 1];
    const unsigned int i2 = aTri[it * 3 + 2];
    dfm2::CVec3d p0(aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]);
    dfm2::CVec3d p1(aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]);
    dfm2::CVec3d p2(aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]);
    double area = delfem2::SolidAngleTri(p0, p1, p2);
    area_sum += area;
    double a0[N];
    dfm2::makeArray_SphericalHarmonics(a0, norder, p0.p[0], p0.p[1], p0.p[2]);
    double a1[N];
    dfm2::makeArray_SphericalHarmonics(a1, norder, p1.p[0], p1.p[1], p1.p[2]);
    double a2[N];
    dfm2::makeArray_SphericalHarmonics(a2, norder, p2.p[0], p2.p[1], p2.p[2]);
    for (int ish = 0; ish < N; ++ish) {
      for (int jsh = 0; jsh < N; ++jsh) {
        double val = 2 * a0[ish] * a0[jsh] + 2 * a1[ish] * a1[jsh] + 2 * a2[ish] * a2[jsh];
        val += a0[ish] * a1[jsh] + a0[ish] * a2[jsh];
        val += a1[ish] * a0[jsh] + a1[ish] * a2[jsh];
        val += a2[ish] * a0[jsh] + a2[ish] * a1[jsh];
        val *= area / 12.0;
        A[N * ish + jsh] += val;
      }
    }
  }
  EXPECT_NEAR(area_sum, M_PI * 4, 1.0e-5);
  for (int ish = 0; ish < N; ++ish) {
    for (int jsh = 0; jsh < N; ++jsh) {
      if (ish == jsh) { continue; }
      EXPECT_NEAR(A[N * ish + jsh], 0.0, 3.0e-3);
    }
  }
  for (int iorder = 0; iorder < norder; ++iorder) {
    for (int ish = iorder * iorder; ish < (iorder + 1) * (iorder + 1); ++ish) {
      if (ish == iorder * (iorder + 1)) {
        EXPECT_NEAR(A[N * ish + ish], 1.0, 1.5e-2);
      } else {
        EXPECT_NEAR(A[N * ish + ish], 0.5, 1.5e-2);
      }
    }
  }
}
