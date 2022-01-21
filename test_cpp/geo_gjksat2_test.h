//
// Created by Nobuyuki Umetani on 2022/01/18.
//

#ifndef GEO_GJKSAT2_TEST_H_
#define GEO_GJKSAT2_TEST_H_

#include <random>

#include "delfem2/geo_gjksat2.h"
#include "delfem2/geo_convhull2.h"

template<typename VEC, typename SCALAR = typename VEC::Scalar>
void TestGjkSat2Test0(unsigned int nitr) {
  std::mt19937 rngeng(std::random_device{}());
  std::uniform_real_distribution<SCALAR> dist_m1p1(-1, +1);
  for (int iitr = 0; iitr < nitr; iitr++) {
    std::vector<VEC> vtxA_xy(10);
    for (auto &p: vtxA_xy) {
      p = {dist_m1p1(rngeng), dist_m1p1(rngeng)};
    }
    std::vector<unsigned int> edgeA_vtx;
    delfem2::ConvexHull2<VEC>(edgeA_vtx, vtxA_xy);
    //
    std::vector<VEC> vtxB_xy0(10);
    for (auto &p: vtxB_xy0) {
      p = {dist_m1p1(rngeng), dist_m1p1(rngeng)};
    }
    std::vector<unsigned int> edgeB_vtx;
    delfem2::ConvexHull2<VEC>(edgeB_vtx, vtxB_xy0);

    for (SCALAR t = 0; t < 3; t += 0.1) {
      std::vector<VEC> vtxB_xy(vtxB_xy0.size());
      for (unsigned int ivtxB = 0; ivtxB < vtxB_xy0.size(); ++ivtxB) { // translate and rotate points
        const SCALAR x0 = vtxB_xy0[ivtxB][0];
        const SCALAR y0 = vtxB_xy0[ivtxB][1];
        const SCALAR x1 = x0 + 2 * std::sin(3 * t);
        const SCALAR y1 = y0;
        const SCALAR x2 = x1 * std::cos(t) - y1 * std::sin(t);
        const SCALAR y2 = x1 * std::sin(t) + y1 * std::cos(t);
        vtxB_xy[ivtxB] = {x2, y2};
      }
      const bool is_intersect_gjk = delfem2::IsIntersect_Points2_Points2_Gjk(vtxA_xy, vtxB_xy);
      EXPECT_EQ(is_intersect_gjk,
                delfem2::IsIntersect_Points2_Points2_Sat(vtxA_xy, vtxB_xy));

      if (!is_intersect_gjk) { continue; }

      VEC normalA;
      delfem2::Penetration_Points2_Points2_Epa(normalA, vtxA_xy, vtxB_xy, 1.0e-5);

      std::vector<VEC> vtxB_xy_translated;

      for (unsigned int ivtxB = 0; ivtxB < vtxB_xy.size(); ++ivtxB) {
        vtxB_xy_translated.push_back(vtxB_xy[ivtxB] + normalA * 1.001);
      }
      EXPECT_EQ(delfem2::IsIntersect_Points2_Points2_Gjk(vtxA_xy, vtxB_xy_translated), false);

      vtxB_xy_translated.clear();
      for (unsigned int ivtxB = 0; ivtxB < vtxB_xy.size(); ++ivtxB) {
        vtxB_xy_translated.push_back(vtxB_xy[ivtxB] + normalA * 0.999);
      }
      EXPECT_EQ(delfem2::IsIntersect_Points2_Points2_Gjk(vtxA_xy, vtxB_xy_translated), true);

    }
  }
}

#endif //GEO_GJKSAT2_TEST_H_
