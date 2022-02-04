#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/fem_rod3_straight.h"
#include "delfem2/sampling.h"
#include "delfem2/geo_vec3.h"

TEST(fem_rod3_straight, check_WdWddW) {
  namespace dfm2 = delfem2;
  std::mt19937 reng(std::random_device{}());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  const double eps = 1.0e-6;
  for (int nitr = 0; nitr < 100; ++nitr) {
    const dfm2::CVec3d vec_pos0[3] = {
        dfm2::RandomVec<3>(dist01, reng),
        dfm2::RandomVec<3>(dist01, reng),
        dfm2::RandomVec<3>(dist01, reng) };
    if ((vec_pos0[1] - vec_pos0[0]).norm() < 0.1) { continue; }
    if ((vec_pos0[2] - vec_pos0[1]).norm() < 0.1) { continue; }
    if ((vec_pos0[1] - vec_pos0[0]).normalized().dot(
        (vec_pos0[2] - vec_pos0[1]).normalized()) < -0.5) {
      continue;
    }
    dfm2::CVec3d dw_dp0[3];
    dfm2::CMat3d ddw_ddp0[3][3];
    double w0 = dfm2::WdWddW_Rod3BendStraight(dw_dp0, ddw_ddp0, vec_pos0);
    for (int ino = 0; ino < 3; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        dfm2::CVec3d vec_pos1[3] = {vec_pos0[0], vec_pos0[1], vec_pos0[2]};
        vec_pos1[ino](idim) += eps;
        dfm2::CVec3d dw_dp1[3];
        dfm2::CMat3d ddw_ddp1[3][3];
        double w1 = dfm2::WdWddW_Rod3BendStraight(dw_dp1, ddw_ddp1, vec_pos1);
//      std::cout << ino << " " << idim << " ## " << (w1 - w0) / eps << " " << dw_dp0[ino](idim) << std::endl;
        double f0 = (w1 - w0) / eps;
        double f1 = dw_dp0[ino](idim);
        EXPECT_NEAR(f0, f1, 1.0e-3 * (1. + fabs(f1)));
        for (int jno = 0; jno < 3; ++jno) {
          for (int jdim = 0; jdim < 3; ++jdim) {
            double d0 = (dw_dp1[jno](jdim) - dw_dp0[jno](jdim)) / eps;
            double d1 = ddw_ddp0[ino][jno].Get(idim, jdim);
            EXPECT_NEAR(d0, d1, 2.e-3 * (1. + fabs(d1)));
//          std::cout << "   " << jno << " " << jdim << " ## " << d0 << " " << d1 << std::endl;
          }
        }
      }
    }
  }
}

TEST(fem_rod3_straight, check_CdC) {
  namespace dfm2 = delfem2;
  std::mt19937 reng(std::random_device{}());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  const double eps = 1.0e-6;
  for (int nitr = 0; nitr < 100; ++nitr) {
    double vec_pos0[3][3];
    delfem2::Fill2dArrayWithRandomValue<3,3>(vec_pos0, dist01, reng);
    if (dfm2::Distance3(vec_pos0[1], vec_pos0[0]) < 0.1) { continue; }
    if (dfm2::Distance3(vec_pos0[2], vec_pos0[1]) < 0.1) { continue; }
    if ((dfm2::CVec3d(vec_pos0[1]) - dfm2::CVec3d(vec_pos0[0])).normalized().dot(
        (dfm2::CVec3d(vec_pos0[2]) - dfm2::CVec3d(vec_pos0[1])).normalized()) < -0.5) {
      continue;
    }
    double c0[3], dc_dp0[3][3][3];
    dfm2::CdC_Rod3BendStraight(c0, dc_dp0, vec_pos0);
    for (int ino = 0; ino < 3; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        double vec_pos1[3][3];
        std::copy_n(&vec_pos0[0][0], 9, &vec_pos1[0][0]);
        vec_pos1[ino][idim] += eps;
        double c1[3], dc_dp1[3][3][3];
        dfm2::CdC_Rod3BendStraight<double>(c1, dc_dp1, vec_pos1);
        for (int jdim = 0; jdim < 3; ++jdim) {
          double f0 = (c1[jdim] - c0[jdim]) / eps;
          double f1 = dc_dp0[jdim][ino][idim];
          EXPECT_NEAR(f0, f1, 1.0e-3 * (1. + fabs(f1)));
        }
      }
    }
  }
}
