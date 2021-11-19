//
// Created by Nobuyuki Umetani on 2021-11-20.
//

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/fem_quadratic_bending.h"
#include "delfem2/vec3.h"

bool RandomTriangle(double P[4][3]){
  namespace dfm2 = delfem2;
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1., +1.);
  for(unsigned int ip=0;ip<4;++ip){
    for(unsigned int idim=0;idim<3;++idim){
      P[ip][idim] = dist_m1p1(randomEng);
    }
  }
  if (dfm2::Distance3(P[0], P[1]) < 0.1) { return false; }
  if (dfm2::Distance3(P[0], P[2]) < 0.1) { return false; }
  if (dfm2::Distance3(P[0], P[3]) < 0.1) { return false; }
  if (dfm2::Distance3(P[1], P[2]) < 0.1) { return false; }
  if (dfm2::Distance3(P[1], P[3]) < 0.1) { return false; }
  if (dfm2::Distance3(P[2], P[3]) < 0.1) { return false; }
  if (dfm2::Area_Tri3(P[0], P[2], P[3]) < 0.01) { return false; }
  if (dfm2::Area_Tri3(P[1], P[2], P[3]) < 0.01) { return false; }
  return true;
}

TEST(fem_quadratic_bending, Check_CdC) {
  namespace dfm2 = delfem2;
  const double eps = 1.0e-5;
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    double P[4][3];
    if( !RandomTriangle(P) ){ continue; }
    double p0[4][3];
    if( !RandomTriangle(p0) ){ continue; }
    double C0[3], dC[3][4][3];
    dfm2::CdC_QuadBend(C0, dC, P, p0);
    for (int ino = 0; ino < 4; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        double p1[4][3];
        std::copy_n(&p0[0][0], 12, &p1[0][0]);
        p1[ino][idim] += eps;
        double C1[3];
        {
          double dC1[3][4][3];
          dfm2::CdC_QuadBend(C1, dC1, P, p1);
        }
        for(int jdim=0;jdim<3;++jdim) {
          const double val0 = (C1[jdim] - C0[jdim]) / eps;
          const double val1 = dC[jdim][ino][idim];
          EXPECT_NEAR(val0, val1, (1.0 + fabs(val1)) * 1.5e-2);
        }
      }
    }
  }
}