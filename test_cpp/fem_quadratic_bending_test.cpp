//
// Created by Nobuyuki Umetani on 2021-11-20.
//

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/fem_quadratic_bending.h"
#include "delfem2/vec3.h"

bool RandomTriangle(
  double P[4][3],
  bool is_avoid_degenerated) {
  namespace dfm2 = delfem2;
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1., +1.);
  for (unsigned int ip = 0; ip < 4; ++ip) {
    for (unsigned int idim = 0; idim < 3; ++idim) {
      P[ip][idim] = dist_m1p1(randomEng);
    }
  }
  if (!is_avoid_degenerated) { return true; }
  if (dfm2::Distance3(P[0], P[1]) < 0.1) { return false; }
  if (dfm2::Distance3(P[0], P[2]) < 0.1) { return false; }
  if (dfm2::Distance3(P[0], P[3]) < 0.1) { return false; }
  if (dfm2::Distance3(P[1], P[2]) < 0.1) { return false; }
  if (dfm2::Distance3(P[1], P[3]) < 0.1) { return false; }
  if (dfm2::Distance3(P[2], P[3]) < 0.1) { return false; }
  if (dfm2::Area_Tri3(P[0], P[2], P[3]) < 0.01) { return false; }
  if (dfm2::Area_Tri3(P[1], P[2], P[3]) < 0.01) { return false; }
  double A0, UN0[3];
  dfm2::UnitNormalAreaTri3(UN0, A0, P[0], P[2], P[3]);
  double A1, UN1[3];
  dfm2::UnitNormalAreaTri3(UN1, A1, P[1], P[3], P[2]);
  const double L0 = dfm2::Distance3(P[2], P[3]);
  const double H0 = A0 * 2.0 / L0;
  const double H1 = A1 * 2.0 / L0;
  if( H0 < 0.4 ){ return false; }
  if( H1 < 0.4 ){ return false; }
  return true;
}

TEST(fem_quadratic_bending, Check_CdC) {
  namespace dfm2 = delfem2;
  const double eps = 1.0e-5;
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    double P[4][3];
    if (!RandomTriangle(P,true)) { continue; }
    double p0[4][3];
    if (!RandomTriangle(p0, false)) { continue; }
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
        for (int jdim = 0; jdim < 3; ++jdim) {
          const double val0 = (C1[jdim] - C0[jdim]) / eps;
          const double val1 = dC[jdim][ino][idim];
          EXPECT_NEAR(val0, val1, (1.0 + fabs(val1)) * 1.5e-2);
        }
      }
    }
  }
}


TEST(fem_quadratic_bending, Check_dRdC) {
  namespace dfm2 = delfem2;
  const double eps = 1.0e-5;
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    double P0[4][3];
    if (!RandomTriangle(P0, true)) { continue; }
    double p[4][3];
    if (!RandomTriangle(p, true)) { continue; }
    double stiff = 1.;
    double Res0[4][3];
    double dRdC0[4][4][3][3];
    {
      double Kmat[4][4][3][3];
      dfm2::WdWddW_QuadraticBending_Sensitivity(Kmat,Res0,dRdC0,P0,p,stiff);
    }
    for(unsigned int ip=0;ip<4;++ip){
      for(unsigned int idim=0;idim<3;++idim) {
        double P1[4][3];
        std::copy_n(&P0[0][0], 12, &P1[0][0]);
        P1[ip][idim] += eps;
        double Res1[4][3];
        {
          double dRdC1[4][4][3][3];
          double Kmat[4][4][3][3];
          dfm2::WdWddW_QuadraticBending_Sensitivity(Kmat,Res1,dRdC1,P1,p,stiff);
        }
        for(unsigned int jp=0;jp<4;++jp){
          for(unsigned int jdim=0;jdim<3;++jdim) {
            double v0 = (Res1[jp][jdim]  - Res0[jp][jdim]) / eps;
            double v1 = dRdC0[jp][ip][jdim][idim];
            EXPECT_NEAR(v0,v1,4.0e-2*(1+fabs(v1)));
          }
        }
      }
    }
  }
}
