//
// Created by Nobuyuki Umetani on 2021-11-19.
//

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/vec2.h"
#include "delfem2/dtri.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/lsmats.h"
#include "delfem2/femmitc3.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/vecxitrsol.h"

TEST(fem_mitc3, check_emat) {
  namespace dfm2 = delfem2;
  std::mt19937 mt(std::random_device{}());
  std::uniform_real_distribution<> dist0(-0.5, +0.5);
  std::uniform_real_distribution<> dist1(+1.0e-10, +1.0);
  for (int itr = 0; itr < 200; ++itr) {
    double C[3][2];
    for (int i = 0; i < 6; ++i) {
      (&C[0][0])[i] = 10.0 * dist0(mt);
    }
    double a0 = dfm2::Area_Tri2(C[0], C[1], C[2]);
    if (a0 < 0.1) continue;
    double u0[3][3];
    for (int i = 0; i < 9; ++i) {
      (&u0[0][0])[i] = 1.0 * dist0(mt);
    }
    double thickness1 = dist1(mt);
    double lambda1 = dist1(mt);
    double myu1 = dist1(mt);
    const double eps = 1.0e-5;
// -------------------
    double W0, dW0[3][3], ddW0[3][3][3][3];
    W0 = 0.0;
    for (int i = 0; i < 9; ++i) { (&dW0[0][0])[i] = 0.0; }
    for (int i = 0; i < 81; ++i) { (&ddW0[0][0][0][0])[i] = 0.0; }
    dfm2::WdWddW_PlateBendingMITC3(W0, dW0, ddW0,
                                   C, u0,
                                   thickness1, lambda1, myu1);
    for (int ino = 0; ino < 3; ++ino) {
      for (int idof = 0; idof < 3; ++idof) {
        double u1[3][3];
        for (int i = 0; i < 9; ++i) { (&u1[0][0])[i] = (&u0[0][0])[i]; }
        u1[ino][idof] += eps;
        double W1, dW1[3][3], ddW1[3][3][3][3];
        W1 = 0.0;
        for (int i = 0; i < 9; ++i) { (&dW1[0][0])[i] = 0.0; }
        for (int i = 0; i < 81; ++i) { (&ddW1[0][0][0][0])[i] = 0.0; }
        dfm2::WdWddW_PlateBendingMITC3(W1, dW1, ddW1,
                                       C, u1,
                                       thickness1, lambda1, myu1);
        EXPECT_NEAR((W1 - W0) / eps,
                    dW0[ino][idof],
                    6.0e-3 * (1.0 + fabs(dW0[ino][idof])));
        for (int jno = 0; jno < 3; ++jno) {
          for (int jdof = 0; jdof < 3; ++jdof) {
            EXPECT_NEAR((dW1[jno][jdof] - dW0[jno][jdof]) / eps,
                        ddW0[ino][jno][idof][jdof],
                        1.0e-3 * (1.0 + fabs(ddW0[ino][jno][idof][jdof])));
          }
        }
      }
    }
  }
}

TEST(fem_mitc3, check_solution_cantilever) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_12(1, +2);
//
  const double lambda = 0.0;
  const double lenx0 = 1.0;
  const double leny0 = 0.2;
  const double thickness0 = 0.05;
  const double myu0 = 10000.0;
  const double rho0 = 1.0;
  const double gravity_z0 = -10.0;
  const double elen0 = 0.03;
  for (int itr = 0; itr < 10; ++itr) {
    const double lenx = lenx0 * dist_12(rndeng);
    const double leny = leny0 * dist_12(rndeng);
    const double thickness = thickness0 * dist_12(rndeng);
    const double myu = myu0 * dist_12(rndeng);
    const double rho = rho0 * dist_12(rndeng);
    const double gravity_z = gravity_z0 * dist_12(rndeng);
    const double elen = elen0 * dist_12(rndeng);
    std::vector<unsigned int> aTri;
    std::vector<double> aXY0;
    {
      std::vector<std::vector<double> > aaXY;
      {
        aaXY.resize(1);
        aaXY[0].push_back(-lenx * 0.5);
        aaXY[0].push_back(-leny * 0.5);
        aaXY[0].push_back(+lenx * 0.5);
        aaXY[0].push_back(-leny * 0.5);
        aaXY[0].push_back(+lenx * 0.5);
        aaXY[0].push_back(+leny * 0.5);
        aaXY[0].push_back(-lenx * 0.5);
        aaXY[0].push_back(+leny * 0.5);
      }
// ---------------------
      std::vector<dfm2::CDynPntSur> aPo2D;
      std::vector<dfm2::CDynTri> aETri;
      std::vector<dfm2::CVec2d> aVec2;
      GenMesh(aPo2D, aETri, aVec2,
              aaXY, elen, elen);
      MeshTri2D_Export(aXY0, aTri,
                       aVec2, aETri);
    }
    std::vector<int> aBCFlag; // boundary condition flag
    {
      const int np = (int) aXY0.size() / 2;
      aBCFlag.assign(np * 3, 0);
      for (int ip = 0; ip < np; ++ip) {
        const double px = aXY0[ip * 2 + 0];
//    const double py = aXY0[ip*2+1];
        if (fabs(px - (-lenx * 0.5)) < 0.0001) {
          aBCFlag[ip * 3 + 0] = 1;
          aBCFlag[ip * 3 + 1] = 1;
          aBCFlag[ip * 3 + 2] = 1;
        }
      }
    }
    dfm2::CMatrixSparse<double> mat_A;
    dfm2::CPreconditionerILU<double> ilu_A;
    {
      std::vector<unsigned int> psup_ind, psup;
      dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                                 aTri.data(), aTri.size() / 3, 3,
                                 (int) aXY0.size() / 2);
      dfm2::JArray_Sort(psup_ind, psup);
//
      const int np = (int) aXY0.size() / 2;
      mat_A.Initialize(np, 3, true);
      mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
      ilu_A.SetPattern0(mat_A);
    }
    std::vector<double> aVal;
    aVal.assign(aXY0.size() / 2 * 3, 0.0);
    std::vector<double> vec_b;
    {
      const size_t np = aXY0.size() / 2;
      const size_t nDoF = np * 3;
// -------------------
      mat_A.setZero();
      vec_b.assign(nDoF, 0.0);
      dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(
        mat_A, vec_b.data(),
        thickness, lambda, myu,
        rho, gravity_z,
        aXY0.data(), aXY0.size() / 2,
        aTri.data(), aTri.size() / 3,
        aVal.data());
      mat_A.SetFixedBC(aBCFlag.data());
      dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
// --------------------------
      std::vector<double> vec_x;
      {
        ilu_A.CopyValue(mat_A);
        ilu_A.Decompose();
        vec_x.resize(vec_b.size());
        std::vector<double> conv;
        {
          const std::size_t n = vec_b.size();
          std::vector<double> tmp0(n), tmp1(n);
          auto vr = dfm2::ViewAsVectorXd(vec_b);
          auto vu = dfm2::ViewAsVectorXd(vec_x);
          auto vs = dfm2::ViewAsVectorXd(tmp0);
          auto vt = dfm2::ViewAsVectorXd(tmp1);
          conv = Solve_PCG(
            vr, vu, vs, vt,
            1.0e-5, 1000, mat_A, ilu_A);
        }
//        std::cout << "convergence   nitr:" << conv.size() << "    res:" << conv[conv.size()-1] << std::endl;
        EXPECT_LT(conv.size(), 1000);
        EXPECT_LT(conv[conv.size() - 1], 1.0e-5);
      }
// -------------------------
      dfm2::XPlusAY(
        aVal,
        static_cast<unsigned int>(nDoF),
        aBCFlag,
        1.0,
        vec_x);
    }
    {
      assert(fabs(lambda) < 1.0e-10);
      const double E = myu * 2.0;
      const double I = thickness * thickness * thickness * leny / 12.0;
      const double W = thickness * lenx * leny * rho * gravity_z;
      const double w = W / lenx;
      const double disp = w * (lenx * lenx * lenx * lenx) / (8.0 * E * I);
//      std::cout << "disp:" << disp << std::endl;
      for (unsigned int ip = 0; ip < aXY0.size() / 2; ++ip) {
        const double px = aXY0[ip * 2 + 0];
        if (fabs(px - (+lenx * 0.5)) > 0.0001) { continue; }
        EXPECT_LE(fabs(aVal[ip * 3 + 0] - disp), 0.002 * fabs(disp));
//        std::cout << aVal[ip*3+0] << " " << disp << "  " << fabs(aVal[ip*3+0] - disp)/disp << std::endl;
      }
    }
  }
}