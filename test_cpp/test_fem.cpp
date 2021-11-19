/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>
#include <cstring>

#include "gtest/gtest.h" // need to be defiend in the beginning
//
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/femmips_geo3.h"
#include "delfem2/defarapenergy_geo3.h"
#include "delfem2/pbd_geo3.h"
#include "delfem2/jagarray.h"
#include "delfem2/vec2.h"

#include "delfem2/vecxitrsol.h"
#include "delfem2/lsmats.h"

#include "delfem2/fempoisson.h"
#include "delfem2/fem_distance3.h"

#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/mshprimitive.h"

namespace dfm2 = delfem2;

// --------------------------------------

TEST(femem2, poisson_quad) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist(0.01, 3);
  double lx = dist(randomEng);
  double ly = dist(randomEng);
  double em0[4][4];
  dfm2::EMat_Poission2_QuadOrth(em0, lx, ly);
  for (unsigned int ngauss = 1; ngauss < 3; ++ngauss) {
    double em1[4][4];
    dfm2::EMat_Poisson2_QuadOrth_GaussInt(em1, lx, ly, ngauss);
    double diff = 0.0;
    for (unsigned int i = 0; i < 16; ++i) {
      double v0 = (&em0[0][0])[i];
      double v1 = (&em1[0][0])[i];
      diff += (v0 - v1) * (v0 - v1);
    }
    EXPECT_NEAR(diff, 0.0, 1.0e-10);
  }
}

TEST(objfunc_v23, Check_CdC_TriStrain) {
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m5p5(-5, 5);
  // -----
  for (int itr = 0; itr < 200; ++itr) {
    const double P[3][2] = {
        {dist_m5p5(randomEng), dist_m5p5(randomEng)},
        {dist_m5p5(randomEng), dist_m5p5(randomEng)},
        {dist_m5p5(randomEng), dist_m5p5(randomEng)}};
    if (dfm2::Distance2(P[0], P[1]) < 0.1) { continue; }
    if (dfm2::Distance2(P[1], P[2]) < 0.1) { continue; }
    if (dfm2::Distance2(P[0], P[2]) < 0.1) { continue; }
    double a0 = dfm2::Area_Tri2(P[0], P[1], P[2]);
    if (fabs(a0) < 0.2) continue;
    const double p[3][3] = {
        {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng)},
        {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng)},
        {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng)}};
    const double eps = 1.0e-5;
    double C[3], dCdp[3][9];
    dfm2::PBD_CdC_TriStrain2D3D(C, dCdp, P, p);
    for (int ine = 0; ine < 3; ++ine) {
      for (int idim = 0; idim < 3; ++idim) {
        double p1[3][3];
        for (int i = 0; i < 9; ++i) { (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1[3], dCdp1[3][9];
        dfm2::PBD_CdC_TriStrain2D3D(C1, dCdp1, P, p1);
        EXPECT_NEAR((C1[0] - C[0]) / eps, dCdp[0][ine * 3 + idim], 1.0e-2);
        EXPECT_NEAR((C1[1] - C[1]) / eps, dCdp[1][ine * 3 + idim], 1.0e-2);
        EXPECT_NEAR((C1[2] - C[2]) / eps, dCdp[2][ine * 3 + idim], 1.0e-2);
      }
    }
  }
}

TEST(pbd, Check_CdC_DiscreteShell) {
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  const double eps = 1.0e-5;
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    const dfm2::CVec3d p[4] = {
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng)};
    if (dfm2::Distance(p[0], p[1]) < 0.1) { continue; }
    if (dfm2::Distance(p[0], p[2]) < 0.1) { continue; }
    if (dfm2::Distance(p[0], p[3]) < 0.1) { continue; }
    if (dfm2::Distance(p[1], p[2]) < 0.1) { continue; }
    if (dfm2::Distance(p[1], p[3]) < 0.1) { continue; }
    if (dfm2::Distance(p[2], p[3]) < 0.1) { continue; }
    if (dfm2::Area_Tri(p[0], p[2], p[3]) < 0.01) { continue; }
    if (dfm2::Area_Tri(p[1], p[2], p[3]) < 0.01) { continue; }
    double C;
    dfm2::CVec3d dC[4];
    CdC_DiscreteShell(C, dC, p[0], p[1], p[2], p[3]);
    for (int ino = 0; ino < 4; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        dfm2::CVec3d p1[4] = {p[0], p[1], p[2], p[3]};
        p1[ino][idim] += eps;
        double C1;
        dfm2::CVec3d dC1[4];
        CdC_DiscreteShell(C1, dC1, p1[0], p1[1], p1[2], p1[3]);
        const double val0 = (C1 - C) / eps;
        const double val1 = dC[ino][idim];
        EXPECT_NEAR(val0, val1, (1.0 + fabs(val1)) * 1.5e-2);
      }
    }
  }
}

TEST(objfunc_v23, MIPS) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0, 1);
  const double eps = 1.0e-5;
  //
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    const double P[3][3] = {
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)}};
    if (dfm2::Distance3(P[0], P[1]) < 0.1) { continue; }
    if (dfm2::Distance3(P[0], P[2]) < 0.1) { continue; }
    if (dfm2::Distance3(P[1], P[2]) < 0.1) { continue; }
    if (dfm2::Area_Tri3(P[0], P[1], P[2]) < 0.01) { continue; }
    double p[3][3];
    {
      dfm2::CMat3d m;
      m.SetRotMatrix_Cartesian(0.3, 1.0, 0.5);
      for (int ino = 0; ino < 3; ++ino) {
        auto vec = m.MatVec(P[ino]);
        p[ino][0] = vec[0];
        p[ino][1] = vec[1];
        p[ino][2] = vec[2];
      }
    }
    if (dfm2::Area_Tri3(p[0], p[1], p[2]) < 0.01) { continue; }
    double E, dE[3][3], ddE[3][3][3][3];
    dfm2::WdWddW_MIPS(E, dE, ddE,
                      p, P);
    for (int ino = 0; ino < 3; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        double c1[3][3] = {
            {p[0][0], p[0][1], p[0][2]},
            {p[1][0], p[1][1], p[1][2]},
            {p[2][0], p[2][1], p[2][2]}};
        c1[ino][idim] += eps;
        double E1, dE1[3][3], ddE1[3][3][3][3];
        dfm2::WdWddW_MIPS(E1, dE1, ddE1,
                          c1, P);
        {
          const double val0 = (E1 - E) / eps;
          const double val1 = dE[ino][idim];
          EXPECT_NEAR(val0, val1, 1.0e-2 * (1 + fabs(val1)));
        }
        for (int jno = 0; jno < 3; ++jno) {
          for (int jdim = 0; jdim < 3; ++jdim) {
            const double val0 = (dE1[jno][jdim] - dE[jno][jdim]) / eps;
            const double val1 = ddE[jno][ino][jdim][idim];
            EXPECT_NEAR(val0, val1, 3.e-2 * (1 + fabs(val1)));
          }
        }
      }
    }
  }
}

TEST(objfunc_v23, distancetri2d3d) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0, 1);
  //
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    const double P[3][2] = {  // undeformed triangle vertex positions
        {dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng)}};
    if (dfm2::Distance3(P[0], P[1]) < 0.1) { continue; }
    if (dfm2::Distance3(P[0], P[2]) < 0.1) { continue; }
    if (dfm2::Distance3(P[1], P[2]) < 0.1) { continue; }
    const double p[3][3] = { //  deformed triangle vertex positions)
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)}};
    // --------------
    double C[3], dCdp[3][9];
    dfm2::PBD_ConstraintProjection_DistanceTri2D3D(C, dCdp, P, p);
    for (int ine = 0; ine < 3; ++ine) {
      for (int idim = 0; idim < 3; ++idim) {
        double eps = 1.0e-6;
        double p1[3][3];
        for (int i = 0; i < 9; ++i) { (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1[3], dCdp1[3][9];
        dfm2::PBD_ConstraintProjection_DistanceTri2D3D(C1, dCdp1, P, p1);
        for (int jdim = 0; jdim < 3; ++jdim) {
          const double val0 = (C1[jdim] - C[jdim]) / eps;
          const double val1 = dCdp[jdim][ine * 3 + idim];
          EXPECT_NEAR(val0, val1, 1.0e-2 * (fabs(val1) + 1.0));
        }
      }
    }
  }
}

TEST(objfunc_v23, pbd_energy_stvk) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0, 1);
  std::uniform_real_distribution<double> dist_12(1, 2);

  //
  for (int itr = 0; itr < 1000; ++itr) {
    const double P[3][2] = {  // undeformed triangle vertex positions
        {dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng)}};
    const double A0 = dfm2::Area_Tri2(P[0], P[1], P[2]);
    if (A0 < 0.1) continue;
    const double p[3][3] = { //  deformed triangle vertex positions)
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
        {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)}};
    const double lambda = dist_12(randomEng);
    const double myu = dist_12(randomEng);
    // ---------------------
    double C, dCdp[9];
    dfm2::PBD_ConstraintProjection_EnergyStVK(C, dCdp, P, p, lambda, myu);
    for (int ine = 0; ine < 3; ++ine) {
      for (int idim = 0; idim < 3; ++idim) {
        double eps = 1.0e-8;
        double p1[3][3];
        for (int i = 0; i < 9; ++i) { (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1, dCdp1[9];
        dfm2::PBD_ConstraintProjection_EnergyStVK(C1, dCdp1, P, p1, lambda, myu);
        EXPECT_NEAR((C1 - C) / eps, dCdp[ine * 3 + idim], 1.0e-3 * (1 + fabs(dCdp[ine * 3 + idim])));
      }
    }
  }
}

TEST(objfunc_v23, WdWddW_SquareLengthLineseg3D) {
  std::random_device rd;
  std::mt19937 rndeng(rd());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  std::uniform_real_distribution<double> dist_12(+1, +2);
  std::uniform_real_distribution<double> dist_01(+0, +1);
  double eps = 1.0e-5;
  //
  for (int itr = 0; itr < 10000; ++itr) {
    const double stiff_stretch = dist_12(rndeng);
    const dfm2::CVec3d P[2] = {
        dfm2::CVec3d::Random(dist_01, rndeng),
        dfm2::CVec3d::Random(dist_01, rndeng)};
    if ((P[0] - P[1]).norm() < 0.1) { continue; }
    dfm2::CVec3d dW_dP[2];
    dfm2::CMat3d ddW_ddP[2][2];
    const double L0 = 1.0;
    double W = WdWddW_SquareLengthLineseg3D(
        dW_dP, ddW_ddP,
        stiff_stretch, P, L0);
    // -----
    const dfm2::CVec3d dP[2] = {
        dfm2::CVec3d::Random(dist_01, rndeng) * eps,
        dfm2::CVec3d::Random(dist_01, rndeng) * eps};
    const dfm2::CVec3d p[2] = {P[0] + dP[0], P[1] + dP[1]};
    double w;
    dfm2::CVec3d dw_dP[2];
    {
      dfm2::CMat3d ddw_ddP[2][2];
      w = WdWddW_SquareLengthLineseg3D(
          dw_dP, ddw_ddP,
          stiff_stretch, p, L0);
    }
    {
      const double val0 = (w - W) / eps;
      const double val1 = (+dW_dP[0].dot(dP[0]) + dW_dP[1].dot(dP[1])) / eps;
      EXPECT_NEAR(val0, val1, 1.0e-3 * (1 + fabs(val1)));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[0] - dW_dP[0]) / eps;
      const dfm2::CVec3d val1 = (+ddW_ddP[0][0] * dP[0] + ddW_ddP[0][1] * dP[1]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-3 * (1 + val1.norm()));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[1] - dW_dP[1]) / eps;
      const dfm2::CVec3d val1 = (+ddW_ddP[1][0] * dP[0] + ddW_ddP[1][1] * dP[1]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-3 * (1 + val1.norm()));
    }
  }
}

TEST(objfunc_v23, arap) {
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3D_Cube(aXYZ0, aTri, 10);
  const size_t np = aXYZ0.size() / 3;

  std::vector<unsigned int> psup_ind, psup;
  {
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        aTri.data(), aTri.size() / 3, 3,
        aXYZ0.size() / 3);
    dfm2::JArray_Sort(psup_ind, psup);
  }

  // base mesh
  // ------------------
  // deformed mesh

  std::vector<double> aXYZ1 = aXYZ0;
  dfm2::Translate_Points3(aXYZ1, 0.1, 0.2, 0.3);
  dfm2::Rotate_Points3(aXYZ1, 0.1, 0.2, 0.3);

  std::vector<double> aQuat1;
  { // initialize rotation
    aQuat1.resize(np * 4);
    for (unsigned int ip = 0; ip < np; ++ip) {
      dfm2::Quat_Identity(aQuat1.data() + 4 * ip);
    }
    for (int itr = 0; itr < 40; ++itr) {
      dfm2::UpdateRotationsByMatchingCluster_Linear(aQuat1,
                                                    aXYZ0, aXYZ1, psup_ind, psup);
    }
  }

  // ------------------

  std::random_device rd;
  std::mt19937 rndeng(rd());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  std::vector<double> dXYZ12(np * 3);
  for (unsigned int i = 0; i < np * 3; ++i) {
    dXYZ12[i] = dist_m1p1(rndeng) * 0.01;
  }
  assert(aXYZ1.size() == np * 3);
  assert(aQuat1.size() == np * 4);
  double w1 = dfm2::W_ArapEnergy(aXYZ0, aXYZ1, aQuat1, psup_ind, psup);

  const double eps = 1.0e-5;
  std::vector<double> aXYZ2 = aXYZ1;
  for (int i = 0; i < aXYZ2.size(); ++i) { aXYZ2[i] += eps * dXYZ12[i]; }

  std::vector<double> aQuat2 = aQuat1;
  dfm2::UpdateRotationsByMatchingCluster_Linear(
      aQuat2,
      aXYZ0, aXYZ2, psup_ind, psup);

  double w2 = dfm2::W_ArapEnergy(
      aXYZ0, aXYZ2, aQuat2, psup_ind, psup);
  std::vector<double> aRes2;
  dfm2::dW_ArapEnergy(
      aRes2,
      aXYZ0, aXYZ2, aQuat2, psup_ind, psup);

  // ---------------------------------

  std::vector<double> aRes1;
  dfm2::dW_ArapEnergy(aRes1,
                      aXYZ0, aXYZ1, aQuat1, psup_ind, psup);

  {
    double dw = dfm2::DotX(aRes1.data(), dXYZ12.data(), aRes1.size());
    double val2 = (w2 - w1) / eps;
    EXPECT_NEAR(dw, val2, 1.0e-5 * (fabs(dw) + 1.0));
  }

  // ---------------------------------

  dfm2::CMatrixSparse<double> Mat;
  {
    Mat.Initialize(static_cast<unsigned int>(np), 3, true);
    std::vector<unsigned int> psup_ind1, psup1;
    dfm2::JArray_Extend(psup_ind1, psup1,
                        psup_ind.data(), psup_ind.size(), psup.data());
    dfm2::JArray_Sort(psup_ind1, psup1);
    assert(psup_ind1.size() == np + 1);
    Mat.SetPattern(psup_ind1.data(), psup_ind1.size(), psup1.data(), psup1.size());
    Mat.setZero();
    std::vector<unsigned int> tmp_buffer;
    for (unsigned int ip = 0; ip < np; ++ip) {
      std::vector<unsigned int> aIP;
      for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
        const unsigned int jp = psup[ipsup];
        aIP.push_back(jp);
      }
      aIP.push_back(ip);
      std::vector<double> eM;
      dfm2::ddW_ArapEnergy(eM,
                           aIP, aXYZ0, aQuat1);
      Mearge(
          Mat,
          aIP.size(), aIP.data(),
          aIP.size(), aIP.data(),
          9, eM.data(), tmp_buffer);
    }
  }
  EXPECT_LT(dfm2::CheckSymmetry(Mat), 1.0e-10);

  {
    std::vector<double> aRes12(np * 3);
    Mat.MatVec(aRes12.data(),
               1.0, dXYZ12.data(), 0.0);

    for (int i = 0; i < np * 3; ++i) {
      double val0 = (aRes2[i] - aRes1[i]) / eps;
      double val1 = aRes12[i];
//      std::cout << i << " " << val0 << " " << val1 << std::endl;
      EXPECT_NEAR(val0, val1, 1.0e-5 * (1 + fabs(val1)));
    }
  }

}

// ------------------------------------------------------------


