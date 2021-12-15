/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning
//
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/defarapenergy_geo3.h"
#include "delfem2/jagarray.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/lsmats.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/mshprimitive.h"

namespace dfm2 = delfem2;

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
      dfm2::UpdateRotationsByMatchingCluster_Linear(
        aQuat1,
        aXYZ0, aXYZ1, psup_ind, psup);
    }
  }

  // ------------------
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  std::vector<double> dXYZ12(np * 3);
  for (unsigned int i = 0; i < np * 3; ++i) {
    dXYZ12[i] = dist_m1p1(rndeng) * 0.01;
  }
  assert(aXYZ1.size() == np * 3);
  assert(aQuat1.size() == np * 4);
  double w1 = dfm2::W_ArapEnergy(
    aXYZ0, aXYZ1, aQuat1, psup_ind, psup);

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
  dfm2::dW_ArapEnergy(
    aRes1,
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


