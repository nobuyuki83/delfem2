/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"
#include <cstring>
#include <random>
#include "delfem2/lsvecx.h"
#include "delfem2/lsmats.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/dtri.h"
#include "delfem2/eigen/ls_dense.h"
#include "delfem2/eigen/ls_sparse.h"
#include "delfem2/eigen/ls_ilu_sparse.h"
#include "delfem2/lsitrsol.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Core>

// ------------------------------------------

void MakeMesh(
    std::vector<double>& aXY1,
    std::vector<unsigned int>& aTri1,
    std::vector<int>& aBCFlag,
    double elen,
    unsigned int ndimval)
{
  std::vector< std::vector<double> > aaXY;
  const double len = 1.0;
  {
    aaXY.resize(1);
    aaXY[0].push_back(-len); aaXY[0].push_back(-len);
    aaXY[0].push_back(-len); aaXY[0].push_back(+len);
    aaXY[0].push_back(+len); aaXY[0].push_back(+len);
    aaXY[0].push_back(+len); aaXY[0].push_back(-len);
  }
  std::vector<delfem2::CDynPntSur> aPo2D;
  std::vector<delfem2::CDynTri> aETri;
  std::vector<delfem2::CVec2d> aVec2;
  delfem2::GenMesh(aPo2D,aETri,aVec2,
                   aaXY,elen,elen);
  MeshTri2D_Export(
      aXY1,aTri1,
      aVec2,aETri);
  const unsigned int np = aXY1.size()/2;
  aBCFlag.assign(np*ndimval, 0);
  for(unsigned int ip=0;ip<np;++ip){
//    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip*2+1];
    if( fabs(py-len) > 0.0001 ){ continue; }
    for(unsigned int idim=0;idim<ndimval;++idim) {
      aBCFlag[ip * 2 + idim] = 1;
    }
  }
//  std::cout<<"  ntri;"<<aTri1.size()/3<<"  nXY:"<<aXY1.size()/2<<std::endl;
}

TEST(ls,test1)
{
  namespace dfm2 = delfem2;
  const double epsilon = 1.0e-5;
  //
  std::vector<unsigned int> aTri1;
  std::vector<double> aXY1;
  std::vector<int> aBCFlag;
  MakeMesh(
      aXY1,aTri1,aBCFlag,
      0.015, 2);
  const unsigned int np = aXY1.size()/2;
  const unsigned int nDoF = np*2;
  const std::vector<double> aVal(nDoF,0);
  // -----------
  delfem2::CMatrixSparseBlock<Eigen::Matrix2d,Eigen::aligned_allocator<Eigen::Matrix2d>> mA0;
  Eigen::VectorXd vb0(nDoF);
  dfm2::CMatrixSparse<double> mA1;
  std::vector<double> vb1(nDoF);
  {
    std::vector<unsigned int> psup_ind0, psup0;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind0, psup0,
        aTri1.data(), aTri1.size() / 3, 3,
        aXY1.size() / 2);
    dfm2::JArray_Sort(psup_ind0, psup0);
    mA0.Initialize(np);
    mA0.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
    const double myu = 10.0, lambda = 10.0, rho = 1.0, g_x = 0.0, g_y = -3.0;
    mA0.setZero();
    vb0.setZero();
    dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
        mA0, vb0.data(),
        myu, lambda, rho, g_x, g_y,
        aXY1.data(), aXY1.size() / 2,
        aTri1.data(), aTri1.size() / 3,
        aVal.data());
    SetFixedBC_Dia(mA0, aBCFlag.data(), 1.f);
    SetFixedBC_Col(mA0, aBCFlag.data());
    SetFixedBC_Row(mA0, aBCFlag.data());
    delfem2::setZero_Flag(vb0, aBCFlag, 0);
    // ----------
    mA1.Initialize(np, 2, true);
    mA1.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
    mA1.setZero();
    vb1.assign(nDoF,0.0);
    dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
        mA1, vb1.data(),
        myu, lambda, rho, g_x, g_y,
        aXY1.data(), aXY1.size() / 2,
        aTri1.data(), aTri1.size() / 3,
        aVal.data());
    mA1.SetFixedBC(aBCFlag.data());
    dfm2::setRHS_Zero(vb1, aBCFlag, 0);
  }
  // ---------------
  double conv_ratio = 1.0e-6;
  int iteration = 1000;
  {
    const auto time0 = std::chrono::system_clock::now();
    unsigned int nitr1 = 0;
    for(int itr=0;itr<10;++itr){ // CG method std::vector
      std::vector<double> vx1(vb1.size());
      const std::size_t n = vb1.size();
      std::vector<double> tmp0(n), tmp1(n);
      std::vector<double> aConv = Solve_CG(
          dfm2::CVecXd(vb1),
          dfm2::CVecXd(vx1),
          dfm2::CVecXd(tmp0),
          dfm2::CVecXd(tmp1),
          conv_ratio, iteration, mA1);
      nitr1 = aConv.size();
    }
    const auto time1 = std::chrono::system_clock::now();
    // ---------------
    unsigned int nitr0 = 0;
    for(int itr=0;itr<10;++itr){ // CG method with eigen
      Eigen::VectorXd vx0(vb0.size());
      const std::size_t n = vb0.size();
      Eigen::VectorXd tmp0(n), tmp1(n);
      std::vector<double> aConv = delfem2::Solve_CG(
          vb0, vx0, tmp0, tmp1,
          conv_ratio, iteration, mA0);
      nitr0 = aConv.size();
    }
    EXPECT_EQ(nitr0,nitr1);
    const auto time2 = std::chrono::system_clock::now();
    double elapsed01 = std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count();
    double elapsed12 = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
    std::cout << "cg std::vector: " << elapsed01 << "   cg eigen: " << elapsed12 << std::endl;
  }
  delfem2::CILU_SparseBlock<Eigen::Matrix2d,Eigen::aligned_allocator<Eigen::Matrix2d>> ilu0;
  delfem2::CPreconditionerILU<double> ilu1;
  { // LU
    ilu1.Initialize_ILU0(mA1);
    ilu1.SetValueILU(mA1);
    ilu1.DoILUDecomp();
    delfem2::ILU_SetPattern0(ilu0, mA0);
    delfem2::ILU_CopyValue(ilu0, mA0);
    delfem2::ILU_Decompose(ilu0);
    // check if the entry is the same
    for(unsigned int icrs=0;icrs<ilu0.mat.valCrs.size();++icrs) {
      EXPECT_NEAR(ilu0.mat.valCrs[icrs](0, 0), ilu1.mat.valCrs[icrs * 4 + 0], epsilon);
      EXPECT_NEAR(ilu0.mat.valCrs[icrs](0, 1), ilu1.mat.valCrs[icrs * 4 + 1], epsilon);
      EXPECT_NEAR(ilu0.mat.valCrs[icrs](1, 0), ilu1.mat.valCrs[icrs * 4 + 2], epsilon);
      EXPECT_NEAR(ilu0.mat.valCrs[icrs](1, 1), ilu1.mat.valCrs[icrs * 4 + 3], epsilon);
    }
    for(unsigned int iblk=0;iblk<ilu0.mat.valDia.size();++iblk) {
      EXPECT_NEAR(ilu0.mat.valDia[iblk](0, 0), ilu1.mat.valDia[iblk * 4 + 0], epsilon);
      EXPECT_NEAR(ilu0.mat.valDia[iblk](0, 1), ilu1.mat.valDia[iblk * 4 + 1], epsilon);
      EXPECT_NEAR(ilu0.mat.valDia[iblk](1, 0), ilu1.mat.valDia[iblk * 4 + 2], epsilon);
      EXPECT_NEAR(ilu0.mat.valDia[iblk](1, 1), ilu1.mat.valDia[iblk * 4 + 3], epsilon);
    }
  }
  {
    const auto time0 = std::chrono::system_clock::now();
    unsigned int nitr1 = 0;
    for (int itr = 0; itr < 1000; ++itr) { // solve with ILU-CG std::vector
      ilu1.SetValueILU(mA1);
      ilu1.DoILUDecomp();
      std::vector<double> vx1(vb1.size());
      const std::size_t n = vb1.size();
      std::vector<double> tmp0(n), tmp1(n);
      std::vector<double> aConv = Solve_PCG(
          dfm2::CVecXd(vb1), dfm2::CVecXd(vx1), dfm2::CVecXd(tmp0), dfm2::CVecXd(tmp1),
          conv_ratio, iteration, mA1, ilu1);
      nitr1 = aConv.size();
    }
    const auto time1 = std::chrono::system_clock::now();
    unsigned int nitr0 = 0;
    for (int itr = 0; itr < 1000; ++itr) { // solve with ILU-CG Eigen
      delfem2::ILU_CopyValue(ilu0, mA0);
      delfem2::ILU_Decompose(ilu0);
      Eigen::VectorXd vx0(vb0.size());
      const std::size_t n = vb1.size();
      Eigen::VectorXd tmp0(n), tmp1(n);
      std::vector<double> aConv = Solve_PCG(
          vb0, vx0, tmp0, tmp1,
          conv_ratio, iteration, mA0, ilu0);
      nitr0 = aConv.size();
    }
    const auto time2 = std::chrono::system_clock::now();
    double elapsed01 = std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count();
    double elapsed12 = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
    std::cout << "ilu-cg std::vector: " << elapsed01 << "   ilu-cg eigen: " << elapsed12 << std::endl;
    EXPECT_EQ(nitr0,nitr1);
  }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
