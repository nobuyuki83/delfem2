/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Core>

#include "gtest/gtest.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/lsmats.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/mshuni.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/dtri.h"
#include "delfem2/eigen/ls_dense.h"
#include "delfem2/eigen/ls_sparse.h"
#include "delfem2/eigen/ls_ilu_sparse.h"
#include "delfem2/lsitrsol.h"


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
  delfem2::CMatrixSparseBlock<Eigen::Matrix2d,Eigen::aligned_allocator<Eigen::Matrix2d>> Aeig;
  Eigen::VectorXd Veig0(nDoF);
  Eigen::Matrix<double,-1,2,Eigen::RowMajor> Veig1(np,2);
  dfm2::CMatrixSparse<double> Astd;
  std::vector<double> Vstd(nDoF);
  {
    std::vector<unsigned int> psup_ind0, psup0;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind0, psup0,
        aTri1.data(), aTri1.size() / 3, 3,
        aXY1.size() / 2);
    Aeig.Initialize(np);
    Aeig.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
    const double myu = 10.0, lambda = 10.0, rho = 1.0, g_x = 0.0, g_y = -3.0;
    Aeig.setZero();
    Veig0.setZero();
    dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
        Aeig, Veig0.data(),
        myu, lambda, rho, g_x, g_y,
        aXY1.data(), aXY1.size() / 2,
        aTri1.data(), aTri1.size() / 3,
        aVal.data());
    SetFixedBC_Dia(Aeig, aBCFlag.data(), 1.f);
    SetFixedBC_Col(Aeig, aBCFlag.data());
    SetFixedBC_Row(Aeig, aBCFlag.data());
    delfem2::setZero_Flag(Veig0, aBCFlag, 0);
    // ----------
    for(unsigned int ip=0;ip<np;++ip){
      Veig1(ip,0) = Veig0(ip*2+0);
      Veig1(ip,1) = Veig0(ip*2+1);
    }
    // ----------
    Astd.Initialize(np, 2, true);
    Astd.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
    Astd.setZero();
    Vstd.assign(nDoF,0.0);
    dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
        Astd, Vstd.data(),
        myu, lambda, rho, g_x, g_y,
        aXY1.data(), aXY1.size() / 2,
        aTri1.data(), aTri1.size() / 3,
        aVal.data());
    Astd.SetFixedBC(aBCFlag.data());
    dfm2::setRHS_Zero(Vstd, aBCFlag, 0);
  }
  // ---------------
  double conv_ratio = 1.0e-6;
  int iteration = 10000;
  {
    const unsigned int max_itr = 10;
    const auto time0 = std::chrono::system_clock::now();
    unsigned int nitr1 = 0;
    for(int itr=0;itr<max_itr;++itr){ // CG method std::vector
      const std::size_t n = Vstd.size();
      std::vector<double> tmp0(n), tmp1(n), vx1(n), vB1(n);
      vB1 = Vstd;
      std::vector<double> aConv = Solve_CG(
          dfm2::ViewAsVectorXd(vB1),
          dfm2::ViewAsVectorXd(vx1),
          dfm2::ViewAsVectorXd(tmp0),
          dfm2::ViewAsVectorXd(tmp1),
          conv_ratio, iteration, Astd);
      nitr1 = aConv.size();
    }
    const auto time1 = std::chrono::system_clock::now();
    // ---------------
    unsigned int nitr0 = 0;
    for(int itr=0;itr<max_itr;++itr){ // CG method with eigen
      const std::size_t n = Veig0.size();
      Eigen::VectorXd tmp0(n), tmp1(n), vx0(n), vB0(n);
      vB0 = Veig0;
      std::vector<double> aConv = delfem2::Solve_CG(
          vB0, vx0, tmp0, tmp1,
          conv_ratio, iteration, Aeig);
      nitr0 = aConv.size();
    }
    const auto time2 = std::chrono::system_clock::now();
    // ---------------
    unsigned int nitr2 = 0;
    for(int itr=0;itr<max_itr;++itr){ // CG method with eigen
      Eigen::Matrix<double,-1,2,Eigen::RowMajor> tmp0(np,2), tmp1(np,2), vx0(np,2), vB0(np,2);
      vB0 = Veig1;
      std::vector<double> aConv = delfem2::Solve_CG(
          vB0, vx0, tmp0, tmp1,
          conv_ratio, iteration, Aeig);
      nitr2 = aConv.size();
    }
    const auto time3 = std::chrono::system_clock::now();
    EXPECT_NEAR(nitr0,nitr1,10);
    EXPECT_NEAR(nitr1,nitr2,10);
    double elapsed01 = std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count();
    double elapsed12 = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
    double elapsed23 = std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count();
    std::cout << "cg std::vector: " << elapsed01 << " " << nitr0 << std::endl;
    std::cout << "cg eigen vec: " << elapsed12 << " " << nitr1 << std::endl;
    std::cout << "cg eigen mat: " << elapsed23 << " " << nitr2 << std::endl;
  }
  delfem2::CILU_SparseBlock<Eigen::Matrix2d,Eigen::aligned_allocator<Eigen::Matrix2d>> ilu0;
  delfem2::CPreconditionerILU<double> ilu1;
  { // LU
    ilu1.SetPattern0(Astd);
    ilu1.CopyValue(Astd);
    ilu1.Decompose();
    delfem2::ILU_SetPattern0(ilu0, Aeig);
    delfem2::ILU_CopyValue(ilu0, Aeig);
    delfem2::ILU_Decompose(ilu0);
    // check if the entry is the same
    for(unsigned int icrs=0;icrs<ilu0.valCrs.size();++icrs) {
      EXPECT_NEAR(ilu0.valCrs[icrs](0, 0), ilu1.valCrs[icrs * 4 + 0], epsilon);
      EXPECT_NEAR(ilu0.valCrs[icrs](0, 1), ilu1.valCrs[icrs * 4 + 1], epsilon);
      EXPECT_NEAR(ilu0.valCrs[icrs](1, 0), ilu1.valCrs[icrs * 4 + 2], epsilon);
      EXPECT_NEAR(ilu0.valCrs[icrs](1, 1), ilu1.valCrs[icrs * 4 + 3], epsilon);
    }
    for(unsigned int iblk=0;iblk<ilu0.valDia.size();++iblk) {
      EXPECT_NEAR(ilu0.valDia[iblk](0, 0), ilu1.valDia[iblk * 4 + 0], epsilon);
      EXPECT_NEAR(ilu0.valDia[iblk](0, 1), ilu1.valDia[iblk * 4 + 1], epsilon);
      EXPECT_NEAR(ilu0.valDia[iblk](1, 0), ilu1.valDia[iblk * 4 + 2], epsilon);
      EXPECT_NEAR(ilu0.valDia[iblk](1, 1), ilu1.valDia[iblk * 4 + 3], epsilon);
    }
  }
  {
    const auto time0 = std::chrono::system_clock::now();
    const unsigned int max_itr = 10;
    unsigned int nitr1 = 0;
    for (int itr = 0; itr < max_itr; ++itr) { // solve with ILU-CG std::vector
      ilu1.CopyValue(Astd);
      ilu1.Decompose();
      const std::size_t n = Vstd.size();
      std::vector<double> tmp0(n), tmp1(n), vx1(n), vB1(n);
      vB1 = Vstd;
      std::vector<double> aConv = Solve_PCG(
          dfm2::ViewAsVectorXd(vB1),
          dfm2::ViewAsVectorXd(vx1),
          dfm2::ViewAsVectorXd(tmp0),
          dfm2::ViewAsVectorXd(tmp1),
          conv_ratio, iteration, Astd, ilu1);
      nitr1 = aConv.size();
    }
    const auto time1 = std::chrono::system_clock::now();
    unsigned int nitr0 = 0;
    for (int itr = 0; itr < max_itr; ++itr) { // solve with ILU-CG Eigen
      delfem2::ILU_CopyValue(ilu0, Aeig);
      delfem2::ILU_Decompose(ilu0);
      const std::size_t n = Vstd.size();
      Eigen::VectorXd tmp0(n), tmp1(n), vx0(n), vB0(n);
      vB0 = Veig0;
      std::vector<double> aConv = Solve_PCG(
          vB0, vx0, tmp0, tmp1,
          conv_ratio, iteration, Aeig, ilu0);
      nitr0 = aConv.size();
    }
    const auto time2 = std::chrono::system_clock::now();
    unsigned int nitr2 = 0;
    for (int itr = 0; itr < max_itr; ++itr) { // solve with ILU-CG Eigen
      delfem2::ILU_CopyValue(ilu0, Aeig);
      delfem2::ILU_Decompose(ilu0);
      Eigen::Matrix<double,-1,2,Eigen::RowMajor> tmp0(np,2), tmp1(np,2), vx0(np,2), vB0(np,2);
      vB0 = Veig1;
      std::vector<double> aConv = Solve_PCG(
          vB0, vx0, tmp0, tmp1,
          conv_ratio, iteration, Aeig, ilu0);
      nitr2 = aConv.size();
    }
    const auto time3 = std::chrono::system_clock::now();
    EXPECT_NEAR(nitr0,nitr1,10);
    EXPECT_NEAR(nitr1,nitr2,10);
    double elapsed01 = std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count();
    double elapsed12 = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
    double elapsed23 = std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count();
    std::cout << "ilu-cg std::vector: " << elapsed01 << " " << nitr0 << std::endl;
    std::cout << "ilu-cg eigen vec: " << elapsed12 << " " << nitr1 << std::endl;
    std::cout << "ilu-cg eigen mat: " << elapsed23 << " " << nitr2 << std::endl;
  }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
