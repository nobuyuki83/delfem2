/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>

#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/emat.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"

#include "delfem2/v23m3q.h"
#include "delfem2/objfunc_v23.h"
#include "delfem2/dtri_v2.h"
#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"

namespace dfm2 = delfem2;

// --------------------------------------

TEST(objfunc_v23, Check_CdC_TriStrain){
  for(int itr=0;itr<200;++itr){
    const double P[3][2] = {
      { 10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5) },
      { 10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5) },
      { 10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5) },
    };
    double a0 = dfm2::Area_Tri2(P[0], P[1], P[2]);
    if( fabs(a0) < 0.1 ) continue;
    const double p[3][3] = {
      { 10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5) },
      { 10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5) },
      { 10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5),
        10.0*(rand()/(RAND_MAX+1.0)-0.5) },
    };
    double diff = dfm2::Check_CdC_TriStrain(P, p, 1.0e-5);
    EXPECT_LT(diff, 0.2);
  }
}


TEST(objfunc_v23, MIPS)
{
  double C[3][3];
  for(auto & ino : C){
    ino[0] = (double)rand()/(RAND_MAX+1.0);
    ino[1] = (double)rand()/(RAND_MAX+1.0);
    ino[2] = (double)rand()/(RAND_MAX+1.0);
  }
  double c[3][3];
  dfm2::CMat3d m;
  m.SetRotMatrix_Cartesian(0.3, 1.0, 0.5);
  for(int ino=0;ino<3;++ino){
    m.MatVec(C[ino], c[ino]);
  }
  double E, dE[3][3], ddE[3][3][3][3];
  dfm2::Energy_MIPS(E, dE, ddE,
                    c, C);
  for(int ino=0;ino<3;++ino){
    for(int idim=0;idim<3;++idim){
      double c1[3][3] = {
        {c[0][0],c[0][1],c[0][2]},
        {c[1][0],c[1][1],c[1][2]},
        {c[2][0],c[2][1],c[2][2]} };
      double eps = 1.0e-5;
      c1[ino][idim] += eps;
      double E1, dE1[3][3], ddE1[3][3][3][3];
      dfm2::Energy_MIPS(E1, dE1, ddE1,
                        c1, C);
      EXPECT_NEAR( (E1-E)/eps, dE[ino][idim], 1.0e-3);
      for(int jno=0;jno<3;++jno){
        for(int jdim=0;jdim<3;++jdim){
          EXPECT_NEAR( (dE1[jno][jdim]-dE[jno][jdim])/eps, ddE[jno][ino][jdim][idim], 1.e-3);
        }
      }
    }
  }
}

TEST(objfunc_v23, distancetri2d3d)
{
  double P[3][2];  // undeformed triangle vertex positions
  for(auto & P0 : P){
    P0[0] = (double)rand()/(RAND_MAX+1.0);
    P0[1] = (double)rand()/(RAND_MAX+1.0);
  }
  double p[3][3]; //  deformed triangle vertex positions)
  for(auto & p0 : p){
    p0[0] = (double)rand()/(RAND_MAX+1.0);
    p0[1] = (double)rand()/(RAND_MAX+1.0);
    p0[2] = (double)rand()/(RAND_MAX+1.0);
  }
   // --------------
  double C[3], dCdp[3][9];
  dfm2::PBD_ConstraintProjection_DistanceTri2D3D(C, dCdp, P, p);
  for(int ine=0;ine<3;++ine){
    for(int idim=0;idim<3;++idim){
      double eps = 1.0e-6;
      double p1[3][3]; for(int i=0;i<9;++i){ (&p1[0][0])[i] = (&p[0][0])[i]; }
      p1[ine][idim] += eps;
      double C1[3], dCdp1[3][9];
      dfm2::PBD_ConstraintProjection_DistanceTri2D3D(C1, dCdp1, P, p1);
      EXPECT_NEAR( (C1[0]-C[0])/eps, dCdp[0][ine*3+idim], 1.0e-5 );
      EXPECT_NEAR( (C1[1]-C[1])/eps, dCdp[1][ine*3+idim], 1.0e-5 );
      EXPECT_NEAR( (C1[2]-C[2])/eps, dCdp[2][ine*3+idim], 1.0e-5 );
    }
  }
}

// ------------------------------------------------------------

TEST(fem,plate_bending_mitc3_emat)
{
  for(int itr=0;itr<200;++itr){
    double C[3][2];
    for(int i=0;i<6;++i){
      (&C[0][0])[i] = 10.0*(rand()/(RAND_MAX+1.0)-0.5);
    }
    double a0 = dfm2::Area_Tri2(C[0], C[1], C[2]);
    if( a0 < 0.1 ) continue;
    double u[3][3];
    for(int i=0;i<9;++i){
      (&u[0][0])[i] = 1.0*(rand()/(RAND_MAX+1.0)-0.5);
    }
    double thickness = (rand()+1.0)/(RAND_MAX+1.0);
    double lambda = (rand()+1.0)/(RAND_MAX+1.0);
    double myu = (rand()+1.0)/(RAND_MAX+1.0);
    double diff = dfm2::Check_WdWddW_PlateBendingMITC3(C, u,
                                                       thickness,lambda,myu, 1.0e-6);
    EXPECT_LT(diff,2.0e-3);
  }
}

TEST(fem,plate_bending_mitc3_cantilever)
{
  const double lambda = 0.0;
  const double lenx0 = 1.0;
  const double leny0 = 0.2;
  const double thickness0 = 0.05;
  const double myu0 = 10000.0;
  const double rho0 = 1.0;
  const double gravity_z0 = -10.0;
  const double elen0 = 0.03;
  for(int itr = 0;itr<10;++itr ){
    const double lenx = lenx0*(1.0+rand()/(RAND_MAX+1.0));
    const double leny = leny0*(1.0+rand()/(RAND_MAX+1.0));
    const double thickness = thickness0*(1.0+rand()/(RAND_MAX+1.0));
    const double myu = myu0*(1.0+rand()/(RAND_MAX+1.0));;
    const double rho = rho0*(1.0+rand()/(RAND_MAX+1.0));
    const double gravity_z = gravity_z0*(1.0+rand()/(RAND_MAX+1.0));
    const double elen = elen0*(1.0+rand()/(RAND_MAX+1.0));
    std::vector<unsigned int> aTri;
    std::vector<double> aXY0;
    {
      std::vector< std::vector<double> > aaXY;
      {
        aaXY.resize(1);
        aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(-leny*0.5);
        aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(-leny*0.5);
        aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(+leny*0.5);
        aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(+leny*0.5);
      }
      // ---------------------
      std::vector<dfm2::CDynPntSur> aPo2D;
      std::vector<dfm2::CDynTri> aETri;
      std::vector<dfm2::CVec2d> aVec2;
      GenMesh(aPo2D, aETri, aVec2,
              aaXY, elen, elen);
      MeshTri2D_Export(aXY0,aTri,
                       aVec2,aETri);
    }
    std::vector<int> aBCFlag; // boundary condition flag
    {
      const int np = (int)aXY0.size()/2;
      aBCFlag.assign(np*3, 0);
      for(int ip=0;ip<np;++ip){
        const double px = aXY0[ip*2+0];
        //    const double py = aXY0[ip*2+1];
        if( fabs(px-(-lenx*0.5)) < 0.0001 ){
          aBCFlag[ip*3+0] = 1;
          aBCFlag[ip*3+1] = 1;
          aBCFlag[ip*3+2] = 1;
        }
      }
    }
    dfm2::CMatrixSparse<double> mat_A;
    dfm2::CPreconditionerILU<double> ilu_A;
    {
      std::vector<unsigned int> psup_ind, psup;
      dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                                 aTri.data(), aTri.size()/3, 3,
                                 (int)aXY0.size()/2);
      dfm2::JArray_Sort(psup_ind, psup);
      ////
      const int np = (int)aXY0.size()/2;
      mat_A.Initialize(np, 3, true);
      mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
      ilu_A.Initialize_ILU0(mat_A);
    }
    std::vector<double> aVal;
    aVal.assign(aXY0.size()/2*3, 0.0);
    std::vector<double> vec_b;
    {
      const int np = (int)aXY0.size()/2;
      const int nDoF = np*3;
      // -------------------
      mat_A.SetZero();
      vec_b.assign(nDoF, 0.0);
      dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(mat_A,vec_b.data(),
                                                               thickness,lambda,myu,
                                                               rho,gravity_z,
                                                               aXY0.data(), aXY0.size()/2,
                                                               aTri.data(), aTri.size()/3,
                                                               aVal.data());
      mat_A.SetFixedBC(aBCFlag.data());
      dfm2::setRHS_Zero(vec_b, aBCFlag,0);
      // --------------------------
      std::vector<double> vec_x;
      {
        ilu_A.SetValueILU(mat_A);
        ilu_A.DoILUDecomp();
        vec_x.resize(vec_b.size());
        std::vector<double> conv = Solve_PCG(vec_b.data(), vec_x.data(),
                                             vec_b.size(),
                                             1.0e-5, 1000,
                                             mat_A, ilu_A);
//        std::cout << "convergence   nitr:" << conv.size() << "    res:" << conv[conv.size()-1] << std::endl;
        EXPECT_LT( conv.size(), 1000 );
        EXPECT_LT( conv[conv.size()-1], 1.0e-5);
      }
      // -------------------------
      dfm2::XPlusAY(aVal,nDoF,aBCFlag,
                    1.0,vec_x);
    }
    {
      assert( fabs(lambda)<1.0e-10 );
      const double E = myu*2.0;
      const double I = thickness*thickness*thickness*leny/12.0;
      const double W = thickness*lenx*leny*rho*gravity_z;
      const double w = W/lenx;
      const double disp = w*(lenx*lenx*lenx*lenx)/(8.0*E*I);
//      std::cout << "disp:" << disp << std::endl;
      for(int ip=0;ip<aXY0.size()/2;++ip){
        const double px = aXY0[ip*2+0];
        if( fabs(px-(+lenx*0.5)) > 0.0001 ){ continue; }
        EXPECT_LE( fabs(aVal[ip*3+0] - disp), 0.002*fabs(disp) );
//        std::cout << aVal[ip*3+0] << " " << disp << "  " << fabs(aVal[ip*3+0] - disp)/disp << std::endl;
      }
    }
  }
}
