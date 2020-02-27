/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <random>
#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/quat.h"
#include "delfem2/voxel.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/specialfuncs.h"
#include "delfem2/funcs.h"
#include "delfem2/evalmathexp.h"
#include "delfem2/primitive.h"
#include "delfem2/slice.h"

#include "delfem2/v23m3q.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(slice,test1){
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CSliceTriMesh> aCS;
  std::vector< std::set<unsigned int> > ReebGraphCS;
  // ----------------------
  delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
           aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<int> aTriSurRel;
  dfm2::ElSuEl_MeshElem(aTriSurRel,
                        aTri.data(), aTri.size()/3, dfm2::MESHELEM_TRI, aXYZ.size()/3);
  // ----------------------
  std::vector<double> aHeight;
  aHeight.push_back(-0.3);
  aHeight.push_back(-0.2);
  aHeight.push_back(-0.1);
  aHeight.push_back(-0.0);
  aHeight.push_back(+0.1);
  aHeight.push_back(+0.2);
  aHeight.push_back(+0.3);
  const double nrm[3] = {0,1,0};
  const double org[3] = {0,0,0};
  std::vector<double> aHeightVtx(aXYZ.size()/3);
  for(unsigned int ip=0;ip<aXYZ.size()/3;++ip){
    double x0 = aXYZ[ip*3+0] - org[0];
    double y0 = aXYZ[ip*3+1] - org[1];
    double z0 = aXYZ[ip*3+2] - org[2];
    aHeightVtx[ip] = nrm[0]*x0 + nrm[1]*y0 + nrm[2]*z0;
  }
  delfem2::Slice_MeshTri3D_Heights(aCS,
                                   aHeight,
                                   aHeightVtx,
                                   aTri,aTriSurRel);
  MakeReebGraph(ReebGraphCS,
                aCS, aTri, aTriSurRel);
  EXPECT_EQ( aCS.size(), ReebGraphCS.size() );
  for(int ics=0;ics<ReebGraphCS.size();++ics){
    for(auto itr = ReebGraphCS[ics].begin();itr!=ReebGraphCS[ics].end();++itr){
      const unsigned int jcs1 = *itr;
      EXPECT_LT( jcs1, aCS.size());
      EXPECT_EQ( abs(aCS[ics].IndHeight() - aCS[jcs1].IndHeight()), 1 );
    }
  }

}


TEST(mathexpeval,test1){
  delfem2::CMathExpressionEvaluator e;
  e.SetKey("x", 3.0);
  e.SetExp("x+3.0");
  EXPECT_EQ(e.Eval(),6);
  e.SetKey("x", 5.0);
  EXPECT_EQ(e.Eval(),8);
}

TEST(funcs,numpy_load_2df){
  std::string path = std::string(PATH_INPUT_DIR)+"/numpy_array4x4_float.npy";
  int ndim0,ndim1;
  std::vector<float> aData;
  bool res = dfm2::LoadNumpy_2DimF(ndim0,ndim1,aData,
                                   path);
  EXPECT_TRUE(res);
  EXPECT_EQ(ndim0,4);
  EXPECT_EQ(ndim1,4);
  EXPECT_EQ(aData.size(),ndim0*ndim1);
}

TEST(funcs,numpy_load_2dd){
  std::string path = std::string(PATH_INPUT_DIR)+"/numpy_array4x4_double.npy";
  int ndim0,ndim1;
  std::vector<double> aData;
  bool res = dfm2::LoadNumpy_2DimD(ndim0,ndim1,aData,
                                   path);
  EXPECT_TRUE(res);
  EXPECT_EQ(ndim0,4);
  EXPECT_EQ(ndim1,4);
  EXPECT_EQ(aData.size(),ndim0*ndim1);
}

TEST(funcs,numpy_load_1df){
  std::string path = std::string(PATH_INPUT_DIR)+"/numpy_array4_float.npy";
  int ndim0;
  std::vector<float> aData;
  bool res = dfm2::LoadNumpy_1DimF(ndim0,aData,
                                   path);
  EXPECT_TRUE(res);
  EXPECT_EQ(ndim0,4);
  EXPECT_EQ(aData.size(),ndim0);
}

TEST(funcs,split_parentheses){
  {
    std::string str = "(a,b),c,(d,e)";
    std::vector<std::string> aS = Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"(a,b)");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"(d,e)");
  }
  {
    std::string str = "(a,b),c";
    std::vector<std::string> aS = Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 2);
    EXPECT_EQ(aS[0],"(a,b)");
    EXPECT_EQ(aS[1],"c");
  }
  {
    std::string str = "a,(b,c)";
    std::vector<std::string> aS = Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 2);
    EXPECT_EQ(aS[0],"a");
    EXPECT_EQ(aS[1],"(b,c)");
  }
}

TEST(funcs,split_quote){
  {
    std::string str = "\"a,b\",c,\"d,e\"";
    std::vector<std::string> aS = Split_Quote(str, ',', '\"' );
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"\"a,b\"");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"\"d,e\"");
  }
}

TEST(mat3, eigen3)
{

  for(int itr=0;itr<10000;itr++){
    double sm[6];
    for(int i=0;i<6;i++){
      sm[i] = ((double)std::rand()/(RAND_MAX+1.0))*100-50;
    }
    double l[3];
    dfm2::CMat3d U;
    dfm2::eigenSym3(U.mat, l,
              sm,20);
    {
      double diffU = (U.Trans()*U-dfm2::CMat3d::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diffU, 0.0, 1.0e-10);
    }
    {
      double L[9] = {l[0],0,0, 0,l[1],0, 0,0,l[2]};
      dfm2::CMat3d UL; dfm2::MatMat3(UL.mat,U.mat,L);
      dfm2::CMat3d ULUt; dfm2::MatMatT3(ULUt.mat, UL.mat, U.mat);
      dfm2::CMat3d SM; SM.SetSymetric(sm);
      double diff = (ULUt-SM).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-6);
    }
  }
  // -----------------------------
  for(int itr=0;itr<100;itr++){
    double sm[6];
    for(int i=0;i<6;i++){
      sm[i] = ((double)std::rand()/(RAND_MAX+1.0))*100-50;
    }
    sm[5] = -sm[4];
    double l[3];
    dfm2::CMat3d U;
    dfm2::eigenSym3(U.mat, l,
              sm,20);
    {
      double diffU = (U.Trans()*U-dfm2::CMat3d::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diffU, 0.0, 1.0e-10);
    }
    {
      double L[9] = {l[0],0,0, 0,l[1],0, 0,0,l[2]};
      dfm2::CMat3d UL; dfm2::MatMat3(UL.mat,U.mat,L);
      dfm2::CMat3d ULUt; dfm2::MatMatT3(ULUt.mat, UL.mat, U.mat);
      dfm2::CMat3d SM; SM.SetSymetric(sm);
      double diff = (ULUt-SM).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-6);
    }
  }
}


TEST(mat3, svd3)
{
  for(int itr=0;itr<10000;itr++){
    dfm2::CMat3d M; M.SetRandom();
    double g[3];
    dfm2::CMat3d U,V;
    dfm2::svd3(U.mat,g,V.mat,
         M.mat,20);
    {
      double diffU = (U.Trans()*U-dfm2::CMat3d::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diffU, 0.0, 1.0e-6);
    }
    {
      double diffV = (V.Trans()*V-dfm2::CMat3d::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diffV, 0.0, 1.0e-10);
    }
    {
      const double G[9] = {g[0],0,0, 0,g[1],0, 0,0,g[2]};
      dfm2::CMat3d UG;   dfm2::MatMat3(UG.mat, U.mat,G);
      dfm2::CMat3d UGVt; dfm2::MatMatT3(UGVt.mat, UG.mat,V.mat);
      double diff = (UGVt - M).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-10);
    }
  }
}


TEST(mat3, rot_comp)
{
  for(int itr=0;itr<10000;itr++){
    dfm2::CMat3d M; M.SetRandom();
    dfm2::CMat3d R; dfm2::GetRotPolarDecomp(R.mat, M.mat, 40);
    {
      double diff = (R.Trans()*R-dfm2::CMat3d::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-5);
    }
    {
      dfm2::CMat3d MR = M.MatMat(R.Trans());
      double diff0 = (MR-MR.Sym()).SqNorm_Frobenius();
      EXPECT_NEAR(diff0, 0.0, 1.0e-5);
    }
    {
      dfm2::CMat3d RM = (R.Trans()).MatMat(M);
      double diff1 = (RM-RM.Sym()).SqNorm_Frobenius();
      EXPECT_NEAR(diff1, 0.0, 1.0e-5);
    }
  }
}

TEST(mat3, quat)
{
  std::uniform_real_distribution<double> dist(-50.0, +50.0);
  std::mt19937 mtd;
  for(int itr=0;itr<10000;itr++){
    double quat[4] = { dist(mtd), dist(mtd), dist(mtd), dist(mtd) };
    dfm2::Normalize_Quat(quat);
    dfm2::CMat3d R;
    R.SetRotMatrix_Quaternion(quat);
    {
      double diff = (R.Trans()*R-dfm2::CMat3d::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-14);
    }
    double puat[4];
    R.GetQuat_RotMatrix(puat);
    dfm2::CMat3d P;
    P.SetRotMatrix_Quaternion(puat);
    double diff = (P-R).SqNorm_Frobenius();
    EXPECT_NEAR(diff, 0.0, 1.0e-20);
  }
  
  for(int itr=0;itr<10000;itr++){
    dfm2::CQuat<double> q0(dist(mtd),dist(mtd),dist(mtd),dist(mtd) );
    dfm2::CQuat<double> q1(dist(mtd),dist(mtd),dist(mtd),dist(mtd) );
    dfm2::CQuat<double> q2 = q0 + q1;
  }
}

TEST(mshio,load_obj)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Obj(std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
                    aXYZ, aTri);
  EXPECT_EQ(aTri.size(),1000*3);
}

TEST(vec2,second_moment_of_area)
{
  for(int itr=0;itr<100;itr++){
    double r0 = (double)rand()/(RAND_MAX+1.0);
    double r1 = (double)rand()/(RAND_MAX+1.0);
    double r2 = (double)rand()/(RAND_MAX+1.0);
    double r3 = (double)rand()/(RAND_MAX+1.0);
    double r4 = (double)rand()/(RAND_MAX+1.0);
    double a = 10*r0;
    double b = a*(3*r1+1);
    std::vector<dfm2::CVec2d> aVec2;
    {
      aVec2.push_back( dfm2::CVec2d(-a*0.5,-b*0.5) );
      aVec2.push_back( dfm2::CVec2d(+a*0.5,-b*0.5) );
      aVec2.push_back( dfm2::CVec2d(+a*0.5,+b*0.5) );
      aVec2.push_back( dfm2::CVec2d(-a*0.5,+b*0.5) );
      double theta0 = r4*3.1415*2.0;
      dfm2::Rotate(aVec2,theta0);
      dfm2::Translate(aVec2, r2*10-5,r3*10-5);
    }
    dfm2::CVec2d cg,pa1,pa2;
    double area,I1,I2;
    dfm2::SecondMomentOfArea_Polygon(cg,area, pa1,I1, pa2,I2,
                               aVec2);
    EXPECT_NEAR(area, a*b, 1.0e-10);
    EXPECT_NEAR(pa1*pa2, 0.0, 1.0e-10 );
    EXPECT_TRUE(I1>=I2);
//    EXPECT_NEAR(pa1.x*pa1.x,  1.0,          1.0e-10);
//    EXPECT_NEAR(pa1.y,        0.0,          1.0e-10);
    EXPECT_NEAR(I1,           a*b*b*b/12.0, 1.0e-10 );
    ///
//    EXPECT_NEAR(pa2.x,        0.0,          1.0e-10);
//    EXPECT_NEAR(pa2.y*pa2.y,  1.0,          1.0e-10);
    EXPECT_NEAR(I2,           b*a*a*a/12.0, 1.0e-10 );
  }
}

TEST(meshtopo,quad_subdiv0)
{
  dfm2::CVoxelGrid3D vg;
  vg.Add(0,0,0);
  {
    std::vector<double> aXYZ0;
    std::vector<unsigned int> aQuad0;
    vg.GetQuad(aXYZ0, aQuad0);
    EXPECT_EQ(aXYZ0.size(),8*3);
    EXPECT_EQ(aQuad0.size(),6*4);
    {
      std::vector<unsigned int> aTri0;
      delfem2::convert2Tri_Quad(aTri0, aQuad0);
      EXPECT_EQ(aTri0.size(),12*3);
    }
  }
}

TEST(meshtopo,quad_subdiv1)
{
  dfm2::CVoxelGrid3D vg;
  vg.Add(0,0,0);
  vg.Add(1,0,0);
  vg.Add(1,1,0);
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aQuad0;
  vg.GetQuad(aXYZ0, aQuad0);
  EXPECT_EQ(aXYZ0.size(),18*3);
  EXPECT_EQ(aQuad0.size(),14*4);
  std::vector<double> aXYZ0a;
  std::vector<unsigned int> aQuad0a;
  std::vector<int> mapInOut;
  dfm2::RemoveUnreferencedPoints_MeshElem(aXYZ0a,aQuad0a, mapInOut,
                                          3,aXYZ0,aQuad0);

  EXPECT_EQ(aXYZ0a.size(),16*3);
  EXPECT_EQ(aQuad0a.size(),14*4);
}


TEST(mathfunc,sherical_harmonics_orthgonality)
{
  std::vector<double> aXYZ,aVal;
  std::vector<unsigned int> aTri;
  delfem2::MeshTri3D_Cube(aXYZ,aTri, 50);
  for(int ip=0;ip<aXYZ.size()/3;ip++){
    double x = aXYZ[ip*3+0];
    double y = aXYZ[ip*3+1];
    double z = aXYZ[ip*3+2];
    double invlen = 1.0/sqrt(x*x+y*y+z*z);
    aXYZ[ip*3+0] *= invlen;
    aXYZ[ip*3+1] *= invlen;
    aXYZ[ip*3+2] *= invlen;
  }
  const int norder = 9;
  const int N = (norder+1)*(norder+1);
  double area_sum = 0.0;
  std::vector<double> A(N*N,0.0);
  for(int it=0;it<aTri.size()/3;++it){
    const int i0 = aTri[it*3+0];
    const int i1 = aTri[it*3+1];
    const int i2 = aTri[it*3+2];
    dfm2::CVec3d p0(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
    dfm2::CVec3d p1(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    dfm2::CVec3d p2(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    double area = SolidAngleTri(p0, p1, p2);
    area_sum += area;
    double a0[N]; makeArray_SphericalHarmonics(a0, norder, p0.p[0],p0.p[1],p0.p[2]);
    double a1[N]; makeArray_SphericalHarmonics(a1, norder, p1.p[0],p1.p[1],p1.p[2]);
    double a2[N]; makeArray_SphericalHarmonics(a2, norder, p2.p[0],p2.p[1],p2.p[2]);
    for(int ish=0;ish<N;++ish){
      for(int jsh=0;jsh<N;++jsh){
        double val = 2*a0[ish]*a0[jsh] + 2*a1[ish]*a1[jsh] + 2*a2[ish]*a2[jsh];
        val += a0[ish]*a1[jsh] + a0[ish]*a2[jsh];
        val += a1[ish]*a0[jsh] + a1[ish]*a2[jsh];
        val += a2[ish]*a0[jsh] + a2[ish]*a1[jsh];
        val *= area/12.0;
        A[N*ish+jsh] += val;
      }
    }
  }
  EXPECT_NEAR(area_sum, M_PI*4, 1.0e-5);
  for(int ish=0;ish<N;++ish){
    for(int jsh=0;jsh<N;++jsh){
      if( ish == jsh ){ continue; }
      EXPECT_NEAR(A[N*ish+jsh], 0.0, 3.0e-3);
    }
  }
  for(int iorder=0;iorder<norder;++iorder){
    for(int ish=iorder*iorder;ish<(iorder+1)*(iorder+1);++ish){
      if( ish == iorder*(iorder+1) ){
        EXPECT_NEAR(A[N*ish+ish], 1.0, 1.5e-2);
      }
      else{
        EXPECT_NEAR(A[N*ish+ish], 0.5, 1.5e-2);
      }
    }
  }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
