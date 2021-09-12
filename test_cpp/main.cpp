/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"
//
#include "delfem2/thread/th.h"
#include "delfem2/filenpy_str.h"
#include "delfem2/str.h"
#include "delfem2/vec3.h"
#include "delfem2/geosolidelm_v3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/specialfuncs.h"
#include "delfem2/evalmathexp.h"
#include <cstring>
#include <random>

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------


TEST(mathexpeval,test1){
  delfem2::CMathExpressionEvaluator e;
  e.SetKey("x", 3.0);
  e.SetExp("x+3.0");
  EXPECT_DOUBLE_EQ(e.Eval(),6);
  e.SetKey("x", 5.0);
  EXPECT_DOUBLE_EQ(e.Eval(),8.0);
  //
  e.SetKey("x", 1.0);
  e.SetKey("y", 2.0);
  e.SetExp("x+y");
  EXPECT_DOUBLE_EQ(e.Eval(),3.0);
  //
  e.SetExp("sin(PI*0.5*x)");
  EXPECT_DOUBLE_EQ(e.Eval(),1.0);
}

TEST(funcs,numpy_load_2df){
  std::string path = std::string(PATH_INPUT_DIR)+"/numpy_array4x4_float.npy";
  int ndim0,ndim1;
  std::vector<float> aData;
  bool res = dfm2::LoadNumpy_2Dim(
      ndim0,ndim1,aData,
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
  bool res = dfm2::LoadNumpy_2Dim(
      ndim0,ndim1,aData,
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
  bool res = dfm2::LoadNumpy_1DimF(
      ndim0,aData,
      path);
  EXPECT_TRUE(res);
  EXPECT_EQ(ndim0,4);
  EXPECT_EQ(aData.size(),ndim0);
}

TEST(funcs,split_parentheses){
  {
    std::string str = "(a,b),c,(d,e)";
    std::vector<std::string> aS = dfm2::Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"(a,b)");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"(d,e)");
  }
  {
    std::string str = "(a,b),c";
    std::vector<std::string> aS = dfm2::Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 2);
    EXPECT_EQ(aS[0],"(a,b)");
    EXPECT_EQ(aS[1],"c");
  }
  {
    std::string str = "a,(b,c)";
    std::vector<std::string> aS = dfm2::Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 2);
    EXPECT_EQ(aS[0],"a");
    EXPECT_EQ(aS[1],"(b,c)");
  }
}

TEST(funcs,split_quote){
  {
    std::string str = R"("a,b",c,"d,e")";
    std::vector<std::string> aS = dfm2::Split_Quote(str, ',', '\"' );
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"\"a,b\"");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"\"d,e\"");
  }
  {
    std::string str = R"("a,b",,c,"d,e")";
    std::vector<std::string> aS = dfm2::Split_Quote(str, ',', '\"' );
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"\"a,b\"");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"\"d,e\"");
  }
}

TEST(funcs,split){
  std::vector<std::string> aToken;
  aToken = dfm2::Split("chr count=80"," =");
  EXPECT_EQ(aToken.size(), 3);
  EXPECT_EQ(aToken[0],"chr");
  EXPECT_EQ(aToken[1],"count");
  EXPECT_EQ(aToken[2],"80");
  //
  aToken = dfm2::Split("chr = 80"," =");
  EXPECT_EQ(aToken.size(), 2);
  EXPECT_EQ(aToken[0],"chr");
  EXPECT_EQ(aToken[1],"80");
  //
  aToken = dfm2::Split("=chr = 80="," =");
  EXPECT_EQ(aToken.size(), 2);
  EXPECT_EQ(aToken[0],"chr");
  EXPECT_EQ(aToken[1],"80");
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
    const unsigned int i0 = aTri[it*3+0];
    const unsigned int i1 = aTri[it*3+1];
    const unsigned int i2 = aTri[it*3+2];
    dfm2::CVec3d p0(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
    dfm2::CVec3d p1(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    dfm2::CVec3d p2(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    double area = SolidAngleTri(p0, p1, p2);
    area_sum += area;
    double a0[N]; dfm2::makeArray_SphericalHarmonics(a0, norder, p0.p[0],p0.p[1],p0.p[2]);
    double a1[N]; dfm2::makeArray_SphericalHarmonics(a1, norder, p1.p[0],p1.p[1],p1.p[2]);
    double a2[N]; dfm2::makeArray_SphericalHarmonics(a2, norder, p2.p[0],p2.p[1],p2.p[2]);
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

TEST(thread,parallel_for0)
{
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_int_distribution<unsigned int> dist0(0,100);
  std::uniform_int_distribution<unsigned int> dist1(0,5);
  for(unsigned int itr=0;itr<100;++itr) {
    const unsigned int N = dist0(rdeng);
    std::vector<int> aIn(N); for(unsigned int i=0;i<N;++i){ aIn[i] = i; }
    std::vector<int> aOut(aIn.size());
    auto func0 = [&aIn, &aOut](int i) { aOut[i] = aIn[i] * aIn[i]; }; // square
    const unsigned int nthread = dist1(rdeng);
    dfm2::thread::parallel_for(
        static_cast<unsigned int>(aIn.size()),
        func0,
        nthread);
    std::vector<int> aTrg(N); for(unsigned int i=0;i<N;++i){ aTrg[i] = i*i; }
    EXPECT_TRUE(0 == std::memcmp(aOut.data(), aTrg.data(), aIn.size() * sizeof(int)));
  }
}


int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
