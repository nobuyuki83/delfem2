#include <iostream>

#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/quat.h"
#include "delfem2/voxel.h"
#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/mshio.h"
#include "delfem2/mathfuncs.h"

TEST(mat3, eigen3)
{
  for(int itr=0;itr<10000;itr++){
    double sm[6];
    for(int i=0;i<6;i++){
      sm[i] = ((double)std::rand()/(RAND_MAX+1.0))*100-50;
    }
    double l[3];
    CMatrix3 U;
    eigenSym3(U.mat, l,
              sm,20);
    {
      double diffU = (U.Trans()*U-CMatrix3::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diffU, 0.0, 1.0e-10);
    }
    {
      double L[9] = {l[0],0,0, 0,l[1],0, 0,0,l[2]};
      CMatrix3 UL; MatMat3(UL.mat,U.mat,L);
      CMatrix3 ULUt; MatMatTrans3(ULUt.mat, UL.mat, U.mat);
      CMatrix3 SM; SM.SetSymetric(sm);
      double diff = (ULUt-SM).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-6);
    }
  }
}


TEST(mat3, svd3)
{
  for(int itr=0;itr<10000;itr++){
    CMatrix3 M; M.SetRandom();
    double g[3];
    CMatrix3 U,V;
    svd3(U.mat,g,V.mat,
         M.mat,20);
    {
      double diffU = (U.Trans()*U-CMatrix3::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diffU, 0.0, 1.0e-10);
    }
    {
      double diffV = (V.Trans()*V-CMatrix3::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diffV, 0.0, 1.0e-10);
    }
    {
      const double G[9] = {g[0],0,0, 0,g[1],0, 0,0,g[2]};
      CMatrix3 UG;   MatMat3(UG.mat, U.mat,G);
      CMatrix3 UGVt; MatMatTrans3(UGVt.mat, UG.mat,V.mat);
      double diff = (UGVt - M).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-10);
    }
  }
}


TEST(mat3, rot_comp)
{
  for(int itr=0;itr<10000;itr++){
    CMatrix3 M; M.SetRandom();
    CMatrix3 R; GetRotPolarDecomp(R.mat, M.mat, 40);
    {
      double diff = (R.Trans()*R-CMatrix3::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-14);
    }
    {
      CMatrix3 MR = M.MatMat(R.Trans());
      double diff0 = (MR-MR.Sym()).SqNorm_Frobenius();
      EXPECT_NEAR(diff0, 0.0, 1.0e-14);
    }
    {
      CMatrix3 RM = (R.Trans()).MatMat(M);
      double diff1 = (RM-RM.Sym()).SqNorm_Frobenius();
      EXPECT_NEAR(diff1, 0.0, 1.0e-14);
    }
  }
}

TEST(mat3, quat)
{
  for(int itr=0;itr<10000;itr++){
    double quat[4];
    quat[0] = (rand()/(RAND_MAX+1.0))*100-50;
    quat[1] = (rand()/(RAND_MAX+1.0))*100-50;
    quat[2] = (rand()/(RAND_MAX+1.0))*100-50;
    quat[3] = (rand()/(RAND_MAX+1.0))*100-50;
    QuatNormalize(quat);
    CMatrix3 R;
    R.SetRotMatrix_Quaternion(quat);
    {
      double diff = (R.Trans()*R-CMatrix3::Identity()).SqNorm_Frobenius();
      EXPECT_NEAR(diff, 0.0, 1.0e-14);
    }
    double puat[4];
    R.GetQuat_RotMatrix(puat);
    CMatrix3 P;
    P.SetRotMatrix_Quaternion(puat);
    double diff = (P-R).SqNorm_Frobenius();
    EXPECT_NEAR(diff, 0.0, 1.0e-20);
  }
}

TEST(mshio,load_obj)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  Read_Obj(std::string(PATH_INPUT_DIR)+"/bunny_1k.obj", aXYZ, aTri);
  EXPECT_EQ(aTri.size(),1000*3);
}

TEST(vec2,second_moment_of_area)
{
  for(int itr=0;itr<10;itr++){
    double r0 = (double)rand()/(RAND_MAX+1.0);
    double r1 = (double)rand()/(RAND_MAX+1.0);
    double r2 = (double)rand()/(RAND_MAX+1.0);
    double r3 = (double)rand()/(RAND_MAX+1.0);
    double r4 = (double)rand()/(RAND_MAX+1.0);
    double a = 10*r0;
    double b = a*(3*r1+1);
    std::vector<CVector2> aVec2;
    {
      aVec2.push_back( CVector2(-a*0.5,-b*0.5) );
      aVec2.push_back( CVector2(+a*0.5,-b*0.5) );
      aVec2.push_back( CVector2(+a*0.5,+b*0.5) );
      aVec2.push_back( CVector2(-a*0.5,+b*0.5) );
      double theta0 = r4*3.1415*2.0;
      Rotate(aVec2,theta0);
      Translate(aVec2, r2*10-5,r3*10-5);
    }
    CVector2 cg,pa1,pa2;
    double area,I1,I2;
    SecondMomentOfArea_Polygon(cg,area, pa1,I1, pa2,I2,
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
  CVoxelGrid3D vg;
  vg.Add(0,0,0);
  {
    std::vector<double> aXYZ0;
    std::vector<unsigned int> aQuad0;
    vg.GetQuad(aXYZ0, aQuad0);
    EXPECT_EQ(aXYZ0.size(),8*3);
    EXPECT_EQ(aQuad0.size(),6*4);
    {
      std::vector<unsigned int> aTri0;
      convert2Tri_Quad(aTri0, aQuad0);
      EXPECT_EQ(aTri0.size(),12*3);
    }
  }
}

TEST(meshtopo,quad_subdiv1)
{
  CVoxelGrid3D vg;
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
  RemoveUnreferencedPoints_MeshElem(aXYZ0a,aQuad0a, mapInOut,
                                    3,aXYZ0,aQuad0);
  EXPECT_EQ(aXYZ0a.size(),16*3);
  EXPECT_EQ(aQuad0a.size(),14*4);
}


TEST(mathfunc,sherical_harmonics_orthgonality)
{
  std::vector<double> aXYZ,aVal;
  std::vector<unsigned int> aTri;
  MeshTri3D_Cube(aXYZ,aTri, 50);
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
    CVector3 p0(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
    CVector3 p1(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    CVector3 p2(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    double area = SolidAngleTri(p0, p1, p2);
    area_sum += area;
    double a0[N]; makeArray_SphericalHarmonics(a0, norder, p0.x,p0.y,p0.z);
    double a1[N]; makeArray_SphericalHarmonics(a1, norder, p1.x,p1.y,p1.z);
    double a2[N]; makeArray_SphericalHarmonics(a2, norder, p2.x,p2.y,p2.z);
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
