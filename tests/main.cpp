#include <iostream>

#include "gtest/gtest.h"

#include "delfem2/mat3.h"
#include "delfem2/quat.h"
#include "delfem2/msh.h"
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


