/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief SMPL model
 * @details skinning
 */

#include <cstdlib>
#include <random>
#include <time.h>
#include <GLFW/glfw3.h>
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
//
#include "delfem2/v23m3q.h"
#include "delfem2/objfunc_v23.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

#include "delfem2/rig_v3q.h"
#include "delfem2/opengl/glold_rig_v23q.h"

#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -------------------

void Check_SensitivityRigSkin(
    std::vector<dfm2::CRigBone> aBone1,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aW)
{
  const unsigned int nb = aBone1.size();
  assert( aW.size() == aXYZ0.size()/3*nb );
  
  // ------------
  std::vector<double> Lx, Ly, Lz;  // [ nsns, nb*4 ]
  for(int ibs=0;ibs<aBone1.size();++ibs){
    for(int idims=0;idims<3;++idims){
      dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                               ibs,idims,true,
                                               aBone1);
    }
  }
  for(int idims=0;idims<3;++idims){
    dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                             0,idims,false,
                                             aBone1);
  }
  // ---------------
  
  std::vector<double> aRefPos; // [ np, nBone*4 ]
  Rig_SkinReferncePositionsBoneWeighted(aRefPos,
                                        aBone1,aXYZ0,aW);
  
  const double eps = 1.0e-4;
  const unsigned int nsns = Lx.size()/(nb*4);
  assert( nsns==(nb+1)*3 );
  for(int isns=0;isns<nsns;++isns){
    unsigned int ib_s = isns/3;
    bool is_rot = true;
    unsigned int idim_s = isns - ib_s*3;
    if( ib_s == nb ){ ib_s = 0; is_rot = false; }
    std::vector<dfm2::CRigBone> aBone2 = aBone1;
    if( is_rot ){
      dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(eps*dfm2::CVec3d::Axis(idim_s));
      dfm2::CQuatd q0 = dq*dfm2::CQuatd(aBone2[ib_s].quatRelativeRot);
      q0.CopyTo(aBone2[ib_s].quatRelativeRot);
    }
    else{
      aBone2[ib_s].transRelative[idim_s] += eps;
    }
    std::vector<double> aXYZ1;
    dfm2::UpdateBoneRotTrans(aBone1);
    dfm2::Skinning_LBS(aXYZ1,
                       aXYZ0, aBone1, aW);
    // ----------------
    std::vector<double> aXYZ2;
    dfm2::UpdateBoneRotTrans(aBone2);
    dfm2::Skinning_LBS(aXYZ2,
                       aXYZ0, aBone2, aW);
    {
      const unsigned int np = aXYZ0.size()/3;
      // ----------------
      double max_ratio = 0.0;
      for(int ip=0;ip<np;++ip){
        const double val0[3] = {
          (aXYZ2[ip*3+0] - aXYZ1[ip*3+0])/eps,
          (aXYZ2[ip*3+1] - aXYZ1[ip*3+1])/eps,
          (aXYZ2[ip*3+2] - aXYZ1[ip*3+2])/eps };
        double val1[3] =  { 0, 0, 0 };
        for(int j=0;j<nb*4;++j){
          val1[0] += aRefPos[ip*(nb*4)+j]*Lx[isns*(nb*4)+j];
          val1[1] += aRefPos[ip*(nb*4)+j]*Ly[isns*(nb*4)+j];
          val1[2] += aRefPos[ip*(nb*4)+j]*Lz[isns*(nb*4)+j];
        }
        for(int i=0;i<3;++i){
          double ratio = fabs(val0[i]-val1[i])/(fabs(val1[i])+1.0);
          max_ratio = (ratio >max_ratio ) ? ratio : max_ratio;
        }
      }
      std::cout << "check sensitivity skin: max error ratio: " << ib_s << " " << idim_s << " " << max_ratio << std::endl;
    }
  }
}

// -------------------------------

void Solve_MinEnergyArap
(std::vector<double>& aXYZ1,
 std::vector<dfm2::CRigBone>& aBone,
 std::vector<double>& aQuat1,
 const std::vector<double>& KMat, // [3*nb*4, 3*nb*4]
 //
 const std::vector<dfm2::CTarget>& aTarget,
 const std::vector<double>& aXYZ0,
 const std::vector<double>& aW,
 const std::vector<unsigned int>& psup_ind,
 const std::vector<unsigned int>& psup)
{
  
  const double weight_rig = 1;
  
  class CSystemMatrix_Reduced{
  public:
    CSystemMatrix_Reduced(const std::vector<double>& Lx_, // [ nsns, nb*4 ]
                          const std::vector<double>& Ly_, // [ nsns, nb*4 ]
                          const std::vector<double>& Lz_, // [ nsns, nb*4 ]
                          unsigned int nb_,
                          const std::vector<double>& K_, // [ 3*nb*4, 3*nb*4 ]
                          unsigned int nMode_,
                          unsigned int np_,
                          const std::vector<double>& adC_,
                          unsigned int nC_,
                          double weight_rig_) :
    Lx(Lx_), Ly(Ly_), Lz(Lz_), nb(nb_), // [np, nsns]
    K(K_),
    nsns(nMode_), np(np_), adC(adC_), nc(nC_), weight_rig(weight_rig_)
    {
//      std::cout << "constructor reduced system matrix " << std::endl;
      assert( Lx.size() == nsns*nb*4 );
      assert( Ly.size() == nsns*nb*4 );
      assert( Lz.size() == nsns*nb*4 );
      assert( K.size() == 3*nb*4 * 3*nb*4 );
      tmpD0.resize(np*3);
      tmpD1.resize(np*3);
      tmpC0.resize(nc);
      tmpM0.resize(nsns);
    }
  public:
    void MatVec(double* y,
                double alpha, const double* x, double beta) const {
      std::vector<double> t2(nsns,0.0);
      {
        std::vector<double> t0(3*nb*4,0.0);
        for(int isns=0;isns<nsns;++isns){
          for(int j=0;j<nb*4;++j){
            t0[0*nb*4+j] += Lx[isns*(nb*4)+j]*x[isns];
            t0[1*nb*4+j] += Ly[isns*(nb*4)+j]*x[isns];
            t0[2*nb*4+j] += Lz[isns*(nb*4)+j]*x[isns];
          }
        }
        std::vector<double> t1(3*nb*4,0.0);
        dfm2::MatVec(t1.data(), K.data(), 3*nb*4, 3*nb*4, t0.data());
        for(int isns=0;isns<nsns;++isns){
          t2[isns] = 0.0;
          for(int j=0;j<nb*4;++j){
            t2[isns]
            += Lx[isns*(nb*4)+j]*t1[0*nb*4+j]
            +  Ly[isns*(nb*4)+j]*t1[1*nb*4+j]
            +  Lz[isns*(nb*4)+j]*t1[2*nb*4+j];
          }
        }
      }
      // ----
      dfm2::MatVec(tmpC0.data(),
                   adC.data(), nc, nsns, x);
      dfm2::MatTVec(tmpM0.data(),
                    adC.data(), nc, nsns, tmpC0.data());
      for(int i=0;i<nsns;++i){ y[i] = beta*x[i] + alpha*(weight_rig*tmpM0[i]+t2[i]); }
    }
  public:
    const std::vector<double>& Lx;
    const std::vector<double>& Ly;
    const std::vector<double>& Lz;
    unsigned int nb;
    const std::vector<double>& K;
    unsigned int nsns;
    unsigned int np;
    const std::vector<double>& adC;
    unsigned int nc;
    double weight_rig;
    mutable std::vector<double> tmpD0;
    mutable std::vector<double> tmpD1;
    mutable std::vector<double> tmpC0;
    mutable std::vector<double> tmpM0;
  };
  
  // --------------
  const unsigned int np = aXYZ0.size()/3;
  const unsigned int nb = aBone.size();
  
  // --------------
  std::vector<double> Lx, Ly, Lz; // [ nsns, nbone*4 ]
  for(int ibs=0;ibs<aBone.size();++ibs){
    for(int idims=0;idims<3;++idims){
      dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                               ibs,idims,true,
                                               aBone);
    }
  }
  for(int idims=0;idims<3;++idims){
    dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                             0,idims,false,
                                             aBone);
  }
  // --------------
  
  const unsigned int nsns = Lx.size()/(nb*4);
    
  std::vector<double> aC0; // [nC]
  std::vector<double> adC0; // [nC, nsns ]
  for(int it=0;it<aTarget.size();++it){
    dfm2::Rig_WdW_Target_Eigen(aC0, adC0,
                               aBone, aTarget[it], Lx, Ly, Lz);
  }
  const unsigned int nc = aC0.size();
  
  // ----------

  CSystemMatrix_Reduced mat(Lx, Ly, Lz, nb, KMat,
                            nsns,np,
                            adC0, nc,
                            weight_rig);
  
  std::vector<double> r(nsns,0.0);
  {
    std::vector<double> aResP;
    dfm2::dW_ArapEnergy(aResP,
                        aXYZ0, aXYZ1, aQuat1, psup_ind, psup);
    std::vector<double> aRefPos; // [ np, nBone*4 ]
    Rig_SkinReferncePositionsBoneWeighted(aRefPos,
                                          aBone,aXYZ0,aW);
    std::vector<double> tx(nb*4, 0.0);
    std::vector<double> ty(nb*4, 0.0);
    std::vector<double> tz(nb*4, 0.0);
    for(int ip=0;ip<np;++ip){
      for(int j=0;j<nb*4;++j){
        tx[j] += aRefPos[ip*(nb*4)+j]*aResP[ip*3+0];
        ty[j] += aRefPos[ip*(nb*4)+j]*aResP[ip*3+1];
        tz[j] += aRefPos[ip*(nb*4)+j]*aResP[ip*3+2];
      }
    }
    for(int isns=0;isns<nsns;++isns){
      for(int j=0;j<nb*4;++j){
        r[isns]
        += Lx[isns*(nb*4)+j]*tx[j]
        +  Ly[isns*(nb*4)+j]*ty[j]
        +  Lz[isns*(nb*4)+j]*tz[j];
      }
    }
    std::vector<double> rC(nsns,0.0);
    dfm2::MatTVec(rC.data(),
                  adC0.data(), nc, nsns, aC0.data());
    for(int i=0;i<rC.size();++i){ r[i] += weight_rig*rC[i]; }
  }
  
  std::vector<double> u(nsns,0.0);
  std::vector<double> reshist = dfm2::Solve_CG(r.data(), u.data(),
                                               nsns, 1.0e-3, 100, mat);
  
  std::cout << "convergence: " << reshist.size() << std::endl;
  assert( u.size() == (nb+1)*3 );
  for(int ib=0;ib<aBone.size();++ib){
    dfm2::CVec3d vec_rot(u.data()+ib*3);
    dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(-vec_rot);
    dfm2::CQuatd q0 = dq*dfm2::CQuatd(aBone[ib].quatRelativeRot);
    q0.CopyTo(aBone[ib].quatRelativeRot);
  }
  {
    dfm2::CVec3d vec_trans(u.data()+nb*3);
    aBone[0].transRelative[0] -= vec_trans.x();
    aBone[0].transRelative[1] -= vec_trans.y();
    aBone[0].transRelative[2] -= vec_trans.z();
  }
  dfm2::UpdateBoneRotTrans(aBone);
  dfm2::Skinning_LBS(aXYZ1,
                     aXYZ0, aBone, aW);
  for(int itr=0;itr<5;++itr){
    dfm2::UpdateRotationsByMatchingCluster(aQuat1,
                                           aXYZ0,aXYZ1,psup_ind,psup);
  }
}



void Draw
(const std::vector<double>& aXYZ1,
 const std::vector<unsigned int>& aTri,
 const std::vector<dfm2::CRigBone>& aBone,
 const std::vector<dfm2::CTarget>& aTarget)
{
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_DEPTH_TEST);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
  //    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
  { // draw bone
    ::glDisable(GL_DEPTH_TEST);
    ::glDisable(GL_LIGHTING);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    for(int it=0;it<aTarget.size();++it){
      const unsigned int ib = aTarget[it].ib;
      ::glColor3d(0,1,0);
      dfm2::opengl::myGlVertex(aBone[ib].Pos());
      ::glColor3d(0,0,1);
      dfm2::opengl::myGlVertex(aTarget[it].pos);
    }
    ::glEnd();
  }
  ::glDisable(GL_DEPTH_TEST);
  delfem2::opengl::DrawBone(aBone,
                            -1, -1,
                            0.01, 1.0);
  ::glEnable(GL_DEPTH_TEST);
  //    dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat1,0.02);
}

// --------------------

int main()
{
  std::vector<double> aXYZ0;
  std::vector<double> aW;
  std::vector<dfm2::CRigBone> aBone;
  std::vector<unsigned int> aTri;
  {
    std::vector<int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl(aXYZ0,
                         aW,
                         aTri,
                         aIndBoneParent,
                         aJntRgrs,
                         std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    dfm2::Smpl2Rig(aBone,
                   aIndBoneParent, aXYZ0, aJntRgrs);
    
  }
  std::vector<unsigned int> psup_ind, psup;
  {
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aTri.data(), aTri.size()/3, 3,
                               (int)aXYZ0.size()/3);
    dfm2::JArray_Sort(psup_ind, psup);
  }
    
  std::vector<double> aXYZ1 = aXYZ0;
  { // initalize pose
    dfm2::UpdateBoneRotTrans(aBone);
    dfm2::Skinning_LBS(aXYZ1,
                       aXYZ0, aBone, aW);
  }
  
  std::vector<double> aQuat1;
  { // initialize rotation
    const unsigned int np = aXYZ1.size()/3;
    aQuat1.resize(np*4);
    for(int ip=0;ip<np;++ip){
      dfm2::Quat_Identity(aQuat1.data()+4*ip);
    }
  }
  
  std::vector<double> K; // [3*nb*4, 3*nb*4]
  {
    const unsigned int np = aXYZ1.size()/3;
    dfm2::CMatrixSparse<double> Mat;
    {
      Mat.Initialize(np, 3, true);
      std::vector<unsigned int> psup_ind1, psup1;
      dfm2::JArray_Extend(psup_ind1, psup1,
                          psup_ind, psup);
      dfm2::JArray_Sort(psup_ind1, psup1);
      assert( psup_ind1.size() == np+1 );
      Mat.SetPattern(psup_ind1.data(), psup_ind1.size(), psup1.data(), psup1.size());
      Mat.SetZero();
      std::vector<int> tmp_buffer;
      for(unsigned int ip=0;ip<np;++ip){
        std::vector<unsigned int> aIP;
        for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
          aIP.push_back(psup[ipsup]);
        }
        aIP.push_back(ip);
        std::vector<double> eM;
        dfm2::ddW_ArapEnergy(eM,
                             aIP,aXYZ0,aQuat1);
        Mat.Mearge(aIP.size(), aIP.data(),
                   aIP.size(), aIP.data(),
                   9, eM.data(),
                   tmp_buffer);
      }
    }
    std::vector<double> aRefPos; // [ np, nBone*4 ]
    Rig_SkinReferncePositionsBoneWeighted(aRefPos,
                                          aBone,aXYZ0,aW);
    const unsigned int M = aBone.size()*4;
    K.resize(3*M*3*M);
    for(int i=0;i<3;++i){
    for(int mi=0;mi<M;++mi){
      std::vector<double> vi(np*3);
      for(int ip=0;ip<np;++ip){ vi[ip*3+i] = aRefPos[ip*M+mi]; }
      std::vector<double> t0(np*3);
      Mat.MatVec(t0.data(), 1.0, vi.data(), 0.0);
      // -----
      for(int j=0;j<3;++j){
      for(int mj=0;mj<M;++mj){
        if( j*M+mj > i*M+mi ){ continue; }
        std::vector<double> vj(np*3);
        for(int jp=0;jp<np;++jp){ vj[jp*3+j] = aRefPos[jp*M+mj]; }
        // -----
        K[(i*M+mi)*(3*M)+(j*M+mj)] = dfm2::Dot(t0, vj);
        K[(j*M+mj)*(3*M)+(i*M+mi)] = K[(i*M+mi)*(3*M)+(j*M+mj)];
      }
      }
      
    }
    }
  }
    
  // -----------
  
  Check_SensitivityRigSkin(aBone,aXYZ0,aW);
  
  // -----------

  std::vector<dfm2::CTarget> aTarget;
  {
    {
      dfm2::CTarget t;
      t.ib = 20;
      t.pos = aBone[t.ib].Pos();
      aTarget.push_back(t);
    }
    {
      dfm2::CTarget t;
      t.ib = 10;
      t.pos = aBone[t.ib].Pos();
      aTarget.push_back(t);
    }
  }
  std::vector< dfm2::CVec3d > aTargetOriginPos;
  for(int it=0;it<aTarget.size();++it){
    aTargetOriginPos.push_back(aTarget[it].pos);
  }
 
  // -----------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.camera_rot_mode = dfm2::CAMERA_ROT_TBALL;
  viewer.nav.camera.view_height = 1.2;
  dfm2::opengl::setSomeLighting();

  int iframe = 0;
  while (true)
  {
    iframe++;
    {
      aTarget[0].pos = aTargetOriginPos[0] + 0.4*dfm2::CVec3d(1-cos(0.1*iframe), sin(0.1*iframe), 0.0);
      aTarget[1].pos = aTargetOriginPos[1] + 0.1*dfm2::CVec3d(sin(0.1*iframe), 1-cos(0.1*iframe), 0.0);
      Solve_MinEnergyArap(aXYZ1,aBone,aQuat1,K,
                          aTarget,
                          aXYZ0,aW,psup_ind,psup);
      Skinning_LBS(aXYZ1,
                   aXYZ0, aBone, aW);
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    Draw(aXYZ1,aTri,aBone,aTarget);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
