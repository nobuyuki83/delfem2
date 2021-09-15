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
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/defarapenergy_geo3.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/rigopt.h"
#include "delfem2/lsmats.h"
#include "delfem2/mshuni.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/jagarray.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"

namespace dfm2 = delfem2;

// -------------------------------

void Solve_MinEnergyArap(
    std::vector<double>& aXYZ1,
    std::vector<dfm2::CRigBone>& aBone,
    std::vector<double>& aQuat1,
    const std::vector<double>& KMat, // [nsns, nsns]
    //
    const std::vector<dfm2::CTarget>& aTarget,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aW,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup)
{
  
  const double weight_rig = 1000;
  
  class CSystemMatrix_Reduced{
  public:
    CSystemMatrix_Reduced(
        unsigned int nsns_,
        const std::vector<double>& K_, // [ nsns, nsns ]
        unsigned int nC_,
        const std::vector<double>& adC_, // [ nc, nsns]
        double weight_rig_) :
      nsns(nsns_), K(K_),
      nc(nC_), adC(adC_), weight_rig(weight_rig_)
    {
      assert( K.size() == nsns * nsns );
      assert( adC.size() == nc * nsns );
      tmpC0.resize(nc);
      tmpM0.resize(nsns);
    }
  public:
    void MatVec(
        double* y,
        double alpha, const double* x, double beta) const
    {
      for(unsigned int isns=0;isns<nsns;++isns){
        y[isns] = beta*y[isns];
      }
      // ----
      for(unsigned int isns=0;isns<nsns;++isns){
        for(unsigned int jsns=0;jsns<nsns;++jsns){
          y[isns] += alpha*K[isns*nsns+jsns]*x[jsns];
        }
      }
      // ----
      dfm2::MatVec(
          tmpC0.data(),
          adC.data(), nc, nsns, x);
      dfm2::MatTVec(
          tmpM0.data(),
          adC.data(), nc, nsns, tmpC0.data());
      for(unsigned int isns=0;isns<nsns;++isns){
        y[isns] += alpha*weight_rig*tmpM0[isns];
      }
    }
  public:
    const unsigned int nsns;
    const std::vector<double>& K;
    const unsigned int nc;
    const std::vector<double>& adC;
    const double weight_rig;
    // ------
    mutable std::vector<double> tmpC0;
    mutable std::vector<double> tmpM0;
  };
  
  // --------------
  const unsigned int np = aXYZ0.size()/3;
  const unsigned int nb = aBone.size();
  
  // --------------
  std::vector<double> L; // [ nsns, nbone*4 ]
  Rig_Sensitivity_Skeleton(
      L,
      aBone);
  // --------------
  
  const unsigned int nsns = L.size()/(nb*12);
    
  std::vector<double> aC0; // [nC]
  std::vector<double> adC0; // [nC, nsns ]
  for(const auto & trg : aTarget){
    dfm2::Rig_WdW_Target(aC0, adC0,
        aBone, trg, L);
  }
  const unsigned int nc = aC0.size();
  
  // ----------

  CSystemMatrix_Reduced mat(nsns, KMat,
                            nc, adC0,
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
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int j=0;j<nb*4;++j){
        tx[j] += aRefPos[ip*(nb*4)+j]*aResP[ip*3+0];
        ty[j] += aRefPos[ip*(nb*4)+j]*aResP[ip*3+1];
        tz[j] += aRefPos[ip*(nb*4)+j]*aResP[ip*3+2];
      }
    }
    for(unsigned int isns=0;isns<nsns;++isns){
      for(unsigned int jb=0;jb<nb;++jb){
        for(unsigned int jdim=0;jdim<4;++jdim) {
          r[isns]
              += L[isns * (nb * 12) + jb*12+4*0+jdim] * tx[jb*4+jdim]
              +  L[isns * (nb * 12) + jb*12+4*1+jdim] * ty[jb*4+jdim]
              +  L[isns * (nb * 12) + jb*12+4*2+jdim] * tz[jb*4+jdim];
        }
      }
    }
    std::vector<double> rC(nsns,0.0);
    dfm2::MatTVec(rC.data(),
        adC0.data(), nc, nsns, aC0.data());
    for(unsigned int i=0;i<rC.size();++i){ r[i] += weight_rig*rC[i]; }
  }
  
  std::vector<double> u(nsns,0.0);
  std::vector<double> reshist;
  {
    const std::size_t n = nsns;
    std::vector<double> tmp0(n), tmp1(n);
    auto vr = dfm2::CVecXd(r);
    auto vu = dfm2::CVecXd(u);
    auto vs = dfm2::CVecXd(tmp0);
    auto vt = dfm2::CVecXd(tmp1);
    reshist = dfm2::Solve_CG(
        vr, vu, vs, vt,
        1.0e-3, 100, mat);
  }
  
  std::cout << "convergence: " << reshist.size() << std::endl;
  assert( u.size() == (nb+1)*3 );
  for(unsigned int ib=0;ib<aBone.size();++ib){
    dfm2::CVec3d vec_rot(u.data()+ib*3);
    dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(-vec_rot);
    dfm2::CQuatd q0 = dq*dfm2::CQuatd(aBone[ib].quatRelativeRot);
    q0.CopyTo(aBone[ib].quatRelativeRot);
  }
  {
    dfm2::CVec3d vec_trans(u.data()+nb*3);
    aBone[0].transRelative[0] -= vec_trans.x;
    aBone[0].transRelative[1] -= vec_trans.y;
    aBone[0].transRelative[2] -= vec_trans.z;
  }
  dfm2::UpdateBoneRotTrans(aBone);
  dfm2::Skinning_LBS(aXYZ1,
      aXYZ0, aBone, aW);
  for(int itr=0;itr<5;++itr){
    dfm2::UpdateRotationsByMatchingCluster_Linear(aQuat1,
        aXYZ0,aXYZ1,psup_ind,psup);
  }
}



void Draw(
    const std::vector<double>& aXYZ1,
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
    ::glPointSize(20);
    ::glBegin(GL_POINTS);
    for(const auto & it : aTarget){
      const unsigned int ib = it.ib;
      ::glColor3d(0,1,0);
      dfm2::opengl::myGlVertex(aBone[ib].Pos());
      ::glColor3d(0,0,1);
      dfm2::opengl::myGlVertex(it.pos);
    }
    ::glEnd();
  }
  ::glDisable(GL_DEPTH_TEST);
  delfem2::opengl::DrawBone_Line(
      aBone,
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
    std::vector<unsigned int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl_Bone(
        aXYZ0,
        aW,
        aTri,
        aIndBoneParent,
        aJntRgrs,
        std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    {
      std::vector<double> aJntPos0;
      dfm2::Points3_WeighttranspPosition(
          aJntPos0,
          aJntRgrs, aXYZ0);
      dfm2::InitBones_JointPosition(
          aBone,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
  }
  std::vector<unsigned int> psup_ind, psup;
  {
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
        aTri.data(), aTri.size()/3, 3,
        aXYZ0.size()/3);
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
    for(unsigned int ip=0;ip<np;++ip){
      dfm2::Quat_Identity(aQuat1.data()+4*ip);
    }
  }
  
  std::vector<double> K; // [nsns, nsns]
  {
    const unsigned int np = aXYZ1.size()/3;
    dfm2::CMatrixSparse<double> Mat;
    {
      Mat.Initialize(np, 3, true);
      std::vector<unsigned int> psup_ind1, psup1;
      dfm2::JArray_Extend(psup_ind1, psup1,
          psup_ind.data(), psup_ind.size(), psup.data());
      dfm2::JArray_Sort(psup_ind1, psup1);
      assert( psup_ind1.size() == np+1 );
      Mat.SetPattern(psup_ind1.data(), psup_ind1.size(), psup1.data(), psup1.size());
      Mat.setZero();
      std::vector<unsigned int> tmp_buffer;
      for(unsigned int ip=0;ip<np;++ip){
        std::vector<unsigned int> aIP;
        for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
          aIP.push_back(psup[ipsup]);
        }
        aIP.push_back(ip);
        std::vector<double> eM;
        dfm2::ddW_ArapEnergy(
            eM,
            aIP,aXYZ0,aQuat1);
        Mearge(
            Mat,
            aIP.size(), aIP.data(),
            aIP.size(), aIP.data(),
            9, eM.data(),
            tmp_buffer);
      }
    }
    
    std::vector<double> dSkin; // [nsns, np*3]
    {
      std::vector<double> aRefPos; // [ np, nBone*4 ]
      Rig_SkinReferncePositionsBoneWeighted(
          aRefPos,
          aBone,aXYZ0,aW);
      std::vector<double> L; // [ nsns, nbone*4 ]
      { // make sensitivity of bone transformations
        for(unsigned int ibs=0;ibs<aBone.size();++ibs){
          for(int idims=0;idims<3;++idims){
            dfm2::Rig_SensitivityBoneTransform(
                L,
                ibs,idims,true,
                aBone);
          }
        }
        for(int idims=0;idims<3;++idims){
          dfm2::Rig_SensitivityBoneTransform(
              L,
              0,idims,false,
              aBone);
        }
      }
      const unsigned int nb = aBone.size();
      const unsigned int nsns = L.size() / (nb*12);
      dSkin.resize( nsns * (np*3) );
      for(unsigned int isns=0;isns<nsns;++isns){
       for(unsigned int ip=0;ip<np;++ip){
          double dx = 0.0, dy = 0.0, dz = 0.0;
/*
          for(unsigned int j=0;j<nb*4;++j){
            dx += aRefPos[ip*(nb*4)+j]*Lx[isns*(nb*4)+j];
            dy += aRefPos[ip*(nb*4)+j]*Ly[isns*(nb*4)+j];
            dz += aRefPos[ip*(nb*4)+j]*Lz[isns*(nb*4)+j];
          }
*/

          for(unsigned int jb=0;jb<nb;++jb){
            for(unsigned int jdim=0;jdim<4;++jdim) {
              dx += aRefPos[ip * (nb * 4) + jb*4+jdim] * L[isns * (nb * 12) + jb*12+0*4+jdim];
              dy += aRefPos[ip * (nb * 4) + jb*4+jdim] * L[isns * (nb * 12) + jb*12+1*4+jdim];
              dz += aRefPos[ip * (nb * 4) + jb*4+jdim] * L[isns * (nb * 12) + jb*12+2*4+jdim];
            }
          }

          dSkin[ isns*(np*3) + ip*3+0 ] = dx;
          dSkin[ isns*(np*3) + ip*3+1 ] = dy;
          dSkin[ isns*(np*3) + ip*3+2 ] = dz;
        }
      }
    }
    
    const unsigned int nsns = dSkin.size()/(np*3);
    K.resize( nsns*nsns );
    for(unsigned int isns=0;isns<nsns;++isns){
      const std::vector<double> vi(dSkin.data()+isns*(np*3), dSkin.data()+(isns+1)*(np*3));
      std::vector<double> t0(np*3);
      Mat.MatVec(t0.data(), 1.0, vi.data(), 0.0);
      // -----
      for(unsigned int jsns=0;jsns<nsns;++jsns){
//        if( jsns > isns ){ continue; }
        const std::vector<double> vj(dSkin.data()+jsns*(np*3), dSkin.data()+(jsns+1)*(np*3));
        // -----
        K[isns*nsns+jsns] = dfm2::Dot(t0, vj);
        K[jsns*nsns+isns] = K[isns*nsns+jsns];
      }
    }
    for(unsigned int isns=0;isns<nsns;++isns){
      std::cout << isns << " " << K[isns*nsns+isns] << std::endl;
    }
    for(unsigned int isns=0;isns<nsns;++isns){
      for(unsigned int jsns=0;jsns<nsns;++jsns){
        std::cout << isns << " " << jsns << " " <<  K[jsns*nsns+isns] << " " << K[isns*nsns+jsns] << std::endl;
      }
    }
    std::cout << "nsns: " << nsns << std::endl;
  }
  
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
  for(auto & it : aTarget){
    aTargetOriginPos.push_back(it.pos);
  }

  // -------------
  // below: opengl
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  viewer.camera.view_height = 1.2;
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
