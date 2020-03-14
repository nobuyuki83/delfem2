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
#include "delfem2/vecxitrsol.h"
#include "delfem2/quat.h"
#include "delfem2/mat3.h"
#include "delfem2/v23m3q.h"
//
#include "delfem2/rig_v3q.h"
#include "delfem2/cnpy/smpl_cnpy.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_rig_v23q.h"

namespace dfm2 = delfem2;

class CRigMap {
public:
  CRigMap():ibone_dist(-1), is_adjust_orientation_parent(false){
    dfm2::Quat_Identity(quatOffset);
  }
  CRigMap(int ibone_dist, bool is_adjust_orientation_parent)
  :ibone_dist(ibone_dist), is_adjust_orientation_parent(is_adjust_orientation_parent){
    dfm2::Quat_Identity(quatOffset);
  }
public:
  int ibone_dist;
  bool is_adjust_orientation_parent;
  double quatOffset[4];
};


int main()
{
  std::vector<unsigned int> aTri;
  std::vector<double> aXYZ0;
  std::vector<double> aW;      
  std::vector<int> aIndBoneParent;
  std::vector<double> aJntPos0;
  dfm2::cnpy::LoadSmpl(aXYZ0, aW, aTri, aIndBoneParent, aJntPos0,
                       std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
  
  // ---------------------------------
  std::vector<dfm2::CRigBone> aBone;
  std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
  int nframe = 0;
  std::vector<double> aValueChanelHistoryRotTrans;
  Read_BioVisionHierarchy(aBone,aChannelRotTransBone,nframe,aValueChanelHistoryRotTrans,
//                          std::string(PATH_INPUT_DIR)+"/walk.bvh");
                          std::string(PATH_INPUT_DIR)+"/jump.bvh");
  
  {
    double scale = (1.0/0.45)*2.54/100.0;
    for(int ibone=0;ibone<aBone.size();++ibone){
      aBone[ibone].invBindMat[ 3] *= scale;
      aBone[ibone].invBindMat[ 7] *= scale;
      aBone[ibone].invBindMat[11] *= scale;
      aBone[ibone].trans[0] *= scale;
      aBone[ibone].trans[1] *= scale;
      aBone[ibone].trans[2] *= scale;
    }
    unsigned int nch = aChannelRotTransBone.size();
    unsigned int nfrm = aValueChanelHistoryRotTrans.size()/nch;
    for(int ifrm=0;ifrm<nfrm;++ifrm){
      for(int ich=0;ich<nch;++ich){
        if( aChannelRotTransBone[ich].isrot ){ continue; }
        aValueChanelHistoryRotTrans[ifrm*nch+ich] *= scale;
      }
    }
  }
  UpdateBoneRotTrans(aBone);
  std::cout << "nBone:" << aBone.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  for(int ibone=0;ibone<aBone.size();++ibone){
    std::cout << ibone << " " << aBone[ibone].name << std::endl;
  }
  // ----------------------------------
  
  std::vector<double> aMat4AffineBone; // global affine matrix of a bone
  {
    const unsigned int nbone = aIndBoneParent.size();
    aMat4AffineBone.resize(nbone*16);
    for(unsigned int ibone=0;ibone<nbone;++ibone){
      dfm2::Mat4_Identity(aMat4AffineBone.data()+ibone*16);
    }
  }
  std::vector<double> aQuatRelativeRot;
  {
    const unsigned int nbone = aIndBoneParent.size();
    aQuatRelativeRot.resize(nbone*4);
    for(unsigned int ibone=0;ibone<nbone;++ibone){
      dfm2::Quat_Identity(aQuatRelativeRot.data()+ibone*4);
    }
  }
  std::vector<double> aXYZ1 = aXYZ0;
  std::vector<double> aJntPos1 = aJntPos0;
    
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  dfm2::opengl::setSomeLighting();
  
  std::vector<CRigMap> aMapBoneSrcTrg(aIndBoneParent.size());
  aMapBoneSrcTrg[0] = CRigMap(0,false);
  aMapBoneSrcTrg[1] = CRigMap(2,false);
  aMapBoneSrcTrg[2] = CRigMap(8,false);
  aMapBoneSrcTrg[4] = CRigMap(3,true);
  aMapBoneSrcTrg[5] = CRigMap(9,true);
  aMapBoneSrcTrg[7] = CRigMap(4,true);
  aMapBoneSrcTrg[8] = CRigMap(10,true);
  aMapBoneSrcTrg[16] = CRigMap(21,false);
  aMapBoneSrcTrg[17] = CRigMap(30,false);
  aMapBoneSrcTrg[18] = CRigMap(22,true);
  aMapBoneSrcTrg[19] = CRigMap(31,true);
  aMapBoneSrcTrg[20] = CRigMap(23,true);
  aMapBoneSrcTrg[21] = CRigMap(32,true);
  
  {
    const unsigned int nBone = aIndBoneParent.size();
    std::vector<double> aQuat0(nBone*4);
    for(int ibone=0;ibone<nBone;++ibone){
      dfm2::Quat_Identity(aQuat0.data()+ibone*4);
    }
    for(int ib=0;ib<nBone;++ib){
      if( !aMapBoneSrcTrg[ib].is_adjust_orientation_parent ){ continue; }
      int jb = aMapBoneSrcTrg[ib].ibone_dist;
      if( jb == -1 ) continue;
      assert( jb < aBone.size() );
      int jbp = aBone[jb].ibone_parent;
      dfm2::CVec3d vj = aBone[jb].Pos() - aBone[jbp].Pos();
      int ibp = aIndBoneParent[ib];
      if( ibp == -1 ) continue;
      dfm2::CVec3d vi = dfm2::CVec3d(aJntPos0.data()+ib*3) - dfm2::CVec3d(aJntPos0.data()+ibp*3);
      const double* qip = aQuat0.data()+ibp*4;
      dfm2::CVec3d vi1 = dfm2::QuatConjVec(qip, vi);
      dfm2::CMat3d m0 = dfm2::Mat3_MinimumRotation(vi1,vj);
      double qij[4]; m0.GetQuat_RotMatrix(qij);
      double qt[4];
      dfm2::QuatQuat(qt,
                     qij,
                     aQuat0.data()+ibp*4);
      dfm2::Copy_Quat(aQuat0.data()+ib*4,
                      qt);
      dfm2::Copy_Quat(aMapBoneSrcTrg[ibp].quatOffset,
                      qij);
    }
  }

  while (true)
  {
    { // biovision
      static int iframe = 0;
      const int nch = aChannelRotTransBone.size();
      SetPose_BioVisionHierarchy(aBone,
                                 aChannelRotTransBone,
                                 aValueChanelHistoryRotTrans.data()+iframe*nch);
      iframe = (iframe+1)%nframe;
    }
    
    {
      const unsigned int nBone = aIndBoneParent.size();
      for(int ibone=0;ibone<nBone;++ibone){
        double* q = aQuatRelativeRot.data()+ibone*4;
        q[0] = 1.0;
        q[1] = 0.00;
        q[2] = 0.00;
        q[3] = 0.00;
        dfm2::Normalize_Quat(q);
      }
      for(int ibs=0;ibs<nBone;++ibs){
        int ibt = aMapBoneSrcTrg[ibs].ibone_dist;
        if( ibt == -1 ) continue;
        double q0[4];
        dfm2::QuatQuat(q0,
                       aBone[ibt].quatRelativeRot,
                       aMapBoneSrcTrg[ibs].quatOffset);
        dfm2::Copy_Quat(aQuatRelativeRot.data()+ibs*4, q0);
        dfm2::Normalize_Quat(aQuatRelativeRot.data()+ibs*4);
      }
      const double trans_root[3] = {
        aBone[0].affmat3Global[3],
        aBone[0].affmat3Global[7],
        aBone[0].affmat3Global[11],
      };
      dfm2::SetMat4AffineBone_FromJointRelativeRotation(aMat4AffineBone,
                                                        trans_root, aQuatRelativeRot, aIndBoneParent, aJntPos0);
    }

  
    for(int ibone=0;ibone<aIndBoneParent.size();++ibone){
      dfm2::Vec3_AffMat3Vec3Projection(aJntPos1.data()+ibone*3,
                                       aMat4AffineBone.data()+ibone*16, aJntPos0.data()+ibone*3);
    }
    
    { // surface rigging
      const unsigned int nbone = aIndBoneParent.size();
      for(int ip=0;ip<aXYZ0.size()/3;++ip){
        const double* p0 = aXYZ0.data()+ip*3;
        double* p1 = aXYZ1.data()+ip*3;
        p1[0] = 0.0;  p1[1] = 0.0;  p1[2] = 0.0;
        for(int ibone=0;ibone<nbone;++ibone){
          double p2[3];
          dfm2::Vec3_AffMat3Vec3Projection(p2,
                                           aMat4AffineBone.data()+ibone*16, p0);
          p1[0] += aW[ip*nbone+ibone]*p2[0];
          p1[1] += aW[ip*nbone+ibone]*p2[1];
          p1[2] += aW[ip*nbone+ibone]*p2[2];
        }
      }
    }
    
   
    // -------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_DEPTH_TEST);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
    dfm2::opengl::DrawJoints(aJntPos1, aIndBoneParent);
//    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
//    dfm2::opengl::DrawJoints(aJntPos0, aIndBoneParent);
    dfm2::opengl::DrawBone(aBone,
                           -1, -1,
                           0.02, 1.0);
    for(int ibs=0;ibs<aIndBoneParent.size();++ibs){
      int ibt = aMapBoneSrcTrg[ibs].ibone_dist;
      if( ibt < 0 ){ continue; }
      ::glLineWidth(3);
      ::glColor3d(1,1,1);
      ::glBegin(GL_LINES);
      ::glVertex3d(aJntPos1[ibs*3+0], aJntPos1[ibs*3+1], aJntPos1[ibs*3+2]);
      dfm2::CVec3d p = aBone[ibt].Pos();
      ::glVertex3d(p.x(),p.y(),p.z());
      ::glEnd();
    }
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
