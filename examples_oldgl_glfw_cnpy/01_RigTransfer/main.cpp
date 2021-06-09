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

#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/quat.h"
#include "delfem2/mat3.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/rig_geo3.h"
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <random>


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
  std::vector<dfm2::CRigBone> aBone;
  {
    std::vector<unsigned int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl_Bone(
        aXYZ0, aW, aTri, aIndBoneParent, aJntRgrs,
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
  
  // ---------------------------------
  std::vector<dfm2::CRigBone> aBone_MotionSrc;
  std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
  int nframe = 0;
  std::vector<double> aValueChanelHistoryRotTrans;
  Read_BioVisionHierarchy(
      aBone_MotionSrc,aChannelRotTransBone,nframe,aValueChanelHistoryRotTrans,
//                      std::string(PATH_INPUT_DIR)+"/walk.bvh");
std::string(PATH_INPUT_DIR)+"/jump.bvh");
  
  {
    double scale = (1.0/0.45)*2.54/100.0;
    for(auto & bone : aBone_MotionSrc){
      bone.invBindMat[ 3] *= scale;
      bone.invBindMat[ 7] *= scale;
      bone.invBindMat[11] *= scale;
      bone.transRelative[0] *= scale;
      bone.transRelative[1] *= scale;
      bone.transRelative[2] *= scale;
    }
    unsigned int nch = aChannelRotTransBone.size();
    unsigned int nfrm = aValueChanelHistoryRotTrans.size()/nch;
    for(unsigned int ifrm=0;ifrm<nfrm;++ifrm){
      for(unsigned int ich=0;ich<nch;++ich){
        if( aChannelRotTransBone[ich].isrot ){ continue; }
        aValueChanelHistoryRotTrans[ifrm*nch+ich] *= scale;
      }
    }
  }
  UpdateBoneRotTrans(aBone_MotionSrc);
  std::cout << "nBone:" << aBone_MotionSrc.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  for(unsigned int ibone=0;ibone<aBone_MotionSrc.size();++ibone){
    std::cout << ibone << " " << aBone_MotionSrc[ibone].name << std::endl;
  }
  // ----------------------------------
  
  std::vector<double> aXYZ1 = aXYZ0;

  delfem2::glfw::CViewer3 viewer;
  // -----------
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  dfm2::opengl::setSomeLighting();
  
  std::vector<CRigMap> aMapBoneTrg2Src(aBone.size());
  aMapBoneTrg2Src[0] = CRigMap(0,false);
  aMapBoneTrg2Src[1] = CRigMap(2,false);
  aMapBoneTrg2Src[2] = CRigMap(8,false);
  aMapBoneTrg2Src[4] = CRigMap(3,true);
  aMapBoneTrg2Src[5] = CRigMap(9,true);
  aMapBoneTrg2Src[7] = CRigMap(4,true);
  aMapBoneTrg2Src[8] = CRigMap(10,true);
  aMapBoneTrg2Src[16] = CRigMap(21,false);
  aMapBoneTrg2Src[17] = CRigMap(30,false);
  aMapBoneTrg2Src[18] = CRigMap(22,true);
  aMapBoneTrg2Src[19] = CRigMap(31,true);
  aMapBoneTrg2Src[20] = CRigMap(23,true);
  aMapBoneTrg2Src[21] = CRigMap(32,true);
  
  {
    const unsigned int nBone = aBone.size();
    std::vector<double> aQuat0(nBone*4);
    for(unsigned int ibone=0;ibone<nBone;++ibone){
      dfm2::Quat_Identity(aQuat0.data()+ibone*4);
    }
    for(unsigned int ib=0;ib<nBone;++ib){
      if( !aMapBoneTrg2Src[ib].is_adjust_orientation_parent ){ continue; }
      int jb = aMapBoneTrg2Src[ib].ibone_dist;
      if( jb == -1 ) continue;
      assert( jb < (int)aBone_MotionSrc.size() );
      int jbp = aBone_MotionSrc[jb].ibone_parent;
      dfm2::CVec3d vj = aBone_MotionSrc[jb].Pos() - aBone_MotionSrc[jbp].Pos();
      int ibp = aBone[ib].ibone_parent;
      if( ibp == -1 ) continue;
      dfm2::CVec3d vi(aBone[ib].transRelative);
      const double* qip = aQuat0.data()+ibp*4;
      dfm2::CVec3d vi1 = dfm2::QuatVec(qip, vi);
      dfm2::CMat3d m0 = dfm2::Mat3_MinimumRotation(vi1,vj);
      double qij[4]; m0.GetQuat_RotMatrix(qij);
      double qt[4];
      dfm2::QuatQuat(qt,
                     qij,
                     aQuat0.data()+ibp*4);
      dfm2::Copy_Quat(aQuat0.data()+ib*4,
                      qt);
      dfm2::Copy_Quat(aMapBoneTrg2Src[ibp].quatOffset,
                      qij);
    }
  }

  while ( !glfwWindowShouldClose(viewer.window) )
  {
    { // biovision
      static int iframe = 0;
      const int nch = aChannelRotTransBone.size();
      SetPose_BioVisionHierarchy(aBone_MotionSrc,
                                 aChannelRotTransBone,
                                 aValueChanelHistoryRotTrans.data()+iframe*nch);
      iframe = (iframe+1)%nframe;
    }
    
    {
      const unsigned int nBone = aBone.size();
      for(unsigned int ibone=0;ibone<nBone;++ibone){
        double* q = aBone[ibone].quatRelativeRot;
        q[0] = 1.0;
        q[1] = 0.00;
        q[2] = 0.00;
        q[3] = 0.00;
        dfm2::Normalize_Quat(q);
      }
      for(unsigned int ibs=0;ibs<nBone;++ibs){
        int ibt = aMapBoneTrg2Src[ibs].ibone_dist;
        if( ibt == -1 ) continue;
        double q0[4];
        dfm2::QuatQuat(q0,
                       aBone_MotionSrc[ibt].quatRelativeRot,
                       aMapBoneTrg2Src[ibs].quatOffset);
        dfm2::Copy_Quat(aBone[ibs].quatRelativeRot, q0);
        dfm2::Normalize_Quat(aBone[ibs].quatRelativeRot);
      }
      aBone[0].transRelative[0] = aBone_MotionSrc[0].affmat3Global[3];
      aBone[0].transRelative[1] = aBone_MotionSrc[0].affmat3Global[7];
      aBone[0].transRelative[2] = aBone_MotionSrc[0].affmat3Global[11];
    }

    dfm2::UpdateBoneRotTrans(aBone);
    
    { // surface rigging
      const unsigned int nbone = aBone.size();
      for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip){
        const double* p0 = aXYZ0.data()+ip*3;
        double* p1 = aXYZ1.data()+ip*3;
        p1[0] = 0.0;  p1[1] = 0.0;  p1[2] = 0.0;
        for(unsigned int ibone=0;ibone<nbone;++ibone){
          double p2[3];
          aBone[ibone].DeformSkin(p2, p0);
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
//    dfm2::opengl::DrawJoints(aJntPos1, aIndBoneParent);
//    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
//    dfm2::opengl::DrawJoints(aJntPos0, aIndBoneParent);
    dfm2::opengl::DrawBone_Line(
        aBone,
        -1, -1,
        0.02, 1.0);
    dfm2::opengl::DrawBone_Line(
        aBone_MotionSrc,
        -1, -1,
        0.02, 1.0);
    for(unsigned int ibs=0;ibs<aBone.size();++ibs){
      int ibt = aMapBoneTrg2Src[ibs].ibone_dist;
      if( ibt < 0 ){ continue; }
      ::glLineWidth(3);
      ::glColor3d(1,1,1);
      ::glBegin(GL_LINES);
      dfm2::CVec3d ps = aBone[ibs].Pos();
      dfm2::CVec3d pt = aBone_MotionSrc[ibt].Pos();
      ::glVertex3d(pt.x,pt.y,pt.z);
      ::glVertex3d(ps.x,ps.y,ps.z);
      ::glEnd();
    }
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
