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
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/rigopt.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

// -------------------

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
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    for(const auto & it : aTarget){
      const unsigned int ib = it.ib;
      ::glColor3d(1,0,0);
      dfm2::opengl::myGlVertex(aBone[ib].RootPosition());
      ::glColor3d(1,0,0);
      dfm2::opengl::myGlVertex(it.pos);
    }
    ::glEnd();
  }
  /*
   ::glDisable(GL_DEPTH_TEST);
   delfem2::opengl::DrawBone(aBone,
   -1, -1,
   0.01, 1.0);
   */
  ::glEnable(GL_DEPTH_TEST);
  //    dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat1,0.02);
}

// --------------------

int main()
{
  std::vector<double> vtx_xyz_ini;
  std::vector<double> aW;
  std::vector<dfm2::CRigBone> bones;
  std::vector<unsigned int> tri_vtx;
  {
    std::vector<unsigned int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl_Bone(
        vtx_xyz_ini,
        aW,
        tri_vtx,
        aIndBoneParent,
        aJntRgrs,
        std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    {
      std::vector<double> aJntPos0;
      dfm2::Points3_WeighttranspPosition(
          aJntPos0,
          aJntRgrs, vtx_xyz_ini);
      dfm2::InitBones_JointPosition(
          bones,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
  }
  
  std::vector<double> vtx_xyz = vtx_xyz_ini;
  { // initalize pose
    /*
    for(int ibone=0;ibone<aBone.size();++ibone){
      dfm2::CQuatd::Random(0.2).CopyTo(aBone[ibone].quatRelativeRot);
    }
     */
    dfm2::UpdateBoneRotTrans(bones);
    dfm2::Skinning_LBS(
        vtx_xyz,
        vtx_xyz_ini, bones, aW);
  }
  
  std::vector<dfm2::CTarget> aTarget;
  {
    {
      dfm2::CTarget t;
      t.ib = 20;
      t.pos = bones[t.ib].RootPosition();
      aTarget.push_back(t);
    }
    {
      dfm2::CTarget t;
      t.ib = 10;
      t.pos = bones[t.ib].RootPosition();
      aTarget.push_back(t);
    }
  }
  std::vector< dfm2::CVec3d > aTargetOriginPos;
  aTargetOriginPos.reserve(aTarget.size());
  for(auto & trg : aTarget){
    aTargetOriginPos.push_back(trg.pos);
  }

  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  // viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  dfm2::opengl::setSomeLighting();

  int iframe = 0;
  while ( !glfwWindowShouldClose(viewer.window) )
  {
    iframe++;
    {
      aTarget[0].pos = aTargetOriginPos[0] + 0.4*dfm2::CVec3d(1-cos(0.1*iframe), sin(0.1*iframe), 0.0);
      aTarget[1].pos = aTargetOriginPos[1] + 0.1*dfm2::CVec3d(sin(0.1*iframe), 1-cos(0.1*iframe), 0.0);
      Solve_MinRigging(bones, aTarget);
      Skinning_LBS(vtx_xyz,
                   vtx_xyz_ini, bones, aW);
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    Draw(vtx_xyz, tri_vtx, bones, aTarget);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
