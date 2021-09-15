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

#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/quat.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/rigv3.h"

namespace dfm2 = delfem2;

void Draw(
    const std::vector<dfm2::CRigBone>& aBone,
    const std::vector<double>& aXYZ1,
    const std::vector<unsigned int>& aTri)
{
  // -------------------
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_DEPTH_TEST);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
  ::glDisable(GL_DEPTH_TEST);
  delfem2::opengl::DrawBone_Line(
      aBone,
      -1, -1,
      0.01, 1.0);
//    dfm2::opengl::DrawJoints(aJntPos1, aIndBoneParent);
//    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
//    dfm2::opengl::DrawJoints(aJntPos0, aIndBoneParent);
}

int main()
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CRigBone> aBone;
  //
  std::vector<double> aW;
  std::vector<double> aSkinningSparseWeight;
  std::vector<unsigned int> aSkinningSparseIdBone;
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
      for(unsigned int ij=0;ij<aJntPos0.size()/3;++ij){
        std::cout << ij << " " << aIndBoneParent[ij] << " " << aJntPos0[ij*3+0] << " " << aJntPos0[ij*3+1] << " " << aJntPos0[ij*3+2] << std::endl;
      }
      dfm2::InitBones_JointPosition(
          aBone,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
    dfm2::SparsifyMatrixRow(
        aSkinningSparseWeight, aSkinningSparseIdBone,
        aW.data(), aXYZ0.size()/3, aBone.size(),
        1.0e-5);
  }
  std::vector<double> aXYZ1 = aXYZ0;

  // --------
  delfem2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  dfm2::opengl::setSomeLighting();

  std::random_device rnd_dev;
  std::mt19937 rnd_eng(rnd_dev());
  std::uniform_real_distribution<double> dist_01(0,1);

  while ( !glfwWindowShouldClose(viewer.window) )
  {
    for(auto& bone : aBone){
      dfm2::CQuatd::Random(0.2).CopyTo(bone.quatRelativeRot);
    }
    dfm2::CVec3d::Random(dist_01,rnd_eng).CopyToScale(aBone[0].transRelative, 0.2);
    dfm2::UpdateBoneRotTrans(aBone);
    //
    for(int iframe=0;iframe<20;++iframe) {
      if( iframe % 2 == 0 ){
        dfm2::Skinning_LBS(aXYZ1,
            aXYZ0, aBone, aW);
      }
      else{
        dfm2::SkinningSparse_LBS(aXYZ1,
            aXYZ0, aBone, aSkinningSparseWeight, aSkinningSparseIdBone);
      }
      viewer.DrawBegin_oldGL();
      Draw(aBone,aXYZ1,aTri);
      glfwPollEvents();
      glfwSwapBuffers(viewer.window);
      viewer.ExitIfClosed();
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
