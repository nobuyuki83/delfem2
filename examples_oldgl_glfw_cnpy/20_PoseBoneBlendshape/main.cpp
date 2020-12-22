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

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/quat.h"
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <random>

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
  delfem2::opengl::DrawBone(aBone,
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
  std::vector<double> aSkinningSparseWeight;
  std::vector<unsigned int> aSkinningSparseIdBone;
  std::vector<double> aBlendShape;
  std::vector<double> aBlendPose;
  std::vector<unsigned int> aIndBoneParent;
  std::vector<double> aJntRgrsSparseWeight;
  std::vector<unsigned int> aJntRgrsSparseIdbone;
  {
    std::vector<double> aW;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl_BoneBlendshape(
        aXYZ0,
        aW,
        aTri,
        aIndBoneParent,
        aJntRgrs,
        aBlendShape,
        aBlendPose,
        std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    dfm2::SparsifyMatrixRow(
        aSkinningSparseWeight, aSkinningSparseIdBone,
        aW.data(), aXYZ0.size()/3, aIndBoneParent.size(),
        1.0e-5);
    {
      std::vector<double> aJntRgrsT;
      dfm2::Transpose_Mat(aJntRgrsT,
          aJntRgrs, aXYZ0.size() / 3, aIndBoneParent.size());
      dfm2::SparsifyMatrixRow(
          aJntRgrsSparseWeight, aJntRgrsSparseIdbone,
          aJntRgrsT.data(), aIndBoneParent.size(), aXYZ0.size() / 3,
          1.0e-5);
    }
    {
      std::vector<double> aJntPos0;
      dfm2::Points3_WeightsparsePosition(
          aJntPos0,
          aIndBoneParent.size(), aJntRgrsSparseWeight,aJntRgrsSparseIdbone, aXYZ0);
      dfm2::InitBones_JointPosition(
          aBone,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
  }
  std::vector<double> aXYZ1, aXYZ2;
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  dfm2::opengl::setSomeLighting();
  //
  std::random_device rnd_dev;
  std::mt19937 rnd_eng(rnd_dev());
  std::uniform_real_distribution<double> dist_01(0,1);
  //
  for(unsigned int ipose=0;ipose<10;++ipose){
    // input pose
    for(auto& bone : aBone){
      dfm2::CQuatd::Random(0.2).CopyTo(bone.quatRelativeRot);
    }
    dfm2::CVec3d::Random(dist_01,rnd_eng).CopyToScale(aBone[0].transRelative, 0.2);
    //
    for(unsigned int iframe=0;iframe<100;++iframe) {
      double t0 = iframe * 0.01 * 2 * M_PI;
      const double body_param[10] = {
          sin(t0 * 9),
          sin(t0 * 8),
          sin(t0 * 7),
          sin(t0 * 6),
          sin(t0 * 5),
          sin(t0 * 4),
          sin(t0 * 3),
          sin(t0 * 2),
          sin(t0 * 1)};
      aXYZ1.assign(aXYZ0.begin(), aXYZ0.end());
      for (unsigned int i = 0; i < aXYZ0.size(); ++i) { // shape contribution
        for (unsigned int j = 0; j < 10; ++j) {
          aXYZ1[i] += aBlendShape[i * 10 + j] * body_param[j];
        }
      }
      { // pose contribution
        const unsigned int nBone = aBone.size();
        const unsigned int nDofR = (nBone - 1) * 9;
        assert(aBlendPose.size() == aXYZ0.size() * nDofR);
        std::vector<double> arot(nDofR);
        for (unsigned int ib = 1; ib < nBone; ++ib) {
          double mR[9];
          dfm2::Mat3_Quat(mR, aBone[ib].quatRelativeRot);
          mR[0] -= 1.0;
          mR[4] -= 1.0;
          mR[8] -= 1.0;
          memcpy(arot.data() + (ib - 1) * 9, mR, 9 * sizeof(double));
        }
        for (unsigned int i = 0; i < aXYZ0.size(); ++i) {
          for (unsigned int j = 0; j < nDofR; ++j) {
            aXYZ1[i] += aBlendPose[i * nDofR + j] * arot[j];
          }
        }
      }
      {
        std::vector<double> aJntPos0;
        dfm2::Points3_WeightsparsePosition(
            aJntPos0,
            aIndBoneParent.size(), aJntRgrsSparseWeight,aJntRgrsSparseIdbone, aXYZ1);
        dfm2::InitBones_JointPosition(
            aBone,
            aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
      }
      dfm2::UpdateBoneRotTrans(aBone);
      dfm2::SkinningSparse_LBS(aXYZ2,
          aXYZ1, aBone, aSkinningSparseWeight, aSkinningSparseIdBone);
      //
      viewer.DrawBegin_oldGL();
      Draw(aBone, aXYZ2, aTri);
      glfwPollEvents();
      glfwSwapBuffers(viewer.window);
      viewer.ExitIfClosed();
    }
  }
}
