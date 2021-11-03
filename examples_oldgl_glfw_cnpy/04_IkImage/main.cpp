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
#include "delfem2/rigopt.h"
#include "delfem2/rig_geo3.h"
#include "inputs_imgboneloc.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/tex.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

namespace dfm2 = delfem2;

// -------------------

void Draw(
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &aTri,
    const std::vector<dfm2::CRigBone> &aBone,
    const std::vector<dfm2::CTarget> &aTarget) {
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_DEPTH_TEST);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size() / 3);
  //    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
  ::glDisable(GL_DEPTH_TEST);
  ::glDisable(GL_LIGHTING);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for (const auto &it : aTarget) {
    const unsigned int ib = it.ib;
    ::glColor3d(1, 0, 0);
    dfm2::opengl::myGlVertex(aBone[ib].RootPosition());
  }
  ::glEnd();
  // ------
  ::glEnable(GL_DEPTH_TEST);
  ::glBegin(GL_LINES);
  ::glColor3d(1, 0, 0);
  for (const auto &it : aTarget) {
    dfm2::CVec3d p = it.pos;
    dfm2::opengl::myGlVertex(p + 10.0 * dfm2::CVec3d(0, 0, 1));
    dfm2::opengl::myGlVertex(p - 10.0 * dfm2::CVec3d(0, 0, 1));
  }
  ::glEnd();
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

int main() {
  std::vector<dfm2::CTarget> aTarget;
  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    std::string name_img_in_test_inputs;
    std::vector<std::pair<double, dfm2::CVec2d> > aBoneLoc;
    double scale = 1.0;
    BoneLocs_SmplUglysweater(name_img_in_test_inputs,
                             scale,
                             aBoneLoc);
    //--
    int width, height;
    {
      int channels;
      unsigned char *img = stbi_load((std::string(PATH_INPUT_DIR) + "/" + name_img_in_test_inputs).c_str(),
                                     &width, &height, &channels, 0);
      tex.Initialize(width, height, channels, img);
      stbi_image_free(img);
      tex.max_x = -scale * width * 0.5;
      tex.min_x = +scale * width * 0.5;
      tex.max_y = -scale * height * 0.5;
      tex.min_y = +scale * height * 0.5;
      tex.z = -0.5;
    }
    aTarget.clear();
    for (auto &it : aBoneLoc) {
      dfm2::CTarget t;
      t.ib = it.first;
      int iw = it.second.x;
      int ih = it.second.y;
      t.pos.p[0] = (double) iw / width - 0.5;
      t.pos.p[1] = 0.5 * height / width - (double) ih / height;
      aTarget.push_back(t);
    }
  }

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
        std::string(PATH_INPUT_DIR) + "/smpl_model_f.npz");
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

  std::vector<double> aXYZ1 = aXYZ0;
  { // initalize pose
    dfm2::UpdateBoneRotTrans(aBone);
    dfm2::Skinning_LBS(aXYZ1,
                       aXYZ0, aBone, aW);
  }
  std::vector<std::pair<dfm2::CVec3d, dfm2::CVec3d> > aTargetOriginPos;
  for (auto &it : aTarget) {
    unsigned int ib = it.ib;
    aTargetOriginPos.emplace_back(it.pos,
                                  aBone[ib].RootPosition());
  }

  // -----------
  dfm2::glfw::CViewer3 viewer;
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  // viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::YTOP;
  dfm2::opengl::setSomeLighting();
  tex.InitGL();

  int iframe = 0;
  while (true) {
    {
      double r = iframe * 0.01;
      if (r > 1) { r = 1; }
      for (unsigned int it = 0; it < aTarget.size(); ++it) {
        aTarget[it].pos = r * aTargetOriginPos[it].first + (1 - r) * aTargetOriginPos[it].second;
      }
      Solve_MinRigging(aBone, aTarget);
      Skinning_LBS(aXYZ1,
                   aXYZ0, aBone, aW);
    }
    if (iframe > 200) {
      for (auto &ib : aBone) {
        dfm2::Quat_Identity(ib.quatRelativeRot);
      }
      dfm2::UpdateBoneRotTrans(aBone);
      for (unsigned int it = 0; it < aTarget.size(); ++it) {
        aTarget[it].pos = aTargetOriginPos[it].second;
      }
      Solve_MinRigging(aBone, aTarget);
      Skinning_LBS(aXYZ1,
                   aXYZ0, aBone, aW);
      iframe = 0;
    }

    // -------------------
    viewer.DrawBegin_oldGL();
    Draw(aXYZ1, aTri, aBone, aTarget);
    tex.Draw_oldGL();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    iframe++;
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
